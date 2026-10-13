#!/bin/bash
#
# 06_liftover_sites.sh - LiftOver the RNA-optimized SNP site panel from GRCh38 to
# GRCh37/b37, for the P22 cross-build sample-identity fingerprint (design spec:
# docs/superpowers/specs/2026-10-08-p22-crossbuild-fingerprint-design.md in the
# z_palbo_antiPD1_11044 project checkout, not in this pipeline repo).
#
# Ensures every site has a stable ID (required so downstream genotype VCFs can be
# joined across builds by ID -- see Comparison A in the spec), then runs Picard
# LiftoverVcf hg38->hg19 and aborts if too many sites fail to lift (a high
# rejection rate usually means a chain/reference mismatch, not normal attrition).
#
# Contig naming (never changed on the data files themselves):
#   input sites_vcf          GRCh38, chr-prefixed (chr1..chr22, chrX, chrY)
#   UCSC hg38ToHg19 chain    chr-prefixed on BOTH sides (source chr1 -> dest chr1)
#   target_ref_fasta (b37)   bare names (1..22, X, Y, MT)
# So the chain's DESTINATION side (header field 8) is rewritten chrN -> N for
# 1..22/X/Y into a local copy (lib/chain_strip_chr.awk, col=8) before Picard
# runs; the source side stays chr-prefixed to match the input VCF. chrM is not
# renamed (hg19 chrM != b37 MT sequence) and is rejected like any unmapped contig.
#
# Usage:
#   ./06_liftover_sites.sh <sites_vcf> <chain_file> <target_ref_fasta> <out_dir> [max_reject_pct]
#
# Required:
#   sites_vcf         RNA-optimized GRCh38 site panel (sites.hg38.rna.vcf.gz)
#   chain_file        UCSC hg38ToHg19.over.chain.gz, as downloaded (see Task 1 of the
#                     implementation plan); its b37-side contig names are rewritten
#                     into out_dir/hg38ToHg19.b37_dest.chain, the original is untouched
#   target_ref_fasta  GRCh37/b37 reference FASTA matching the IMPACT BAM's alignment
#                     (e.g. Homo_sapiens_assembly19.fasta); this script creates a
#                     matching .dict if one doesn't already exist alongside it
#   out_dir           Output directory
#
# Optional:
#   max_reject_pct    Abort if Picard LiftoverVcf rejects more than this % of sites
#                     [default: 20]
#
# Exit Codes:
#   0  - Success
#   1  - Invalid arguments or missing files/commands
#   10 - ID-stabilization step failed
#   11 - Picard LiftoverVcf failed
#   12 - Rejection rate exceeded max_reject_pct
#
# Examples:
#   ./06_liftover_sites.sh \
#       /data1/singers/genomes/homo_sapiens/Ensembl/GRCh38.p5/Annotation/Variation/sites.hg38.rna.vcf.gz \
#       /data1/singers/genomes/homo_sapiens/Ensembl/GRCh38.p5/Annotation/Variation/liftover/hg38ToHg19.over.chain.gz \
#       /data1/core006/access/production/resources/reference/versions/hg19/Homo_sapiens_assembly19.fasta \
#       fingerprint/crossbuild_p22
#
set -euo pipefail

## ---- Parse arguments --------------------------------------------------------
SITES_VCF="${1:?Error: sites_vcf is required. Usage: $0 <sites_vcf> <chain_file> <target_ref_fasta> <out_dir> [max_reject_pct]}"
CHAIN_FILE="${2:?Error: chain_file is required. Usage: $0 <sites_vcf> <chain_file> <target_ref_fasta> <out_dir> [max_reject_pct]}"
TARGET_REF_FASTA="${3:?Error: target_ref_fasta is required. Usage: $0 <sites_vcf> <chain_file> <target_ref_fasta> <out_dir> [max_reject_pct]}"
OUT_DIR="${4:?Error: out_dir is required. Usage: $0 <sites_vcf> <chain_file> <target_ref_fasta> <out_dir> [max_reject_pct]}"
MAX_REJECT_PCT="${5:-20}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STABILIZE_AWK="${SCRIPT_DIR}/lib/stabilize_vcf_ids.awk"
CHAIN_STRIP_AWK="${SCRIPT_DIR}/lib/chain_strip_chr.awk"

## ---- Validate file arguments (before activating the conda env, so a bad
## argument fails fast with a clear message instead of a conda activation
## error if the env happens to be unavailable) ---------------------------------
for required_file in "${SITES_VCF}" "${CHAIN_FILE}" "${TARGET_REF_FASTA}" "${STABILIZE_AWK}" "${CHAIN_STRIP_AWK}"; do
    if [ ! -f "${required_file}" ]; then
        echo "Error: Required file not found: ${required_file}" >&2
        exit 1
    fi
done

## ---- Environment paths --------------------------------------------------------
source /usersoftware/singers/miniforge3/bin/activate sample-fingerprinting

## ---- Validate commands on PATH (needs the env activated first) --------------
for required_cmd in picard bgzip tabix; do
    if ! command -v "${required_cmd}" >/dev/null 2>&1; then
        echo "Error: ${required_cmd} not found on PATH after activating sample-fingerprinting" >&2
        exit 1
    fi
done

mkdir -p "${OUT_DIR}"

## ---- Top-level checkpoint ----------------------------------------------------
## Checkpoint files are only ever created by an mv AFTER their step's checks
## pass (outputs are written to *.tmp.vcf.gz first), so a file existing here
## means it passed validation -- an aborted run cannot leave a "done" file behind.
LIFTED_VCF="${OUT_DIR}/sites.hg19_from_hg38.vcf.gz"
STABLE_SITES_VCF="${OUT_DIR}/sites.hg38.rna.stable_ids.vcf.gz"
if [ -f "${LIFTED_VCF}" ] && [ -f "${STABLE_SITES_VCF}" ]; then
    echo "$(date) - [DONE] Checkpoint found, liftover already complete: ${LIFTED_VCF}"
    exit 0
fi

## ---- Step 1: ensure stable site IDs (required for Comparison A's join) ------
if [ ! -f "${STABLE_SITES_VCF}" ]; then
    echo "$(date) - [Step 1/3] Stabilizing site IDs..."
    STABLE_SITES_PLAIN="${OUT_DIR}/sites.hg38.rna.stable_ids.vcf"

    set +e
    gzip -dc "${SITES_VCF}" | awk -f "${STABILIZE_AWK}" > "${STABLE_SITES_PLAIN}"
    stabilize_status=$?
    set -e

    if [ ${stabilize_status} -ne 0 ] || [ ! -s "${STABLE_SITES_PLAIN}" ]; then
        echo "$(date) - [FAILED] ID stabilization failed." >&2
        exit 10
    fi

    STABLE_SITES_TMP="${OUT_DIR}/sites.hg38.rna.stable_ids.tmp.vcf.gz"
    bgzip -c "${STABLE_SITES_PLAIN}" > "${STABLE_SITES_TMP}"
    tabix -f -p vcf "${STABLE_SITES_TMP}"
    rm -f "${STABLE_SITES_PLAIN}"
    mv "${STABLE_SITES_TMP}.tbi" "${STABLE_SITES_VCF}.tbi"
    mv "${STABLE_SITES_TMP}" "${STABLE_SITES_VCF}"
    echo "$(date) - [Step 1/3] Stable-ID site VCF: ${STABLE_SITES_VCF}"
else
    echo "$(date) - [Step 1/3] checkpoint found, skipping: ${STABLE_SITES_VCF}"
fi

## ---- Step 2: build target .dict if missing -----------------------------------
TARGET_DICT="${TARGET_REF_FASTA%.*}.dict"
if [ ! -f "${TARGET_DICT}" ]; then
    echo "$(date) - [Step 2/3] Building sequence dictionary for ${TARGET_REF_FASTA}..."
    picard CreateSequenceDictionary R="${TARGET_REF_FASTA}" O="${TARGET_DICT}"
else
    echo "$(date) - [Step 2/3] Sequence dictionary already exists: ${TARGET_DICT}"
fi

## ---- Step 3: Picard LiftoverVcf hg38 -> hg19 ---------------------------------
REJECTED_VCF="${OUT_DIR}/sites.hg19_from_hg38.rejected.vcf.gz"
LIFTED_TMP="${OUT_DIR}/sites.hg19_from_hg38.tmp.vcf.gz"
B37_CHAIN="${OUT_DIR}/hg38ToHg19.b37_dest.chain"

## UCSC chain dest side is chr-prefixed; b37 .dict is not -- see header.
## gzip -dcf passes an already-uncompressed chain through unchanged.
gzip -dcf "${CHAIN_FILE}" | awk -v col=8 -f "${CHAIN_STRIP_AWK}" > "${B37_CHAIN}"

echo "$(date) - [Step 3/3] Running Picard LiftoverVcf (hg38 -> hg19)..."
rm -f "${LIFTED_TMP}" "${LIFTED_TMP}.tbi"
set +e
picard LiftoverVcf \
    I="${STABLE_SITES_VCF}" \
    O="${LIFTED_TMP}" \
    CHAIN="${B37_CHAIN}" \
    REJECT="${REJECTED_VCF}" \
    R="${TARGET_REF_FASTA}" \
    WARN_ON_MISSING_CONTIG=true \
    CREATE_INDEX=true
liftover_status=$?
set -e

if [ ${liftover_status} -ne 0 ] || [ ! -f "${LIFTED_TMP}" ]; then
    echo "$(date) - [FAILED] Picard LiftoverVcf failed." >&2
    rm -f "${LIFTED_TMP}" "${LIFTED_TMP}.tbi"
    exit 11
fi

## gzip -dc (not zcat) -- BSD zcat (macOS) refuses .gz input without a .Z
## extension, while gzip -dc decompresses .gz identically on GNU and BSD.
## grep -vc exits 1 when the count is 0, which would otherwise trip `set -e`
## on a clean "zero rejects" run -- `|| true` absorbs that without masking a
## real gzip/grep failure, since gzip -dc's own exit status isn't what `||`
## is guarding here.
n_input=$(gzip -dc "${STABLE_SITES_VCF}" | { grep -vc '^#' || true; })
n_rejected=$(gzip -dc "${REJECTED_VCF}" 2>/dev/null | { grep -vc '^#' || true; })
reject_pct=$(awk -v r="${n_rejected}" -v n="${n_input}" 'BEGIN { printf "%.2f", (n > 0) ? (100 * r / n) : 0 }')

echo "$(date) - [Step 3/3] Liftover complete: ${n_rejected}/${n_input} sites rejected (${reject_pct}%)"

if awk -v p="${reject_pct}" -v max="${MAX_REJECT_PCT}" 'BEGIN { exit !(p > max) }'; then
    echo "Error: rejection rate ${reject_pct}% exceeds max_reject_pct=${MAX_REJECT_PCT}%. This usually means a chain/reference mismatch, not normal liftover attrition -- check CHAIN_FILE and TARGET_REF_FASTA before trusting this panel. See: ${REJECTED_VCF}" >&2
    rm -f "${LIFTED_TMP}" "${LIFTED_TMP}.tbi"
    exit 12
fi

## All checks passed -- only now promote to the checkpoint name (index first,
## so the checkpoint file never exists without its index).
mv "${LIFTED_TMP}.tbi" "${LIFTED_VCF}.tbi"
mv "${LIFTED_TMP}" "${LIFTED_VCF}"

echo "$(date) - [DONE] Lifted site panel: ${LIFTED_VCF} (rejected sites logged: ${REJECTED_VCF})"
