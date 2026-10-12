#!/bin/bash
#
# 06_liftover_sites.sh - LiftOver the RNA-optimized SNP site panel from GRCh38 to
# GRCh37/b37, for the P22 cross-build sample-identity fingerprint (see
# docs/superpowers/specs/2026-10-08-p22-crossbuild-fingerprint-design.md).
#
# Ensures every site has a stable ID (required so downstream genotype VCFs can be
# joined across builds by ID -- see Comparison A in the spec), then runs Picard
# LiftoverVcf hg38->hg19 and aborts if too many sites fail to lift (a high
# rejection rate usually means a chain/reference mismatch, not normal attrition).
#
# Usage:
#   ./06_liftover_sites.sh <sites_vcf> <chain_file> <target_ref_fasta> <out_dir> [max_reject_pct]
#
# Required:
#   sites_vcf         RNA-optimized GRCh38 site panel (sites.hg38.rna.vcf.gz)
#   chain_file        UCSC hg38ToHg19.over.chain.gz (see Task 1 of the implementation plan)
#   target_ref_fasta  GRCh37/b37 reference FASTA matching the IMPACT BAM's alignment
#                     (e.g. Homo_sapiens_assembly19.fasta); this script creates a
#                     matching .dict if one doesn't already exist alongside it
#   out_dir           Output directory
#
# Optional:
#   max_reject_pct    Abort if Picard LiftoverVcf rejects more than this %% of sites
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

## ---- Validate file arguments (before activating the conda env, so a bad
## argument fails fast with a clear message instead of a conda activation
## error if the env happens to be unavailable) ---------------------------------
for required_file in "${SITES_VCF}" "${CHAIN_FILE}" "${TARGET_REF_FASTA}" "${STABILIZE_AWK}"; do
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
LIFTED_VCF="${OUT_DIR}/sites.hg19_from_hg38.vcf.gz"
if [ -f "${LIFTED_VCF}" ]; then
    echo "$(date) - [DONE] Checkpoint found, liftover already complete: ${LIFTED_VCF}"
    exit 0
fi

## ---- Step 1: ensure stable site IDs (required for Comparison A's join) ------
STABLE_SITES_VCF="${OUT_DIR}/sites.hg38.rna.stable_ids.vcf.gz"

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

    bgzip -f "${STABLE_SITES_PLAIN}"
    tabix -p vcf "${STABLE_SITES_VCF}"
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

echo "$(date) - [Step 3/3] Running Picard LiftoverVcf (hg38 -> hg19)..."
set +e
picard LiftoverVcf \
    I="${STABLE_SITES_VCF}" \
    O="${LIFTED_VCF}" \
    CHAIN="${CHAIN_FILE}" \
    REJECT="${REJECTED_VCF}" \
    R="${TARGET_REF_FASTA}" \
    WARN_ON_MISSING_CONTIG=true \
    CREATE_INDEX=true
liftover_status=$?
set -e

if [ ${liftover_status} -ne 0 ] || [ ! -f "${LIFTED_VCF}" ]; then
    echo "$(date) - [FAILED] Picard LiftoverVcf failed." >&2
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
    exit 12
fi

echo "$(date) - [DONE] Lifted site panel: ${LIFTED_VCF} (rejected sites logged: ${REJECTED_VCF})"
