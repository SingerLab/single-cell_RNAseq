#!/bin/bash
#
# 06_crossbuild_fingerprint.sh - Wrapper for the P22 cross-build sample-identity
# fingerprint: genotypes s_C_FLCCA3 (GRCh37 IMPACT BAM) and a patient's scRNA BAMs
# (GRCh38) at a shared lifted site panel, reverse-lifts s_C_FLCCA3's genotypes
# back to hg38, and runs the concordance analysis. See the design spec
# (docs/superpowers/specs/2026-10-08-p22-crossbuild-fingerprint-design.md in the
# z_palbo_antiPD1_11044 project checkout) and 06_liftover_sites.sh (run first --
# produces lifted_sites_vcf/original_sites_vcf).
#
# Genotyping uses `bcftools call -m` (NOT the spec diagram's `-mv`) plus a
# FORMAT/DP >= min_dp filter: -v would drop hom-ref calls, so the inner join
# would only see sites where BOTH samples are variant and hom-ref-vs-variant
# disagreements would vanish. This honors the spec's "adequate-depth
# informative sites" intent rather than its diagram's literal flags.
#
# Contig naming (never changed on the data files themselves):
#   s_C_FLCCA3 BAM / its genotype VCF   b37, bare names (1..22, X, Y, MT)
#   UCSC hg19ToHg38 chain               chr-prefixed on BOTH sides
#   hg38_ref_fasta (CellRanger 2020-A)  chr-prefixed (chr1..chr22, chrX, chrY, chrM)
#                                       -- per the spec's Inputs table; same as P22's BAMs
# So the reverse chain's SOURCE side (header field 3) is rewritten chrN -> N for
# 1..22/X/Y into a local copy (lib/chain_strip_chr.awk, col=3); its destination
# side stays chr-prefixed to match hg38_ref_fasta. MT is not mapped (hg19 chrM
# != b37 MT sequence) and is rejected like any unmapped contig.
#
# Usage:
#   ./06_crossbuild_fingerprint.sh <lifted_sites_vcf> <original_sites_vcf> \
#       <s_flcca3_bam> <b37_ref_fasta> <hg38_ref_fasta> <reverse_chain_file> \
#       <sample_manifest_csv> <patient_id> <out_dir> [min_dp] [max_reject_pct]
#
# Required:
#   lifted_sites_vcf    Output of 06_liftover_sites.sh (sites.hg19_from_hg38.vcf.gz)
#   original_sites_vcf  06_liftover_sites.sh's stable-ID hg38 panel
#                        (sites.hg38.rna.stable_ids.vcf.gz) -- NOT the raw upstream
#                        panel, so P22's genotype VCFs carry the same IDs
#                        Comparison A joins on
#   s_flcca3_bam         s_C_FLCCA3_T001_d01.rg.md.abra.printreads.bam (GRCh37)
#   b37_ref_fasta        Homo_sapiens_assembly19.fasta
#   hg38_ref_fasta       refdata-gex-GRCh38-2020-A genome.fa; this script creates a
#                        matching .dict if one doesn't already exist alongside it
#   reverse_chain_file   UCSC hg19ToHg38.over.chain.gz, as downloaded (a separately-built
#                        UCSC chain, not a computed inverse of 06_liftover_sites.sh's
#                        chain); its b37-side names are rewritten into
#                        out_dir/hg19ToHg38.b37_source.chain, the original is untouched
#   sample_manifest_csv  04_fingerprint_manifest.R's manifest CSV (bam_path,sample_id,
#                        patient_id,... columns) -- filtered internally to patient_id
#   patient_id           Patient to filter the manifest to (e.g. P22)
#   out_dir              Output directory
#
# Optional:
#   min_dp               Minimum per-sample FORMAT/DP for a genotype call to be kept,
#                        applied to s_C_FLCCA3 and every scRNA BAM [default: 10]
#   max_reject_pct       Abort if the reverse (hg19 -> hg38) Picard LiftoverVcf rejects
#                        more than this % of s_C_FLCCA3's calls [default: 20]
#
# Exit Codes:
#   0  - Success
#   1  - Invalid arguments or missing files/commands
#   10 - bcftools mpileup/call failed for s_C_FLCCA3
#   11 - bcftools mpileup/call failed for one or more scRNA BAMs
#   12 - Reverse Picard LiftoverVcf failed
#   13 - Downstream R analysis failed
#   14 - Reverse liftover rejection rate exceeded max_reject_pct
#
# Examples:
#   ./06_crossbuild_fingerprint.sh fingerprint/crossbuild_p22/sites.hg19_from_hg38.vcf.gz \
#       fingerprint/crossbuild_p22/sites.hg38.rna.stable_ids.vcf.gz \
#       /data1/singers/share/sandra_dangelo_lab/ResearchIMPACT_20062/17549/1.7.0/20260204_21_11_791149/bam/s_C_FLCCA3_T001_d01.rg.md.abra.printreads.bam \
#       /data1/core006/access/production/resources/reference/versions/hg19/Homo_sapiens_assembly19.fasta \
#       /data1/singers/genomes/homo_sapiens/Ensembl/GRCh38.p5/Sequence/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
#       /data1/singers/genomes/homo_sapiens/Ensembl/GRCh38.p5/Annotation/Variation/liftover/hg19ToHg38.over.chain.gz \
#       fingerprint/sample_manifest.csv P22 fingerprint/crossbuild_p22
#
set -euo pipefail

## ---- Parse arguments --------------------------------------------------------
LIFTED_SITES_VCF="${1:?Error: lifted_sites_vcf is required. See the header comment above for full usage.}"
ORIGINAL_SITES_VCF="${2:?Error: original_sites_vcf is required.}"
S_FLCCA3_BAM="${3:?Error: s_flcca3_bam is required.}"
B37_REF_FASTA="${4:?Error: b37_ref_fasta is required.}"
HG38_REF_FASTA="${5:?Error: hg38_ref_fasta is required.}"
REVERSE_CHAIN_FILE="${6:?Error: reverse_chain_file is required.}"
SAMPLE_MANIFEST_CSV="${7:?Error: sample_manifest_csv is required.}"
PATIENT_ID="${8:?Error: patient_id is required.}"
OUT_DIR="${9:?Error: out_dir is required.}"
MIN_DP="${10:-10}"
MAX_REJECT_PCT="${11:-20}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_SCRIPT="${SCRIPT_DIR}/06_crossbuild_fingerprint_analysis.R"
CHAIN_STRIP_AWK="${SCRIPT_DIR}/lib/chain_strip_chr.awk"

## ---- Validate file arguments (before activating the conda env, so a bad
## argument fails fast with a clear message instead of a conda activation
## error if the env happens to be unavailable) ---------------------------------
for required_file in "${LIFTED_SITES_VCF}" "${ORIGINAL_SITES_VCF}" "${S_FLCCA3_BAM}" \
                      "${B37_REF_FASTA}" "${HG38_REF_FASTA}" "${REVERSE_CHAIN_FILE}" \
                      "${SAMPLE_MANIFEST_CSV}" "${ANALYSIS_SCRIPT}" "${CHAIN_STRIP_AWK}"; do
    if [ ! -f "${required_file}" ]; then
        echo "Error: Required file not found: ${required_file}" >&2
        exit 1
    fi
done

## ---- Environment paths --------------------------------------------------------
source /usersoftware/singers/miniforge3/bin/activate sample-fingerprinting
RSCRIPT="/usersoftware/singers/miniforge3/envs/sample-fingerprinting/bin/Rscript"

## ---- Validate commands on PATH (needs the env activated first) --------------
for required_cmd in bcftools picard; do
    if ! command -v "${required_cmd}" >/dev/null 2>&1; then
        echo "Error: ${required_cmd} not found on PATH after activating sample-fingerprinting" >&2
        exit 1
    fi
done

mkdir -p "${OUT_DIR}"

## Non-header record count of a (b)gzipped VCF. gzip -dc (not zcat) -- BSD zcat
## (macOS) refuses .gz without a .Z extension; grep -vc exits 1 on a zero count,
## which `|| true` absorbs (same pattern as 06_liftover_sites.sh).
count_vcf_records() {
    gzip -dc "$1" 2>/dev/null | { grep -vc '^#' || true; }
}

## Checkpoint files below are only ever created by an mv AFTER their step's
## checks pass (outputs go to *.tmp.vcf.gz first), so an aborted run can't leave
## a half-written file that a rerun would skip over as done.

## ---- Step 1: derive this patient's BAM manifest from the cohort manifest -----
PATIENT_BAM_MANIFEST="${OUT_DIR}/${PATIENT_ID}_bam_manifest.csv"
echo "sample_id,bam_path" > "${PATIENT_BAM_MANIFEST}"
tail -n +2 "${SAMPLE_MANIFEST_CSV}" | awk -F',' -v pid="${PATIENT_ID}" '{gsub(/"/,""); if ($3==pid) print $2","$1}' >> "${PATIENT_BAM_MANIFEST}"

n_patient_bams=$(tail -n +2 "${PATIENT_BAM_MANIFEST}" | wc -l)
if [ "${n_patient_bams}" -eq 0 ]; then
    echo "Error: no rows for patient_id=${PATIENT_ID} found in ${SAMPLE_MANIFEST_CSV}" >&2
    exit 1
fi
echo "$(date) - [Step 1/4] ${n_patient_bams} BAM(s) for ${PATIENT_ID}: ${PATIENT_BAM_MANIFEST}"

## ---- Step 2: genotype s_C_FLCCA3 (native b37) at the lifted site panel -------
S_FLCCA3_A_VCF="${OUT_DIR}/s_C_FLCCA3.genotypes.b37.vcf.gz"

if [ ! -f "${S_FLCCA3_A_VCF}" ]; then
    S_FLCCA3_RAW_VCF="${OUT_DIR}/s_C_FLCCA3.genotypes.b37.raw.vcf.gz"
    S_FLCCA3_A_TMP="${OUT_DIR}/s_C_FLCCA3.genotypes.b37.tmp.vcf.gz"
    echo "$(date) - [Step 2/4] Genotyping s_C_FLCCA3 (b37) at the lifted site panel (FORMAT/DP>=${MIN_DP})..."
    ## bcftools call does not carry the -T sites file's ID field through on its own,
    ## hence the annotate. (LIFTED_SITES_VCF is the "-a" annotations file here and
    ## must be indexed -- it already is, from 06_liftover_sites.sh's CREATE_INDEX=true;
    ## the target file being annotated does not need its own index.)
    set +e
    bcftools mpileup -a FORMAT/DP -f "${B37_REF_FASTA}" -T "${LIFTED_SITES_VCF}" "${S_FLCCA3_BAM}" 2>"${OUT_DIR}/s_C_FLCCA3.mpileup.log" \
        | bcftools call -m -Ou \
        | bcftools view -i "FORMAT/DP>=${MIN_DP}" -Oz -o "${S_FLCCA3_RAW_VCF}" \
        && bcftools annotate -a "${LIFTED_SITES_VCF}" -c ID "${S_FLCCA3_RAW_VCF}" -Oz -o "${S_FLCCA3_A_TMP}"
    call_status=$?
    set -e
    if [ ${call_status} -ne 0 ] || [ ! -f "${S_FLCCA3_A_TMP}" ]; then
        echo "$(date) - [FAILED] bcftools mpileup/call/annotate failed for s_C_FLCCA3; see ${OUT_DIR}/s_C_FLCCA3.mpileup.log" >&2
        rm -f "${S_FLCCA3_A_TMP}"
        exit 10
    fi
    mv "${S_FLCCA3_A_TMP}" "${S_FLCCA3_A_VCF}"
    echo "$(date) - [Step 2/4] s_C_FLCCA3 genotyped (Comparison A input): ${S_FLCCA3_A_VCF}"
else
    echo "$(date) - [Step 2/4] checkpoint found, skipping: ${S_FLCCA3_A_VCF}"
fi

## ---- Step 3: genotype each scRNA BAM (native hg38) at the original panel ----
PATIENT_GENOTYPE_MANIFEST="${OUT_DIR}/${PATIENT_ID}_genotype_manifest.csv"
echo "sample_id,genotype_vcf_path" > "${PATIENT_GENOTYPE_MANIFEST}"

echo "$(date) - [Step 3/4] Genotyping ${PATIENT_ID}'s scRNA BAMs (hg38) at the original site panel..."
CALL_FAILED="${OUT_DIR}/.${PATIENT_ID}_call_failed"
rm -f "${CALL_FAILED}"

tail -n +2 "${PATIENT_BAM_MANIFEST}" | while IFS=',' read -r sample_id bam_path; do
    raw_vcf="${OUT_DIR}/${sample_id}.genotypes.hg38.raw.vcf.gz"
    annotated_vcf="${OUT_DIR}/${sample_id}.genotypes.hg38.vcf.gz"
    annotated_tmp="${OUT_DIR}/${sample_id}.genotypes.hg38.tmp.vcf.gz"

    if [ ! -f "${annotated_vcf}" ]; then
        set +e
        bcftools mpileup -a FORMAT/DP -f "${HG38_REF_FASTA}" -T "${ORIGINAL_SITES_VCF}" "${bam_path}" 2>"${OUT_DIR}/${sample_id}.mpileup.log" \
            | bcftools call -m -Ou \
            | bcftools view -i "FORMAT/DP>=${MIN_DP}" -Oz -o "${raw_vcf}" \
            && bcftools annotate -a "${ORIGINAL_SITES_VCF}" -c ID "${raw_vcf}" -Oz -o "${annotated_tmp}"
        call_status=$?
        set -e
        if [ ${call_status} -ne 0 ] || [ ! -f "${annotated_tmp}" ]; then
            echo "$(date) - [FAILED] bcftools mpileup/call/annotate failed for ${sample_id}; see ${OUT_DIR}/${sample_id}.mpileup.log" >&2
            rm -f "${annotated_tmp}"
            touch "${CALL_FAILED}"
            continue
        fi
        mv "${annotated_tmp}" "${annotated_vcf}"
    fi
    echo "${sample_id},${annotated_vcf}" >> "${PATIENT_GENOTYPE_MANIFEST}"
done

if [ -f "${CALL_FAILED}" ]; then
    echo "$(date) - [FAILED] One or more ${PATIENT_ID} BAM genotype calls failed; see per-sample .mpileup.log files." >&2
    exit 11
fi
echo "$(date) - [Step 3/4] Genotyping complete: ${PATIENT_GENOTYPE_MANIFEST}"

## ---- Step 4: reverse-lift s_C_FLCCA3's genotypes hg19 -> hg38 (Comparison B) -
S_FLCCA3_B_VCF="${OUT_DIR}/s_C_FLCCA3.genotypes.liftedToHg38.vcf.gz"
S_FLCCA3_B_REJECTED="${OUT_DIR}/s_C_FLCCA3.genotypes.liftedToHg38.rejected.vcf.gz"

if [ ! -f "${S_FLCCA3_B_VCF}" ]; then
    S_FLCCA3_B_TMP="${OUT_DIR}/s_C_FLCCA3.genotypes.liftedToHg38.tmp.vcf.gz"
    REVERSE_CHAIN_B37="${OUT_DIR}/hg19ToHg38.b37_source.chain"

    ## CellRanger reference bundles usually ship only .fa + .fai; Picard needs a .dict.
    HG38_DICT="${HG38_REF_FASTA%.*}.dict"
    if [ ! -f "${HG38_DICT}" ]; then
        echo "$(date) - [Step 4/4] Building sequence dictionary for ${HG38_REF_FASTA}..."
        picard CreateSequenceDictionary R="${HG38_REF_FASTA}" O="${HG38_DICT}"
    fi

    ## UCSC chain source side is chr-prefixed; s_C_FLCCA3's b37 calls are not -- see header.
    ## gzip -dcf passes an already-uncompressed chain through unchanged.
    gzip -dcf "${REVERSE_CHAIN_FILE}" | awk -v col=3 -f "${CHAIN_STRIP_AWK}" > "${REVERSE_CHAIN_B37}"

    echo "$(date) - [Step 4/4] Reverse-lifting s_C_FLCCA3's genotypes (hg19 -> hg38)..."
    rm -f "${S_FLCCA3_B_TMP}" "${S_FLCCA3_B_TMP}.tbi"
    set +e
    picard LiftoverVcf \
        I="${S_FLCCA3_A_VCF}" \
        O="${S_FLCCA3_B_TMP}" \
        CHAIN="${REVERSE_CHAIN_B37}" \
        REJECT="${S_FLCCA3_B_REJECTED}" \
        R="${HG38_REF_FASTA}" \
        WARN_ON_MISSING_CONTIG=true \
        CREATE_INDEX=true
    reverse_liftover_status=$?
    set -e
    if [ ${reverse_liftover_status} -ne 0 ] || [ ! -f "${S_FLCCA3_B_TMP}" ]; then
        echo "$(date) - [FAILED] Reverse Picard LiftoverVcf failed." >&2
        rm -f "${S_FLCCA3_B_TMP}" "${S_FLCCA3_B_TMP}.tbi"
        exit 12
    fi

    ## Same >max_reject_pct abort as 06_liftover_sites.sh's forward direction.
    n_input=$(count_vcf_records "${S_FLCCA3_A_VCF}")
    n_rejected=$(count_vcf_records "${S_FLCCA3_B_REJECTED}")
    reject_pct=$(awk -v r="${n_rejected}" -v n="${n_input}" 'BEGIN { printf "%.2f", (n > 0) ? (100 * r / n) : 0 }')
    echo "$(date) - [Step 4/4] Reverse liftover complete: ${n_rejected}/${n_input} calls rejected (${reject_pct}%)"

    if awk -v p="${reject_pct}" -v max="${MAX_REJECT_PCT}" 'BEGIN { exit !(p > max) }'; then
        echo "Error: reverse liftover rejection rate ${reject_pct}% exceeds max_reject_pct=${MAX_REJECT_PCT}%. This usually means a chain/reference mismatch, not normal liftover attrition -- check REVERSE_CHAIN_FILE and HG38_REF_FASTA. See: ${S_FLCCA3_B_REJECTED}" >&2
        rm -f "${S_FLCCA3_B_TMP}" "${S_FLCCA3_B_TMP}.tbi"
        exit 14
    fi

    ## All checks passed -- only now promote to the checkpoint name (index first).
    mv "${S_FLCCA3_B_TMP}.tbi" "${S_FLCCA3_B_VCF}.tbi"
    mv "${S_FLCCA3_B_TMP}" "${S_FLCCA3_B_VCF}"
else
    echo "$(date) - [Step 4/4] checkpoint found, skipping: ${S_FLCCA3_B_VCF}"
fi

## ---- Step 5: run the concordance analysis ------------------------------------
## Rejection counts for both liftover directions, for the report (recomputed
## here so they're right on checkpointed reruns too). The forward reject VCF
## sits next to lifted_sites_vcf, per 06_liftover_sites.sh's naming.
FORWARD_REJECTED_VCF="${LIFTED_SITES_VCF%.vcf.gz}.rejected.vcf.gz"
forward_rejection="NA"
if [ -f "${FORWARD_REJECTED_VCF}" ]; then
    forward_rejection="$(count_vcf_records "${FORWARD_REJECTED_VCF}")/$(count_vcf_records "${ORIGINAL_SITES_VCF}")"
fi
reverse_rejection="$(count_vcf_records "${S_FLCCA3_B_REJECTED}")/$(count_vcf_records "${S_FLCCA3_A_VCF}")"

echo "$(date) - Running crossbuild fingerprint analysis..."
set +e
"${RSCRIPT}" "${ANALYSIS_SCRIPT}" \
    --s_flcca3_genotypes_a "${S_FLCCA3_A_VCF}" \
    --s_flcca3_genotypes_b "${S_FLCCA3_B_VCF}" \
    --p22_manifest "${PATIENT_GENOTYPE_MANIFEST}" \
    --out_dir "${OUT_DIR}" \
    --forward_liftover_rejected "${forward_rejection}" \
    --reverse_liftover_rejected "${reverse_rejection}"
analysis_status=$?
set -e

if [ ${analysis_status} -ne 0 ]; then
    echo "$(date) - [FAILED] Crossbuild fingerprint analysis failed." >&2
    exit 13
fi

echo "$(date) - [DONE] ${PATIENT_ID} cross-build fingerprint complete: ${OUT_DIR}"
