#!/bin/bash
#
# 06_crossbuild_fingerprint.sh - Wrapper for the P22 cross-build sample-identity
# fingerprint: genotypes s_C_FLCCA3 (GRCh37 IMPACT BAM) and a patient's scRNA BAMs
# (GRCh38) at a shared lifted site panel, reverse-lifts s_C_FLCCA3's genotypes
# back to hg38, and runs the concordance analysis. See
# docs/superpowers/specs/2026-10-08-p22-crossbuild-fingerprint-design.md and
# 06_liftover_sites.sh (run first -- produces lifted_sites_vcf/original_sites_vcf).
#
# Usage:
#   ./06_crossbuild_fingerprint.sh <lifted_sites_vcf> <original_sites_vcf> \
#       <s_flcca3_bam> <b37_ref_fasta> <hg38_ref_fasta> <reverse_chain_file> \
#       <sample_manifest_csv> <patient_id> <out_dir>
#
# Required:
#   lifted_sites_vcf    Output of 06_liftover_sites.sh (sites.hg19_from_hg38.vcf.gz)
#   original_sites_vcf  06_liftover_sites.sh's stable-ID hg38 panel
#                        (sites.hg38.rna.stable_ids.vcf.gz) -- NOT the raw upstream
#                        panel, so P22's genotype VCFs carry the same IDs
#                        Comparison A joins on
#   s_flcca3_bam         s_C_FLCCA3_T001_d01.rg.md.abra.printreads.bam (GRCh37)
#   b37_ref_fasta        Homo_sapiens_assembly19.fasta
#   hg38_ref_fasta       refdata-gex-GRCh38-2020-A genome.fa
#   reverse_chain_file   UCSC hg19ToHg38.over.chain.gz (a separately-built UCSC
#                        chain, not a computed inverse of 06_liftover_sites.sh's chain)
#   sample_manifest_csv  04_fingerprint_manifest.R's manifest CSV (bam_path,sample_id,
#                        patient_id,... columns) -- filtered internally to patient_id
#   patient_id           Patient to filter the manifest to (e.g. P22)
#   out_dir              Output directory
#
# Exit Codes:
#   0  - Success
#   1  - Invalid arguments or missing files/commands
#   10 - bcftools mpileup/call failed for s_C_FLCCA3
#   11 - bcftools mpileup/call failed for one or more scRNA BAMs
#   12 - Reverse Picard LiftoverVcf failed
#   13 - Downstream R analysis failed
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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_SCRIPT="${SCRIPT_DIR}/06_crossbuild_fingerprint_analysis.R"

## ---- Validate file arguments (before activating the conda env, so a bad
## argument fails fast with a clear message instead of a conda activation
## error if the env happens to be unavailable) ---------------------------------
for required_file in "${LIFTED_SITES_VCF}" "${ORIGINAL_SITES_VCF}" "${S_FLCCA3_BAM}" \
                      "${B37_REF_FASTA}" "${HG38_REF_FASTA}" "${REVERSE_CHAIN_FILE}" \
                      "${SAMPLE_MANIFEST_CSV}" "${ANALYSIS_SCRIPT}"; do
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
    echo "$(date) - [Step 2/4] Genotyping s_C_FLCCA3 (b37) at the lifted site panel..."
    set +e
    bcftools mpileup -f "${B37_REF_FASTA}" -T "${LIFTED_SITES_VCF}" "${S_FLCCA3_BAM}" 2>"${OUT_DIR}/s_C_FLCCA3.mpileup.log" \
        | bcftools call -mv -Oz -o "${S_FLCCA3_RAW_VCF}"
    call_status=$?
    set -e
    if [ ${call_status} -ne 0 ] || [ ! -f "${S_FLCCA3_RAW_VCF}" ]; then
        echo "$(date) - [FAILED] bcftools mpileup/call failed for s_C_FLCCA3; see ${OUT_DIR}/s_C_FLCCA3.mpileup.log" >&2
        exit 10
    fi
    ## bcftools call does not carry the -T sites file's ID field through on its own.
    ## (LIFTED_SITES_VCF is the "-a" annotations file here and must be indexed --
    ## it already is, from 06_liftover_sites.sh's CREATE_INDEX=true; the target
    ## file being annotated does not need its own index.)
    bcftools annotate -a "${LIFTED_SITES_VCF}" -c ID "${S_FLCCA3_RAW_VCF}" -Oz -o "${S_FLCCA3_A_VCF}"
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

    if [ ! -f "${annotated_vcf}" ]; then
        set +e
        bcftools mpileup -f "${HG38_REF_FASTA}" -T "${ORIGINAL_SITES_VCF}" "${bam_path}" 2>"${OUT_DIR}/${sample_id}.mpileup.log" \
            | bcftools call -mv -Oz -o "${raw_vcf}"
        call_status=$?
        set -e
        if [ ${call_status} -ne 0 ] || [ ! -f "${raw_vcf}" ]; then
            echo "$(date) - [FAILED] bcftools mpileup/call failed for ${sample_id}; see ${OUT_DIR}/${sample_id}.mpileup.log" >&2
            touch "${CALL_FAILED}"
            continue
        fi
        bcftools annotate -a "${ORIGINAL_SITES_VCF}" -c ID "${raw_vcf}" -Oz -o "${annotated_vcf}"
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
    echo "$(date) - [Step 4/4] Reverse-lifting s_C_FLCCA3's genotypes (hg19 -> hg38)..."
    set +e
    picard LiftoverVcf \
        I="${S_FLCCA3_A_VCF}" \
        O="${S_FLCCA3_B_VCF}" \
        CHAIN="${REVERSE_CHAIN_FILE}" \
        REJECT="${S_FLCCA3_B_REJECTED}" \
        R="${HG38_REF_FASTA}" \
        WARN_ON_MISSING_CONTIG=true \
        CREATE_INDEX=true
    reverse_liftover_status=$?
    set -e
    if [ ${reverse_liftover_status} -ne 0 ] || [ ! -f "${S_FLCCA3_B_VCF}" ]; then
        echo "$(date) - [FAILED] Reverse Picard LiftoverVcf failed." >&2
        exit 12
    fi
else
    echo "$(date) - [Step 4/4] checkpoint found, skipping: ${S_FLCCA3_B_VCF}"
fi

## ---- Step 5: run the concordance analysis ------------------------------------
echo "$(date) - Running crossbuild fingerprint analysis..."
set +e
"${RSCRIPT}" "${ANALYSIS_SCRIPT}" \
    --s_flcca3_genotypes_a "${S_FLCCA3_A_VCF}" \
    --s_flcca3_genotypes_b "${S_FLCCA3_B_VCF}" \
    --p22_manifest "${PATIENT_GENOTYPE_MANIFEST}" \
    --out_dir "${OUT_DIR}"
analysis_status=$?
set -e

if [ ${analysis_status} -ne 0 ]; then
    echo "$(date) - [FAILED] Crossbuild fingerprint analysis failed." >&2
    exit 13
fi

echo "$(date) - [DONE] ${PATIENT_ID} cross-build fingerprint complete: ${OUT_DIR}"
