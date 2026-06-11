#!/bin/bash
#
# 03_numbat_cnv.sh - Wrapper for Numbat CNV analysis with pileup and phasing
#
# Usage:
#   ./03_numbat_cnv.sh <sample_dir> [ref_rds] [min_LLR] [ncores]
#
# Required:
#   sample_dir   Path to the 10X sample directory
#
# Optional:
#   ref_rds      Path to reference expression RDS file [default: numbat ref_hca]
#   min_LLR      Minimum log-likelihood ratio [default: 3]
#   ncores       Number of cores [default: 22]
#
# Exit Codes:
#   0  - Success
#   1  - Invalid arguments or missing files
#   10 - Pileup and phasing failed
#   11 - Pileup completed but output files are incomplete
#   20 - Numbat CNV analysis failed
#
# Examples:
#   ./03_numbat_cnv.sh P14_C0_D0
#   ./03_numbat_cnv.sh P14_C0_D0 refs/ref_sarc.rds
#   for dir in P*/; do ./03_numbat_cnv.sh "$dir"; done
#
set -euo pipefail

## ---- Parse arguments --------------------------------------------------------
SAMPLE_DIR="${1:?Error: sample_dir is required. Usage: $0 <sample_dir> [ref_rds] [min_LLR] [ncores]}"
REF_RDS="${2:-}"
MIN_LLR="${3:-3}"
NCORES="${4:-${LSB_MAX_NUM_PROCESSORS:-1}}"

## ---- Remove trailing slash --------------------------------------------------
SAMPLE_DIR="${SAMPLE_DIR%/}"
SAMPLE_NAME=$(basename "${SAMPLE_DIR}")

## ---- Environment paths ------------------------------------------------------
RSCRIPT="/work/singer/opt/miniconda3/envs/single-cell-rnaseq-v5/bin/Rscript"
EAGLE_DIR="/work/singer/opt/src/Eagle_v2.4.1"
GENOME="/work/singer/genomes/homo_sapiens/Ensembl/GRCh38.p5"
GENOME1K_VCF="${GENOME}/Annotation/Variation/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf"
PANEL_DIR="${GENOME}/Annotation/Variation/1000G_hg38"
PILEUP_SCRIPT="${HOME}/opt/numbat/inst/bin/pileup_and_phase.R"

## ---- Validate inputs --------------------------------------------------------
BAM_FILE="${SAMPLE_DIR}/outs/possorted_genome_bam.bam"
BARCODES="${SAMPLE_DIR}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"

for required_file in "${BAM_FILE}" "${BARCODES}" "${PILEUP_SCRIPT}" "${GENOME1K_VCF}"; do
    if [ ! -f "${required_file}" ]; then
        echo "Error: Required file not found: ${required_file}" >&2
        exit 1
    fi
done

## ---- Build Rscript arguments for numbat ------------------------------------
RSCRIPT_ARGS="--sample.dir ${SAMPLE_DIR} --min_LLR ${MIN_LLR} --ncores ${NCORES}"
if [ -n "${REF_RDS}" ]; then
    RSCRIPT_ARGS="${RSCRIPT_ARGS} --ref ${REF_RDS}"
fi

## ---- Create output directory ------------------------------------------------
NUMBAT_DIR="${SAMPLE_DIR}/numbat"
mkdir -p "${NUMBAT_DIR}"

## ---- Step 1: Pileup and Phase (with checkpoint) ----------------------------
ALLELE_FILE="${NUMBAT_DIR}/${SAMPLE_NAME}_allele_counts.tsv.gz"
PHASING_LOG="${NUMBAT_DIR}/phasing.log"
PILEUP_LOG="${NUMBAT_DIR}/pileup.log"

## Check if pileup has already completed successfully
pileup_complete=false

if [ -f "${ALLELE_FILE}" ] && [ -f "${PHASING_LOG}" ] && [ -f "${PILEUP_LOG}" ]; then
    ## Verify log files indicate successful completion
    phasing_ok=false
    pileup_ok=false

    if grep -q "PHASE_CONFIDENCE" "${PHASING_LOG}"; then
        phasing_ok=true
    fi

    if grep -q "All Done!" "${PILEUP_LOG}"; then
        pileup_ok=true
    fi

    if ${phasing_ok} && ${pileup_ok}; then
        pileup_complete=true
        echo "$(date) - [Step 1/2] Pileup checkpoint found for: ${SAMPLE_NAME}"
        echo "  Allele file: ${ALLELE_FILE}"
        echo "  Phasing log: OK"
        echo "  Pileup log:  OK"
        echo "  Skipping to Numbat CNV analysis..."
    else
        echo "$(date) - [Step 1/2] Incomplete pileup detected for: ${SAMPLE_NAME}"
        [ "${phasing_ok}" = false ] && echo "  Phasing log: INCOMPLETE"
        [ "${pileup_ok}" = false ]  && echo "  Pileup log:  INCOMPLETE"
        echo "  Re-running pileup and phasing..."
    fi
fi

## Run pileup if not already complete
if ! ${pileup_complete}; then
    echo "$(date) - [Step 1/2] Starting pileup_and_phase for: ${SAMPLE_NAME}"

    ${RSCRIPT} "${PILEUP_SCRIPT}" \
        --label "${SAMPLE_NAME}" \
        --samples "${SAMPLE_DIR}" \
        --bams "${BAM_FILE}" \
        --barcodes <(zcat "${BARCODES}") \
        --gmap "${EAGLE_DIR}/tables/genetic_map_hg38_withX.txt.gz" \
        --eagle "${EAGLE_DIR}/eagle" \
        --outdir "${NUMBAT_DIR}" \
        --snpvcf "${GENOME1K_VCF}" \
        --paneldir "${PANEL_DIR}" \
        --ncores "${NCORES}"

    if [ $? -ne 0 ]; then
        echo "$(date) - [FAILED] Pileup and phasing failed for: ${SAMPLE_NAME}" >&2
        exit 10
    fi

    ## Verify outputs were generated
    if [ ! -f "${ALLELE_FILE}" ]; then
        echo "$(date) - [FAILED] Allele counts file not generated: ${ALLELE_FILE}" >&2
        exit 11
    fi

    if [ ! -f "${PHASING_LOG}" ] || ! grep -q "PHASE_CONFIDENCE" "${PHASING_LOG}"; then
        echo "$(date) - [FAILED] Phasing did not complete successfully for: ${SAMPLE_NAME}" >&2
        exit 11
    fi

    if [ ! -f "${PILEUP_LOG}" ] || ! grep -q "All Done!" "${PILEUP_LOG}"; then
        echo "$(date) - [FAILED] Pileup did not complete successfully for: ${SAMPLE_NAME}" >&2
        exit 11
    fi

    echo "$(date) - [Step 1/2] Pileup and phasing completed for: ${SAMPLE_NAME}"
fi

## ---- Step 2: Run Numbat -----------------------------------------------------
echo "$(date) - [Step 2/2] Starting Numbat CNV analysis for: ${SAMPLE_NAME}"

${RSCRIPT} src/03_numbat_cnv.R ${RSCRIPT_ARGS}

if [ $? -ne 0 ]; then
    echo "$(date) - [FAILED] Numbat CNV analysis failed for: ${SAMPLE_NAME}" >&2
    exit 20
fi

echo "$(date) - [Step 2/2] Numbat CNV analysis completed for: ${SAMPLE_NAME}"
echo "$(date) - [DONE] All steps completed successfully for: ${SAMPLE_NAME}"
