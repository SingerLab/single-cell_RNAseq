#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANALYSIS_SCRIPT="${SCRIPT_DIR}/../06_crossbuild_fingerprint_analysis.R"
WORK="$(mktemp -d)"
trap 'rm -rf "${WORK}"' EXIT

HDR=$'##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE'

## s_C_FLCCA3, Comparison A input (b37 coords, ID-keyed) -- matches P22 at 2/2 sites
printf '%s\n1\t900\tsiteA\tA\tG\t.\t.\t.\tGT\t0/1\n1\t1900\tsiteB\tC\tT\t.\t.\t.\tGT\t1/1\n' "${HDR}" > "${WORK}/s_flcca3_a.vcf"

## s_C_FLCCA3, Comparison B input (reverse-lifted to hg38 coords, same genotypes)
printf '%s\nchr1\t1000\tsiteA\tA\tG\t.\t.\t.\tGT\t0/1\nchr1\t2000\tsiteB\tC\tT\t.\t.\t.\tGT\t1/1\n' "${HDR}" > "${WORK}/s_flcca3_b.vcf"

## P22 sample, native hg38 -- concordant with s_C_FLCCA3 at both sites, both ways
printf '%s\nchr1\t1000\tsiteA\tA\tG\t.\t.\t.\tGT\t0/1\nchr1\t2000\tsiteB\tC\tT\t.\t.\t.\tGT\t1/1\n' "${HDR}" > "${WORK}/p22_sample1.vcf"

printf 'sample_id,genotype_vcf_path\nP22_TEST_01,%s\n' "${WORK}/p22_sample1.vcf" > "${WORK}/p22_manifest.csv"

Rscript "${ANALYSIS_SCRIPT}" \
    --s_flcca3_genotypes_a "${WORK}/s_flcca3_a.vcf" \
    --s_flcca3_genotypes_b "${WORK}/s_flcca3_b.vcf" \
    --p22_manifest "${WORK}/p22_manifest.csv" \
    --out_dir "${WORK}/out" \
    --min_informative_sites 1 \
    > "${WORK}/stdout.log" 2>&1 || { echo "FAIL: Rscript exited non-zero"; cat "${WORK}/stdout.log"; exit 1; }

grep -q "A=100.0%" "${WORK}/stdout.log" || { echo "FAIL: expected 100% concordance for Comparison A"; cat "${WORK}/stdout.log"; exit 1; }
grep -q "B=100.0%" "${WORK}/stdout.log" || { echo "FAIL: expected 100% concordance for Comparison B"; cat "${WORK}/stdout.log"; exit 1; }
grep -q "agree=TRUE" "${WORK}/stdout.log" || { echo "FAIL: expected A/B agreement"; cat "${WORK}/stdout.log"; exit 1; }
if grep -q "WARNING" "${WORK}/stdout.log"; then echo "FAIL: unexpected A/B disagreement warning"; cat "${WORK}/stdout.log"; exit 1; fi

[ -f "${WORK}/out/p22_crossbuild_fingerprint_report.tsv" ] || { echo "FAIL: report.tsv not written"; exit 1; }
[ -f "${WORK}/out/p22_crossbuild_fingerprint_ab_agreement.tsv" ] || { echo "FAIL: ab_agreement.tsv not written"; exit 1; }
[ -f "${WORK}/out/p22_crossbuild_fingerprint_report.html" ] || { echo "FAIL: report.html not written"; exit 1; }

grep -q "tumor-only" "${WORK}/out/p22_crossbuild_fingerprint_report.html" || { echo "FAIL: tumor-only caveat missing from HTML"; exit 1; }
grep -q "A vs. B agreement" "${WORK}/out/p22_crossbuild_fingerprint_report.html" || { echo "FAIL: A/B agreement section missing from HTML"; exit 1; }
if grep -q "WARNING" "${WORK}/out/p22_crossbuild_fingerprint_report.html"; then echo "FAIL: unexpected A/B warning in HTML"; exit 1; fi

## ---- A/B disagreement branch: B genotypes deliberately flipped at both sites, so
## A=100% but B=0% -- must trip the warning on stdout AND in the HTML (this branch
## shipped the isTRUE()-on-a-vector bug once).
printf '%s\nchr1\t1000\tsiteA\tA\tG\t.\t.\t.\tGT\t1/1\nchr1\t2000\tsiteB\tC\tT\t.\t.\t.\tGT\t0/0\n' "${HDR}" > "${WORK}/s_flcca3_b_mismatch.vcf"

Rscript "${ANALYSIS_SCRIPT}" \
    --s_flcca3_genotypes_a "${WORK}/s_flcca3_a.vcf" \
    --s_flcca3_genotypes_b "${WORK}/s_flcca3_b_mismatch.vcf" \
    --p22_manifest "${WORK}/p22_manifest.csv" \
    --out_dir "${WORK}/out_mismatch" \
    --min_informative_sites 1 \
    --forward_liftover_rejected "3/100" \
    --reverse_liftover_rejected "1/2" \
    > "${WORK}/stdout_mismatch.log" 2>&1 || { echo "FAIL: Rscript exited non-zero (mismatch run)"; cat "${WORK}/stdout_mismatch.log"; exit 1; }

grep -q "B=0.0%" "${WORK}/stdout_mismatch.log" || { echo "FAIL: expected 0% concordance for mismatched Comparison B"; cat "${WORK}/stdout_mismatch.log"; exit 1; }
grep -q "agree=FALSE" "${WORK}/stdout_mismatch.log" || { echo "FAIL: expected A/B disagreement"; cat "${WORK}/stdout_mismatch.log"; exit 1; }
grep -q "^WARNING: Comparison A and B disagree" "${WORK}/stdout_mismatch.log" || { echo "FAIL: A/B disagreement warning not printed"; cat "${WORK}/stdout_mismatch.log"; exit 1; }
MISMATCH_HTML="${WORK}/out_mismatch/p22_crossbuild_fingerprint_report.html"
grep -q "WARNING: Comparison A and B disagree" "${MISMATCH_HTML}" || { echo "FAIL: A/B disagreement warning missing from HTML"; exit 1; }
grep -qF "3 of 100 rejected (3.00%)" "${MISMATCH_HTML}" || { echo "FAIL: forward rejection count missing from HTML"; exit 1; }
grep -qF "1 of 2 rejected (50.00%)" "${MISMATCH_HTML}" || { echo "FAIL: reverse rejection count missing from HTML"; exit 1; }

echo "PASS: test_crossbuild_fingerprint_analysis_cli"
