#!/bin/bash
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TARGET="${SCRIPT_DIR}/../06_crossbuild_fingerprint.sh"

"${TARGET}" 2>/tmp/crossbuild_args_test.err
status=$?
grep -q "lifted_sites_vcf is required" /tmp/crossbuild_args_test.err || { echo "FAIL: missing-arg error message not shown"; exit 1; }
[ "${status}" -ne 0 ] || { echo "FAIL: expected non-zero exit with no arguments"; exit 1; }

"${TARGET}" /nonexistent/sites.vcf.gz /nonexistent/stable.vcf.gz /nonexistent/s.bam \
    /nonexistent/b37.fasta /nonexistent/hg38.fasta /nonexistent/chain.gz \
    /nonexistent/manifest.csv P22 /tmp/crossbuild_args_test_out 2>/tmp/crossbuild_args_test2.err
status=$?
grep -q "Required file not found" /tmp/crossbuild_args_test2.err || { echo "FAIL: missing-file error message not shown"; cat /tmp/crossbuild_args_test2.err; exit 1; }
[ "${status}" -ne 0 ] || { echo "FAIL: expected non-zero exit with nonexistent input files"; exit 1; }

rm -f /tmp/crossbuild_args_test.err /tmp/crossbuild_args_test2.err
echo "PASS: test_crossbuild_fingerprint_args"
