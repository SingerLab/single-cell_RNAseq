#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AWK_SCRIPT="${SCRIPT_DIR}/../lib/stabilize_vcf_ids.awk"

TMP_IN="$(mktemp)"
cat > "${TMP_IN}" <<'EOF'
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	rs123	A	G	.	.	.
chr1	200	.	C	T	.	.	.
EOF

OUTPUT="$(awk -f "${AWK_SCRIPT}" "${TMP_IN}")"
rm -f "${TMP_IN}"

echo "${OUTPUT}" | grep -q "^chr1	100	rs123	" || { echo "FAIL: populated ID was altered"; exit 1; }
echo "${OUTPUT}" | grep -q "^chr1	200	hg38_chr1_200	" || { echo "FAIL: blank ID was not stabilized"; exit 1; }
echo "PASS: stabilize_vcf_ids.awk"
