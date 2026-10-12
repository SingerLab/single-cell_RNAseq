#!/bin/bash
set -euo pipefail
WORK="$(mktemp -d)"
trap 'rm -rf "${WORK}"' EXIT
cd "${WORK}"

cat > target.fasta <<'EOF'
>chr1
ACGTACGTACGTACGTACGTACGTACGTAC
EOF
samtools faidx target.fasta
picard CreateSequenceDictionary R=target.fasta O=target.dict >/dev/null 2>&1

## Single ungapped block, target shifted +5 vs source -- source pos 10 -> target pos 15.
cat > shift.chain <<'EOF'
chain 1000 chr1 30 + 0 25 chr1 30 + 5 30 1
25
EOF

## Two sites: one with a REF allele matching the target base at its lifted position
## (should lift successfully), one deliberately mismatched (should be rejected) --
## target.fasta base at pos 15 is 'G' (1-based).
cat > input.vcf <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=30>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	10	siteA	G	C	.	.	.	GT	0/1
chr1	11	siteB	A	T	.	.	.	GT	1/1
EOF

picard LiftoverVcf I=input.vcf O=lifted.vcf.gz CHAIN=shift.chain REJECT=rejected.vcf.gz \
    R=target.fasta WARN_ON_MISSING_CONTIG=true CREATE_INDEX=true >picard.log 2>&1

## gzip -dc (not zcat) -- BSD zcat (macOS) refuses .gz input without a .Z
## extension; grep -vc exits 1 on a zero count, so `|| true` absorbs that
## without masking a real failure (nothing else in the pipe can fail here).
n_input=$(grep -vc '^#' input.vcf)
n_rejected=$(gzip -dc rejected.vcf.gz 2>/dev/null | { grep -vc '^#' || true; })
reject_pct=$(awk -v r="${n_rejected}" -v n="${n_input}" 'BEGIN { printf "%.2f", (n > 0) ? (100 * r / n) : 0 }')

[ "${n_rejected}" -eq 1 ] || { echo "FAIL: expected 1 rejected site (siteB, mismatched REF), got ${n_rejected}"; cat picard.log; exit 1; }

gzip -dc lifted.vcf.gz | grep -q "^chr1	15	siteA	" || { echo "FAIL: siteA did not lift to chr1:15"; exit 1; }

## Reject-rate abort logic: 1/2 = 50% should trip a 20% max; 1/2 should NOT trip a 60% max.
awk -v p="${reject_pct}" -v max=20 'BEGIN { exit !(p > max) }' \
    && echo "PASS: 50% reject rate correctly exceeds max_reject_pct=20" \
    || { echo "FAIL: 50% reject rate should have exceeded max_reject_pct=20"; exit 1; }

awk -v p="${reject_pct}" -v max=60 'BEGIN { exit (p > max) }' \
    && echo "PASS: 50% reject rate correctly does NOT exceed max_reject_pct=60" \
    || { echo "FAIL: 50% reject rate should not have exceeded max_reject_pct=60"; exit 1; }

echo "PASS: test_liftover_sites_fixture"
