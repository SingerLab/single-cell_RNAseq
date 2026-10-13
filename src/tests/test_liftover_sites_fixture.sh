#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "${WORK}"' EXIT
cd "${WORK}"

## Stand-in for Homo_sapiens_assembly19 (b37): bare contig name "1", not "chr1".
cat > target.fasta <<'EOF'
>1
ACGTACGTACGTACGTACGTACGTACGTAC
EOF
samtools faidx target.fasta
picard CreateSequenceDictionary R=target.fasta O=target.dict >/dev/null 2>&1

## Single ungapped block, target shifted +5 vs source -- source pos 10 -> target pos 15.
## Named like the real UCSC hg38ToHg19 chain: chr-prefixed on BOTH sides. Fed to
## Picard only after the script's own chain rewrite (dest side, field 8 -> "1").
cat > ucsc.chain <<'EOF'
chain 1000 chr1 30 + 0 25 chr1 30 + 5 30 1
25
EOF
awk -v col=8 -f "${SCRIPT_DIR}/../lib/chain_strip_chr.awk" ucsc.chain > shift.chain
grep -q "^chain 1000 chr1 30 + 0 25 1 30 + 5 30 1$" shift.chain \
    || { echo "FAIL: chain dest side not rewritten to bare b37 name"; cat shift.chain; exit 1; }

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

gzip -dc lifted.vcf.gz | grep -q "^1	15	siteA	" || { echo "FAIL: siteA did not lift to b37 1:15"; exit 1; }

## Negative control: the raw UCSC chain (dest "chr1", absent from the b37 dict)
## must lift nothing -- this is the bug the chain rewrite exists to prevent.
picard LiftoverVcf I=input.vcf O=raw.vcf.gz CHAIN=ucsc.chain REJECT=raw.rejected.vcf.gz \
    R=target.fasta WARN_ON_MISSING_CONTIG=true >picard_raw.log 2>&1 || true
n_raw_lifted=$(gzip -dc raw.vcf.gz 2>/dev/null | { grep -vc '^#' || true; })
[ "${n_raw_lifted}" -eq 0 ] || { echo "FAIL: raw chr-named chain unexpectedly lifted ${n_raw_lifted} site(s) onto b37"; exit 1; }
echo "PASS: unrewritten UCSC chain lifts 0 sites onto b37 (rewrite is required)"

## Reject-rate abort logic: 1/2 = 50% should trip a 20% max; 1/2 should NOT trip a 60% max.
awk -v p="${reject_pct}" -v max=20 'BEGIN { exit !(p > max) }' \
    && echo "PASS: 50% reject rate correctly exceeds max_reject_pct=20" \
    || { echo "FAIL: 50% reject rate should have exceeded max_reject_pct=20"; exit 1; }

awk -v p="${reject_pct}" -v max=60 'BEGIN { exit (p > max) }' \
    && echo "PASS: 50% reject rate correctly does NOT exceed max_reject_pct=60" \
    || { echo "FAIL: 50% reject rate should not have exceeded max_reject_pct=60"; exit 1; }

## ---- TARGET_DICT naming: run the script's own line (not a copy of it) for both extensions
dict_line="$(grep '^TARGET_DICT=' "${SCRIPT_DIR}/../06_liftover_sites.sh")"
for ref in foo.fasta foo.fa; do
    TARGET_REF_FASTA="${ref}"
    eval "${dict_line}"
    [ "${TARGET_DICT}" = "foo.dict" ] || { echo "FAIL: ${ref} -> TARGET_DICT=${TARGET_DICT}, expected foo.dict"; exit 1; }
done
echo "PASS: TARGET_DICT naming for .fasta and .fa"

echo "PASS: test_liftover_sites_fixture"
