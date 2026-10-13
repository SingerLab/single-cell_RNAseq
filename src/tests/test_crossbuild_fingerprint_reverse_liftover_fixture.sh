#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "${WORK}"' EXIT
cd "${WORK}"

## Models the real reverse liftover in 06_crossbuild_fingerprint.sh Step 4:
##   - s_C_FLCCA3's calls are on b37 contigs with BARE names ("1"), as in its BAM.
##   - the hg38 target (CellRanger GRCh38-2020-A genome.fa) is chr-prefixed ("chr1").
##   - the downloaded UCSC hg19ToHg38 chain is chr-prefixed on BOTH sides
##     ("chr1" -> "chr1"), so its source side matches nothing in the b37 VCF
##     until the script rewrites header field 3 with lib/chain_strip_chr.awk.

## Stand-in for a GRCh38 reference target: 1 contig, 30bp.
cat > hg38.fasta <<'FA'
>chr1
ACGTACGTACGTACGTACGTACGTACGTAC
FA
samtools faidx hg38.fasta
picard CreateSequenceDictionary R=hg38.fasta O=hg38.dict >/dev/null 2>&1

## UCSC-style hg19ToHg38 chain (chr on both sides), shifted +5.
cat > ucsc.chain <<'CH'
chain 1000 chr1 30 + 0 25 chr1 30 + 5 30 1
25
CH
awk -v col=3 -f "${SCRIPT_DIR}/../lib/chain_strip_chr.awk" ucsc.chain > reverse.chain
grep -q "^chain 1000 1 30 + 0 25 chr1 30 + 5 30 1$" reverse.chain \
    || { echo "FAIL: chain source side not rewritten to bare b37 name"; cat reverse.chain; exit 1; }

## Stand-in for s_C_FLCCA3.genotypes.b37.vcf.gz (Comparison A input) -- REF
## matches hg38.fasta's base at the lifted target position. siteB is a hom-ref
## call with ALT "." as `bcftools call -m` (no -v) emits; it must lift too.
cat > s_flcca3_a.vcf <<'VCF'
##fileformat=VCFv4.2
##contig=<ID=1,length=30>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	10	siteA	G	C	.	.	.	GT	0/1
1	11	siteB	T	.	.	.	.	GT	0/0
VCF

picard LiftoverVcf I=s_flcca3_a.vcf O=s_flcca3_b.vcf.gz CHAIN=reverse.chain \
    REJECT=s_flcca3_b.rejected.vcf.gz R=hg38.fasta WARN_ON_MISSING_CONTIG=true \
    CREATE_INDEX=true >picard.log 2>&1

gzip -dc s_flcca3_b.vcf.gz | grep -q "^chr1	15	siteA	" \
    || { echo "FAIL: siteA did not reverse-lift to chr1:15"; cat picard.log; exit 1; }
gzip -dc s_flcca3_b.vcf.gz | grep -q "^chr1	16	siteB	T	\.	" \
    || { echo "FAIL: hom-ref siteB did not reverse-lift to chr1:16"; cat picard.log; exit 1; }

## Negative control: the raw UCSC chain (source "chr1") must lift nothing from
## bare-named b37 calls -- the bug the chain rewrite exists to prevent.
picard LiftoverVcf I=s_flcca3_a.vcf O=raw.vcf.gz CHAIN=ucsc.chain \
    REJECT=raw.rejected.vcf.gz R=hg38.fasta WARN_ON_MISSING_CONTIG=true >picard_raw.log 2>&1 || true
n_raw_lifted=$(gzip -dc raw.vcf.gz 2>/dev/null | { grep -vc '^#' || true; })
[ "${n_raw_lifted}" -eq 0 ] || { echo "FAIL: raw chr-named chain unexpectedly lifted ${n_raw_lifted} b37 call(s)"; exit 1; }

echo "PASS: test_crossbuild_fingerprint_reverse_liftover_fixture"
