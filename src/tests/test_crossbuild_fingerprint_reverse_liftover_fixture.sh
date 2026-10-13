#!/bin/bash
set -euo pipefail
WORK="$(mktemp -d)"
trap 'rm -rf "${WORK}"' EXIT
cd "${WORK}"

## Stand-in for a GRCh38 reference target: 1 contig, 30bp.
cat > hg38.fasta <<'EOF'
>chr1
ACGTACGTACGTACGTACGTACGTACGTAC
EOF
samtools faidx hg38.fasta
picard CreateSequenceDictionary R=hg38.fasta O=hg38.dict >/dev/null 2>&1

## Reverse chain: b37 chr1 (no "chr" prefix, as s_C_FLCCA3's native build) ->
## hg38 "chr1", shifted +5 (the mirror image of 06_liftover_sites.sh's chain).
cat > reverse.chain <<'EOF'
chain 1000 1 30 + 0 25 chr1 30 + 5 30 1
25
EOF

## Stand-in for s_C_FLCCA3.genotypes.b37.vcf.gz (Comparison A input) -- REF
## matches hg38.fasta's base at the lifted target position (pos 15 = 'G').
cat > s_flcca3_a.vcf <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=1,length=30>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	10	siteA	G	C	.	.	.	GT	0/1
EOF

picard LiftoverVcf I=s_flcca3_a.vcf O=s_flcca3_b.vcf.gz CHAIN=reverse.chain \
    REJECT=s_flcca3_b.rejected.vcf.gz R=hg38.fasta WARN_ON_MISSING_CONTIG=true \
    CREATE_INDEX=true >picard.log 2>&1

gzip -dc s_flcca3_b.vcf.gz | grep -q "^chr1	15	siteA	" \
    || { echo "FAIL: siteA did not reverse-lift to chr1:15"; cat picard.log; exit 1; }

echo "PASS: test_crossbuild_fingerprint_reverse_liftover_fixture"
