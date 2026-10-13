# chain_strip_chr.awk - strip the UCSC "chr" prefix from ONE side of a liftOver
# chain file, so the chain's contig names match the b37 (bare 1..22, X, Y) side
# of a cross-build liftover. Only the chain's own bookkeeping is renamed; the
# VCF/BAM data files are never touched.
#
# Usage: gzip -dcf <chain[.gz]> | awk -v col=<3|8> -f chain_strip_chr.awk > out.chain
#
# Chain header layout (UCSC chain format, verified empirically against Picard
# LiftoverVcf on a synthetic chain):
#   chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
# For liftOver, field 3 (tName) is the SOURCE assembly (the input VCF's contigs)
# and field 8 (qName) is the DESTINATION assembly (the target reference's .dict).
#   col=3: hg19ToHg38 chain, so its source matches s_C_FLCCA3's bare b37 contigs
#   col=8: hg38ToHg19 chain, so its destination matches Homo_sapiens_assembly19's .dict
#
# Only chr1..chr22, chrX, chrY are renamed: those are sequence-identical between
# UCSC hg19 and b37. chrM is deliberately NOT mapped to MT -- UCSC hg19 chrM
# (NC_001807) is a different mitochondrial sequence from b37 MT (rCRS,
# NC_012920), so renaming would place sites at wrong coordinates. chrM and
# alt/random/Un contigs keep their UCSC names, have no b37 match, and are
# counted as liftover rejections (visible in the reject-rate check).
# Alignment-block lines (not starting with "chain") pass through untouched.
$1 == "chain" && $col ~ /^chr([0-9]+|X|Y)$/ { sub(/^chr/, "", $col) }
{ print }
