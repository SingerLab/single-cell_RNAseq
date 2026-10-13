#!/usr/bin/env Rscript
#' 06_crossbuild_fingerprint_analysis.R - Cross-build (GRCh37 <-> GRCh38) genotype
#' concordance for the P22 sample-identity litmus test (see
#' docs/superpowers/specs/2026-10-08-p22-crossbuild-fingerprint-design.md).
#'
#' Computes two independent comparisons per P22 scRNA library against the
#' s_C_FLCCA3 IMPACT BAM:
#'   Comparison A: join by site ID   (s_C_FLCCA3 genotyped in native b37 at the
#'                 lifted-down site panel; joined against P22's native hg38 calls)
#'   Comparison B: join by chrom:pos (s_C_FLCCA3's own b37 genotype calls lifted UP
#'                 to hg38; joined against P22's native hg38 calls by coordinate)
#' A and B use different chain files and different join logic -- their agreement
#' is reported as a check on the liftover itself, separate from the biological
#' concordance result. Reports concordance % and informative-site counts, NOT a
#' calibrated pass/fail threshold -- the cohort's Somalier threshold does not
#' transfer to this one-off cross-platform comparison.
#'
#' @usage
#' Rscript 06_crossbuild_fingerprint_analysis.R --s_flcca3_genotypes_a <vcf> \
#'   --s_flcca3_genotypes_b <vcf> --p22_manifest <csv> --out_dir <dir> \
#'   [--min_informative_sites <int>] [--tolerance_pct <num>]

suppressPackageStartupMessages(library(optparse))

this_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", this_file_arg)))
source(file.path(SCRIPT_DIR, "06_crossbuild_fingerprint_functions.R"))

option_list <- list(
    make_option("--s_flcca3_genotypes_a", type = "character", default = NULL,
                help = "s_C_FLCCA3 genotype VCF, native b37, ID-annotated from the lifted sites panel (Comparison A input) [required]"),
    make_option("--s_flcca3_genotypes_b", type = "character", default = NULL,
                help = "s_C_FLCCA3 genotype VCF, reverse-lifted to hg38 (Comparison B input) [required]"),
    make_option("--p22_manifest", type = "character", default = NULL,
                help = "2-column CSV (sample_id,genotype_vcf_path) for P22's native-hg38 genotype VCFs [required]"),
    make_option("--out_dir", type = "character", default = NULL,
                help = "Output directory [required]"),
    make_option("--min_informative_sites", type = "integer", default = 1366L,
                help = "Minimum informative sites for a conclusive call [default: %default, matching the cohort's existing low-coverage flag]"),
    make_option("--tolerance_pct", type = "double", default = 10,
                help = "Max allowed |Comparison A - Comparison B| concordance % before flagging liftover disagreement [default: %default; calibrate via the known-pair sanity check before trusting this on s_C_FLCCA3]")
)

opt <- parse_args(OptionParser(
    option_list = option_list,
    description = "P22 cross-build sample-identity concordance (LiftOver + genotype concordance)",
    epilogue = "Reports genotype concordance % and informative-site counts, NOT a calibrated pass/fail threshold -- see the design spec for the strong/weak/inconclusive bands."))

for (required in c("s_flcca3_genotypes_a", "s_flcca3_genotypes_b", "p22_manifest", "out_dir")) {
    if (is.null(opt[[required]])) stop("--", required, " is required. Use --help for usage.")
}
for (required_file in c(opt$s_flcca3_genotypes_a, opt$s_flcca3_genotypes_b, opt$p22_manifest)) {
    if (!file.exists(required_file)) stop("File not found: ", required_file)
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Load inputs -----------------------------------------------------------------
s_flcca3_a <- read_genotype_vcf(opt$s_flcca3_genotypes_a)
s_flcca3_b <- read_genotype_vcf(opt$s_flcca3_genotypes_b)

p22_manifest <- read.csv(opt$p22_manifest, stringsAsFactors = FALSE, colClasses = "character")
missing_cols <- setdiff(c("sample_id", "genotype_vcf_path"), names(p22_manifest))
if (length(missing_cols) > 0) {
    stop("p22_manifest is missing column(s): ", paste(missing_cols, collapse = ", "))
}

## ---- Run Comparisons A and B per P22 sample ---------------------------------------
rows <- list()
agreement_rows <- list()

for (i in seq_len(nrow(p22_manifest))) {
    sample_id <- p22_manifest$sample_id[i]
    p22_gt <- read_genotype_vcf(p22_manifest$genotype_vcf_path[i])

    result_a <- compute_concordance(join_by_id(s_flcca3_a, p22_gt), min_informative_sites = opt$min_informative_sites)
    result_b <- compute_concordance(join_by_coord(s_flcca3_b, p22_gt), min_informative_sites = opt$min_informative_sites)

    band_a <- band_concordance(result_a)
    band_b <- band_concordance(result_b)
    agreement <- check_ab_agreement(result_a, result_b, tolerance_pct = opt$tolerance_pct)

    rows[[length(rows) + 1]] <- list(sample = sample_id, comparison = "A (join by ID)",
                                      n_sites = result_a$n_sites, pct_concordant = result_a$pct_concordant,
                                      band = band_a)
    rows[[length(rows) + 1]] <- list(sample = sample_id, comparison = "B (join by chrom:pos)",
                                      n_sites = result_b$n_sites, pct_concordant = result_b$pct_concordant,
                                      band = band_b)
    agreement_rows[[length(agreement_rows) + 1]] <- data.frame(
        sample_id = sample_id, pct_concordant_a = result_a$pct_concordant,
        pct_concordant_b = result_b$pct_concordant, diff = agreement$diff, agree = agreement$agree
    )

    cat(sprintf("%s: A=%s (n=%d, %s)  B=%s (n=%d, %s)  A/B agree=%s\n",
                sample_id,
                if (is.na(result_a$pct_concordant)) "NA" else sprintf("%.1f%%", result_a$pct_concordant),
                result_a$n_sites, band_a,
                if (is.na(result_b$pct_concordant)) "NA" else sprintf("%.1f%%", result_b$pct_concordant),
                result_b$n_sites, band_b, agreement$agree))
}

## ---- Write reports -----------------------------------------------------------------
report_tsv <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
write.table(report_tsv, file.path(opt$out_dir, "p22_crossbuild_fingerprint_report.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

agreement_df <- do.call(rbind, agreement_rows)
write.table(agreement_df, file.path(opt$out_dir, "p22_crossbuild_fingerprint_ab_agreement.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write_html_report(file.path(opt$out_dir, "p22_crossbuild_fingerprint_report.html"),
                   title = "P22 cross-build sample-identity fingerprint (vs. s_C_FLCCA3)",
                   rows = rows)

cat("\nDone. Report: ", file.path(opt$out_dir, "p22_crossbuild_fingerprint_report.tsv"), "\n", sep = "")
## Per-row check: isTRUE() on a vector is always FALSE, so the brief's version warned unconditionally.
if (any(agreement_df$agree %in% FALSE)) {
    cat("WARNING: Comparison A and B disagree beyond tolerance for at least one sample -- treat results as liftover-uncertain, not a biological verdict, until investigated. See: ",
        file.path(opt$out_dir, "p22_crossbuild_fingerprint_ab_agreement.tsv"), "\n", sep = "")
}
