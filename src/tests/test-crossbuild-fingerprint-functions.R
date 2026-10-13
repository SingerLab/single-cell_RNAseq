suppressPackageStartupMessages(library(testthat))

## Run from the repo root: Rscript -e "testthat::test_file('src/tests/test-crossbuild-fingerprint-functions.R')"
## testthat::test_file() changes the working directory to the TEST FILE's own
## directory while sourcing it (verified: getwd() inside test_file() == the
## directory containing this file) -- so the path below is relative to
## src/tests/, not the repo root. This differs from
## 06_crossbuild_fingerprint_analysis.R, which is always invoked directly via
## `Rscript <path>` and uses a commandArgs()-based self-locating path instead.
source("../06_crossbuild_fingerprint_functions.R")

test_that("normalize_gt sorts alleles and nulls out missing calls", {
    expect_equal(normalize_gt("1/0"), "0/1")
    expect_equal(normalize_gt("1|0"), "0/1")
    expect_equal(normalize_gt("0/0"), "0/0")
    expect_true(is.na(normalize_gt("./.")))
    expect_true(is.na(normalize_gt(NA_character_)))
})

test_that("read_genotype_vcf parses CHROM/POS/ID/GT and normalizes genotypes", {
    tmp <- tempfile(fileext = ".vcf")
    writeLines(c(
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
        "chr1\t1000\tsiteA\tA\tG\t.\t.\t.\tGT\t0/1",
        "chr1\t2000\tsiteB\tC\tT\t.\t.\t.\tGT\t1/1",
        "chr1\t3000\tsiteC\tG\tA\t.\t.\t.\tGT\t./."
    ), tmp)

    df <- read_genotype_vcf(tmp)
    unlink(tmp)

    expect_equal(nrow(df), 3)
    expect_equal(df$id, c("siteA", "siteB", "siteC"))
    expect_equal(df$gt, c("0/1", "1/1", NA_character_))
})

test_that("read_genotype_vcf handles a header-only VCF with no data lines", {
    tmp <- tempfile(fileext = ".vcf")
    writeLines(c(
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    ), tmp)

    df <- read_genotype_vcf(tmp)
    unlink(tmp)

    expect_equal(nrow(df), 0)
    expect_setequal(names(df), c("id", "chrom", "pos", "gt"))
})

test_that("join_by_id inner-joins on id and suffixes chrom/pos/gt", {
    ref_df <- data.frame(id = c("A", "B", "C"), chrom = c("1", "1", "1"),
                          pos = c(100L, 200L, 300L), gt = c("0/0", "0/1", "1/1"),
                          stringsAsFactors = FALSE)
    query_df <- data.frame(id = c("B", "C", "D"), chrom = c("1", "1", "1"),
                            pos = c(200L, 300L, 400L), gt = c("0/1", "0/0", "1/1"),
                            stringsAsFactors = FALSE)

    joined <- join_by_id(ref_df, query_df)

    expect_equal(sort(joined$id), c("B", "C"))
    expect_true(all(c("gt_ref", "gt_query") %in% names(joined)))
})

test_that("join_by_coord inner-joins on chrom+pos regardless of id", {
    ref_df <- data.frame(id = c("X", "Y"), chrom = c("1", "1"), pos = c(500L, 600L),
                          gt = c("0/0", "1/1"), stringsAsFactors = FALSE)
    query_df <- data.frame(id = c("Y_lifted", "Z"), chrom = c("1", "1"), pos = c(600L, 700L),
                            gt = c("1/1", "0/1"), stringsAsFactors = FALSE)

    joined <- join_by_coord(ref_df, query_df)

    expect_equal(nrow(joined), 1)
    expect_equal(joined$gt_ref, "1/1")
    expect_equal(joined$gt_query, "1/1")
})

test_that("compute_concordance ignores NA genotypes and flags low informative-site counts", {
    joined <- data.frame(gt_ref = c("0/1", "1/1", NA), gt_query = c("0/1", "0/0", "0/1"),
                          stringsAsFactors = FALSE)

    result <- compute_concordance(joined, min_informative_sites = 1)
    expect_equal(result$n_sites, 2)
    expect_equal(result$n_concordant, 1)
    expect_equal(result$pct_concordant, 50)
    expect_true(result$informative)

    result_strict <- compute_concordance(joined, min_informative_sites = 10)
    expect_false(result_strict$informative)
})

test_that("compute_concordance returns NA pct_concordant, not an error, when there are zero informative sites", {
    joined <- data.frame(gt_ref = character(0), gt_query = character(0), stringsAsFactors = FALSE)
    result <- compute_concordance(joined, min_informative_sites = 1)
    expect_equal(result$n_sites, 0)
    expect_true(is.na(result$pct_concordant))
    expect_false(result$informative)
})

test_that("band_concordance applies the spec's strong/weak/inconclusive bands", {
    expect_equal(band_concordance(list(pct_concordant = 97, informative = TRUE)), "strong_match")
    expect_equal(band_concordance(list(pct_concordant = 10, informative = TRUE)), "clearly_different")
    expect_equal(band_concordance(list(pct_concordant = 80, informative = TRUE)), "inconclusive")
    expect_equal(band_concordance(list(pct_concordant = 99, informative = FALSE)), "inconclusive")
})

test_that("check_ab_agreement flags disagreement beyond tolerance", {
    close <- check_ab_agreement(list(pct_concordant = 92), list(pct_concordant = 95), tolerance_pct = 10)
    expect_true(close$agree)

    far <- check_ab_agreement(list(pct_concordant = 50), list(pct_concordant = 100), tolerance_pct = 10)
    expect_false(far$agree)
})

test_that("write_html_report writes a table containing each row's sample and band", {
    tmp <- tempfile(fileext = ".html")
    write_html_report(tmp, "Test Report", list(
        list(sample = "P22_C1_D0_01", comparison = "A (join by ID)", n_sites = 1500,
             pct_concordant = 97.3, band = "strong_match")
    ))
    html <- paste(readLines(tmp), collapse = "\n")
    unlink(tmp)

    expect_true(grepl("P22_C1_D0_01", html))
    expect_true(grepl("strong_match", html))
    expect_true(grepl("97.3", html))
})
