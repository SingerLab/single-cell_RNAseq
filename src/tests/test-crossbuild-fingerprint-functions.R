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

test_that("join_by_id drops unannotated '.'/NA ids instead of cross-joining them", {
    ## 2 "." rows on one side x 3 on the other + 1 real match: a plain merge on
    ## id would yield 2*3 + 1 = 7 rows (57% concordant); the truth is 1 site, 100%.
    ref_df <- data.frame(id = c(".", ".", "real", NA), chrom = "1",
                          pos = c(1L, 2L, 3L, 4L), gt = c("0/0", "0/1", "0/1", "1/1"),
                          stringsAsFactors = FALSE)
    query_df <- data.frame(id = c(".", ".", ".", "real", NA), chrom = "chr1",
                            pos = c(10L, 20L, 30L, 40L, 50L),
                            gt = c("1/1", "0/1", "1/1", "0/1", "1/1"),
                            stringsAsFactors = FALSE)

    joined <- join_by_id(ref_df, query_df)
    expect_equal(nrow(joined), 1)
    expect_equal(joined$id, "real")

    result <- compute_concordance(joined, min_informative_sites = 1)
    expect_equal(result$n_sites, 1)
    expect_equal(result$pct_concordant, 100)
})

test_that("count_ab_overlap counts informative sites shared by A and B on P22's hg38 chrom:pos", {
    joined_a <- data.frame(id = c("s1", "s2", "s3"), chrom_query = "chr1", pos_query = c(10L, 20L, 30L),
                            gt_ref = c("0/1", "0/0", NA), gt_query = c("0/1", "0/1", "0/1"),
                            stringsAsFactors = FALSE)
    joined_b <- data.frame(chrom = "chr1", pos = c(20L, 30L, 40L),
                            gt_ref = c("0/0", "1/1", "0/1"), gt_query = c("0/1", "0/1", "0/1"),
                            stringsAsFactors = FALSE)
    ## chr1:20 informative in both; chr1:30 not informative in A; 10/40 one side only
    expect_equal(count_ab_overlap(joined_a, joined_b), 1)
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

test_that("write_html_report includes A/B agreement, liftover rejections and notes when given", {
    tmp <- tempfile(fileext = ".html")
    agreement <- data.frame(sample_id = "P22_C1_D0_01", pct_concordant_a = 97.3,
                            pct_concordant_b = 60.0, diff = 37.3, agree = FALSE,
                            n_overlap_ab = 1234L)
    write_html_report(tmp, "Test Report", list(
        list(sample = "P22_C1_D0_01", comparison = "A (join by ID)", n_sites = 1500,
             pct_concordant = 97.3, band = "strong_match")
    ), agreement = agreement, forward_rejection = "100/17266", reverse_rejection = "NA",
    notes = c("WARNING: A/B disagree", "Caveat: tumor-only"))
    html <- paste(readLines(tmp), collapse = "\n")
    unlink(tmp)

    expect_true(grepl("A vs. B agreement", html))
    expect_true(grepl("1234", html))
    expect_true(grepl("FALSE", html))
    expect_true(grepl("100 of 17266 rejected (0.58%)", html, fixed = TRUE))
    expect_true(grepl("Reverse (hg19 -> hg38, s_C_FLCCA3 calls): not available", html, fixed = TRUE))
    expect_true(grepl("WARNING: A/B disagree", html))
    expect_true(grepl("Caveat: tumor-only", html))
})
