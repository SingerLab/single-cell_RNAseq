#!/usr/bin/env Rscript
#' 06_crossbuild_fingerprint_functions.R - Pure helper functions for the P22
#' cross-build (GRCh37 <-> GRCh38) genotype-concordance fingerprint.
#'
#' Sourced by 06_crossbuild_fingerprint_analysis.R and by this file's own
#' testthat suite (src/tests/test-crossbuild-fingerprint-functions.R). No
#' top-level code runs when this file is sourced -- only function definitions.

normalize_gt <- function(gt) {
    if (is.na(gt) || gt %in% c(".", "./.", ".|.")) return(NA_character_)
    alleles <- strsplit(gt, "[/|]")[[1]]
    if (length(alleles) != 2 || any(alleles == ".")) return(NA_character_)
    paste(sort(alleles), collapse = "/")
}

read_genotype_vcf <- function(path) {
    lines <- readLines(path)
    header_idx <- which(startsWith(lines, "#CHROM"))
    if (length(header_idx) != 1) {
        stop("Expected exactly one #CHROM header line in ", path, ", found ", length(header_idx))
    }
    header <- strsplit(sub("^#", "", lines[header_idx]), "\t")[[1]]
    ## lines[-seq_len(header_idx)], not lines[(header_idx+1):length(lines)] --
    ## the latter counts BACKWARDS in R (e.g. 3:2) when the header is the last
    ## line, silently reintroducing the header line plus an out-of-bounds NA
    ## as fake "data". Caught by the header-only-VCF test below.
    data_lines <- lines[-seq_len(header_idx)]
    data_lines <- data_lines[nzchar(data_lines)]

    if (length(data_lines) == 0) {
        return(data.frame(id = character(0), chrom = character(0),
                           pos = integer(0), gt = character(0), stringsAsFactors = FALSE))
    }

    fields <- do.call(rbind, strsplit(data_lines, "\t"))
    colnames(fields) <- header
    sample_col <- header[length(header)]

    gt_raw <- vapply(seq_len(nrow(fields)), function(i) {
        format_fields <- strsplit(fields[i, "FORMAT"], ":")[[1]]
        gt_idx <- match("GT", format_fields)
        if (is.na(gt_idx)) return(NA_character_)
        strsplit(fields[i, sample_col], ":")[[1]][gt_idx]
    }, character(1))

    data.frame(
        id = fields[, "ID"],
        chrom = fields[, "CHROM"],
        pos = as.integer(fields[, "POS"]),
        gt = vapply(gt_raw, normalize_gt, character(1)),
        stringsAsFactors = FALSE
    )
}

join_by_id <- function(ref_df, query_df) {
    ## "." (or NA) means bcftools annotate assigned no ID -- not a key. Merging
    ## on it would pair every "." row with every "." row on the other side.
    has_id <- function(df) df[!is.na(df$id) & df$id != ".", , drop = FALSE]
    merge(has_id(ref_df), has_id(query_df), by = "id", suffixes = c("_ref", "_query"))
}

join_by_coord <- function(ref_df, query_df) {
    merge(ref_df, query_df, by = c("chrom", "pos"), suffixes = c("_ref", "_query"))
}

## Informative sites shared by Comparison A and B, keyed on the P22 (query)
## side's native hg38 chrom:pos -- the one coordinate system both joins share.
count_ab_overlap <- function(joined_a, joined_b) {
    inf_a <- !is.na(joined_a$gt_ref) & !is.na(joined_a$gt_query)
    inf_b <- !is.na(joined_b$gt_ref) & !is.na(joined_b$gt_query)
    length(intersect(paste(joined_a$chrom_query[inf_a], joined_a$pos_query[inf_a], sep = ":"),
                     paste(joined_b$chrom[inf_b], joined_b$pos[inf_b], sep = ":")))
}

compute_concordance <- function(joined_df, min_informative_sites = 1366) {
    informative_rows <- joined_df[!is.na(joined_df$gt_ref) & !is.na(joined_df$gt_query), ]
    n_sites <- nrow(informative_rows)
    n_concordant <- sum(informative_rows$gt_ref == informative_rows$gt_query)
    pct_concordant <- if (n_sites > 0) 100 * n_concordant / n_sites else NA_real_

    list(
        n_sites = n_sites,
        n_concordant = n_concordant,
        pct_concordant = pct_concordant,
        informative = n_sites >= min_informative_sites
    )
}

band_concordance <- function(result) {
    if (!isTRUE(result$informative)) return("inconclusive")
    if (is.na(result$pct_concordant)) return("inconclusive")
    if (result$pct_concordant >= 95) return("strong_match")
    if (result$pct_concordant < 60) return("clearly_different")
    "inconclusive"
}

check_ab_agreement <- function(result_a, result_b, tolerance_pct = 10) {
    if (is.na(result_a$pct_concordant) || is.na(result_b$pct_concordant)) {
        return(list(agree = NA, diff = NA_real_))
    }
    diff <- abs(result_a$pct_concordant - result_b$pct_concordant)
    list(agree = diff <= tolerance_pct, diff = diff)
}

## "rejected/total" (as passed by 06_crossbuild_fingerprint.sh) -> readable text.
format_rejection <- function(x) {
    if (is.null(x) || is.na(x) || !grepl("^[0-9]+/[0-9]+$", x)) return("not available")
    n <- as.numeric(strsplit(x, "/")[[1]])
    sprintf("%.0f of %.0f rejected (%.2f%%)", n[1], n[2], if (n[2] > 0) 100 * n[1] / n[2] else 0)
}

write_html_report <- function(path, title, rows, agreement = NULL,
                              forward_rejection = NULL, reverse_rejection = NULL,
                              notes = character(0)) {
    fmt_pct <- function(x) if (is.na(x)) "NA" else sprintf("%.1f%%", x)
    row_html <- vapply(rows, function(r) {
        sprintf(
            "<tr><td>%s</td><td>%s</td><td>%d</td><td>%s</td><td>%s</td></tr>",
            r$sample, r$comparison, r$n_sites,
            fmt_pct(r$pct_concordant),
            r$band
        )
    }, character(1))

    agreement_html <- ""
    if (!is.null(agreement) && nrow(agreement) > 0) {
        agreement_rows <- vapply(seq_len(nrow(agreement)), function(i) {
            a <- agreement[i, ]
            sprintf("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%d</td></tr>",
                    a$sample_id, fmt_pct(a$pct_concordant_a), fmt_pct(a$pct_concordant_b),
                    if (is.na(a$diff)) "NA" else sprintf("%.1f", a$diff), a$agree, a$n_overlap_ab)
        }, character(1))
        agreement_html <- paste0(
            "<h2>Comparison A vs. B agreement</h2>\n",
            "<table border=\"1\"><tr><th>Sample</th><th>A concordance</th><th>B concordance</th>",
            "<th>|A - B| (pct points)</th><th>Agree</th><th>Informative sites in both A and B</th></tr>\n",
            paste(agreement_rows, collapse = "\n"), "\n</table>\n")
    }

    rejection_html <- paste0(
        "<h2>Liftover rejections</h2>\n<ul>",
        "<li>Forward (hg38 -> hg19, site panel): ", format_rejection(forward_rejection), "</li>",
        "<li>Reverse (hg19 -> hg38, s_C_FLCCA3 calls): ", format_rejection(reverse_rejection), "</li>",
        "</ul>\n")

    notes_html <- paste0("<p>", notes, "</p>\n", collapse = "")

    html <- paste0(
        "<html><head><title>", title, "</title></head><body>\n",
        "<h1>", title, "</h1>\n",
        "<table border=\"1\"><tr><th>Sample</th><th>Comparison</th><th>N sites</th>",
        "<th>Concordance</th><th>Band</th></tr>\n",
        paste(row_html, collapse = "\n"),
        "\n</table>\n",
        agreement_html, rejection_html, notes_html,
        "</body></html>"
    )
    writeLines(html, path)
    invisible(path)
}
