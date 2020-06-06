### Get gene signatures for different cell types from Novershten et. al 2011

## libraries
library(tidyverse)
library(org.Hs.eg.db)

## excel file
excel.file <- "1-s2.0-S0092867411000055-mmc2.xls"

## get sheet names
sheet.names <- readxl::excel_sheets(excel.file)

## loop through sheets to get list of genes
aa <- do.call(rbind, lapply(sheet.names[2:17], function(sht) {
    tt <- readxl::read_excel(excel.file, sheet = sht)
    tt$cell.type <- sht
    return(tt)
}))

aa$updown[grep("up", aa$cell.type)] <- "up"
aa$updown[grep("down", aa$cell.type)] <- "down"

aa$cell.type <- gsub("_(up|down)", "", aa$cell.type)

ssi <- data.frame(aa)

## thesaurus
hs <- org.Hs.eg.db
my.ids <- as.character(unique(ssi$Gene.ID))

signature.genes <- select(hs, 
                   keys = as.character(ssi$Gene.ID),
                   columns = c("ENTREZID", "SYMBOL"),
                   keytype = "ENTREZID")


ssi$hgnc.symbol <- signature.genes$SYMBOL


cell.types <- unique(ssi$cell.type)

sapply(cell.types, function(ct) {
    write.table(ssi[ssi$cell.type == ct, ], file = paste(ct, "signature.txt", sep = "_"), quote = FALSE, row.names = FALSE, sep = "\t")
})
