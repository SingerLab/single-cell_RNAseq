#' read multiple 10x outputs
#' @param sample.dir sample directory containing  `outs/filtered_feature_bc_matrix`
#' @param sample.names names to use for the samples e.g. directory name w/o path
read.m.10x <- function(sample.dir, sample.names) {
    ## read all files
    message(date(), "  start")
    a.data <- lapply(file.path(sample.dir, "outs/filtered_feature_bc_matrix/"),
                  function(i) {
                      out <- Read10X(i)
                      message(date(), "  finished reading: ", i)
                      return(out)
                  })
    
    names(a.data) <- sample.names

    return(a.data)
}


#' build s.data ; list of seurat objects
#' @param a.data read counts from read.m.10x list form
#' @param sample.names name of the samples if different than a.data, default is NULL
#' 
#' @importFrom DoubletFinder dubFinder
#' @export
build.seurat.list <- function(a.data, sample.names = NULL, dbf = TRUE,
                              max.mt.pct = NULL, min.genes = 200) {
                              
    if(is.null(sample.names)) {
        sample.names <- names(a.data)
    }

    s.data <- lapply(sample.names, function(sm) {
        print(sm)
        ## run doublet finder or jsut CreateSeuratObject
        if(dbf) {
            ## create dbf temporary output directory
            dbfPath = "tmp/00.dbf/"
            if(!dir.exists(dbfPath)) { dir.create(dbfPath, recursive = TRUE) }
            ## run dbf
            pso <- dubFinder(a.data[[sm]], sm, png.path = dbfPath,
                             max.mt.pct = max.mt.pct, min.genes = min.genes)
            
        } else {
            
            pso <- CreateSeuratObject(a.data[[sm]], project = sm)
            
        }

        ## parse sample info from orig.ident
        cx <- as.character(pso@meta.data$orig.ident)
        sAnnot <- data.frame(do.call(rbind, strsplit(cx, split = "_")))

        ## names depending if it has 3, 4, or 5 columns
        if(ncol(sAnnot) == 3) {
            names(sAnnot) <- c("patient.id", "cycle", "day")
        } else {
            if(ncol(sAnnot) == 4) {
                names(sAnnot) <- c("patient.id", "cycle", "day", "replicate")
            } else {
                if(ncol(sAnnot) == 5) {
                    names(sAnnot) <- c("patient.id", "cycle", "day",
                                       "region", "replicate")
                }
            }
        }
        if(! "region" %in% names(sAnnot)) {
            sAnnot$region <- NA
        }
        
        columns.remove <- c("percent.mt", "nCount_SCT", "nFeature_SCT",
                            "SCT_snn_res.0.8", "seurat_clusters")
        columns.keep <- setdiff(names(pso@meta.data), columns.remove)
        
        pso@meta.data <- pso@meta.data[, columns.keep]

        pso@meta.data <- cbind(pso@meta.data, sAnnot)

        DefaultAssay(pso) <- "RNA"
        pso <- DietSeurat(pso, assays = "RNA",
                          misc = FALSE)
        
        return(pso)
    })
    names(s.data) <- sample.names

    return(s.data)
    
}

#' prepare a joined seurat object from a list of files
#'
#' @export
prepare_seu <- function(a.data, sample.names = NULL,
                        project = "Palbo_aPD1",
                        dbf = FALSE,
                        nbin = 24,
                        seed = 2021, ...) {
    
    if(is.null(sample.names)) {
        sample.names <- names(a.data)
    }

    message(Sys.time(), "  start")
    s.data <- build.seurat.list(a.data, sample.names = sample.names, dbf = dbf, ...)

    message(Sys.time(), "  merging seurat")
    ## merge data w/o doublet tagging
    pre.seu <- merge(s.data[[1]], s.data[2:length(s.data)],
                     add.cell.ids = sample.names,
                     project = project) %>%
        PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
        PercentageFeatureSet(pattern = "^RP[SL]", col.name = "percent.rb") %>%
        PercentageFeatureSet(pattern = "^HB[^(P)]", col.name = "percent.hb") %>%
        PercentageFeatureSet(pattern = "PECAM1|PF4", col.name = "percent.plat") %>%
        CellCycleScoring(s.features = cc.genes.updated.2019$s.genes,
                         g2m.features = cc.genes.updated.2019$g2m.genes,
                         nbin = nbin,
                         seed = seed)
    pre.seu$CC.Difference <- pre.seu$S.Score - pre.seu$G2M.Score

    message(Sys.time(), "  new metadata")
    pre.seu@meta.data <- pre.seu@meta.data %>%
        mutate(
            bioID = patient.id,
            ## new columns
            biopsy.id = ifelse(is.na(region),
                               paste(patient.id,
                                     cycle, day, sep = "_"),
                               paste(paste(patient.id, cycle,
                                           day, region, sep = "_"))),
            mt.cat = ifelse(as.numeric(percent.mt) < 24, "<24",
                     ifelse(as.numeric(percent.mt) > 60, ">60",
                            ">24")),
            rb.cat = ifelse(as.numeric(percent.rb) > 50, ">50", "<50"),
            hb.cat = ifelse(as.numeric(percent.hb) > 3, ">3", "<3"),
            plat.cat = ifelse(as.numeric(percent.plat) > 1, ">1", "<1"),
            active.cell = ifelse(nFeature_RNA > 200 & nFeature_RNA < 6000,
                               TRUE, FALSE),
            ## class change
            Phase = factor(Phase, c("G1", "S", "G2M")),
            mt.cat = factor(mt.cat, c("<24", ">24", ">60")),
            rb.cat = factor(rb.cat, c(">50", "<50")),
            hb.cat = factor(hb.cat, c(">3", "<3")),
            plat.cat = factor(plat.cat, c(">1", "<1"))
        )
    
    
    return(pre.seu)
} ## end prepare_seu


#' ----build_seurat_object_with_doublet_detection-------------------------------
prepare_seu_dbf <- function(a.data,  ...) {
    pre.seu.out <- prepare_seu(a.data = a.data, dbf = TRUE,
                               ...)
    return(pre.seu.out)
}
                            

#' ----convert_to_monocle, eval=TRUE--------------------------------------------
#'
#' @param seu seurat object
#' @param align logical, align using align_cds
#' @param alignment_group alignment variable
#' @param k  k 
#' @param min_expr minimum expresssion
#' @param round_expr round expression
#' 
#' @export
convert_to_monocle <- function(seu,
                               k = 10, min_expr = 0.01,
                               round_exprs = FALSE,  ...) {
    ## convert seurat object into cds ; cluster cells, and learn_graph
    cds <- as.cell_data_set(seu) %>%
        cluster_cells(k = k) %>%
        learn_graph()
    ## these were missing from the conversion and are required for
    ## monocle to function
    ## detect expressed genes ## by default 0.01
    cds <- detect_genes(cds, min_expr = min_expr)
    ## estimate size factors for samples
    cds <- estimate_size_factors(cds, round_exprs = round_exprs)
    ## row data must have a column named gene_short_name, use rownames
    ## we loose ensembl annotation which is more comprehensive than seurat
    rowData(cds)$gene_short_name <- rownames(rowData(cds))
    ## append cluster and partition to the metadata
    pData(cds)$cluster <- clusters(cds)
    pData(cds)$partition <- partitions(cds)

    return(cds)
}


