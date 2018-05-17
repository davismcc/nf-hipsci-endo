"
Make phenotypes for single-cell QTL mapping

Usage:
make_qtl_phenotypes.R --input_file <SCE_IN> --output_prefix <PRFX_OUT> --pheno_type <PHENO> --npcs_cell <INT> --npcs_line <INT> [-hv]

Options:
-i --input_file <SCE_IN>       input .rds file containing a QC'd SCESet object
-o --output_prefix <PRFX_OUT>  prefix for output files; .csv files produced
-p --pheno_type <PHENO>        string providing the desired phenotype type [default: 'mean', 'var', 'alpha']
--npcs_cell <INT>              number of PCs to regress out at the cell level
--npcs_line <INT>              number of PCs to regress out at the line (summarised) level [default: 10] 
-h --help                      show this
-v --version                   print version and stop

The program returns phenotypes and support metadata as .tsv files.

This program requires the R packages:

Davis McCarthy
July 2017
" -> doc

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(viridis))
library(ggthemes)
library(ggbeeswarm)
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(pryr))
suppressPackageStartupMessages(library(scater))
suppressWarnings(suppressPackageStartupMessages(library(scran)))
library(DT)
library(preprocessCore)

pluri_markers <- c("ENSG00000111704_NANOG", "ENSG00000204531_POU5F1", 
                   "ENSG00000164362_TERT", "ENSG00000181449_SOX2",
                   "ENSG00000121570_DPPA4")
pluri_markers_select <- c("ENSG00000111704_NANOG", "ENSG00000204531_POU5F1", 
                          "ENSG00000181449_SOX2", "ENSG00000121570_DPPA4")
mesendo_markers <- c("ENSG00000164458_T", "ENSG00000163508_EOMES",
                     "ENSG00000185155_MIXL1")
defendo_markers <- c("ENSG00000121966_CXCR4", "ENSG00000136574_GATA4",
                     "ENSG00000141448_GATA6", "ENSG00000125798_FOXA2",
                     "ENSG00000164736_SOX17", "ENSG00000133937_GSC",
                     "ENSG00000147869_CER1", "ENSG00000132130_LHX1")
defendo_markers_select <- c("ENSG00000121966_CXCR4", "ENSG00000136574_GATA4",
                            "ENSG00000164736_SOX17", "ENSG00000141448_GATA6",
                            "ENSG00000147869_CER1", "ENSG00000132130_LHX1",
                           "ENSG00000133937_GSC")
later_markers <- c("ENSG00000152804_HHEX", "ENSG00000170608_FOXA3")
diff_genes <- c(pluri_markers, mesendo_markers, defendo_markers, later_markers)


#' @param sce an SCESet object
filter_donors <-  function(sce, min_donors = 15) {
    ncell_tbl <- table(sce$donor_long_id)
    retained_donors <- names(ncell_tbl)[ncell_tbl >= min_donors]
    sce <- sce[, sce$donor_long_id %in% retained_donors]
    sce
}

#' @param sce an SCESet object
regress_cell_factors <- function(sce, npcs = 0) {
    sce <- plotPCA(sce, ncomponents = 20, return_SCESet = TRUE,
                   draw_plot = FALSE)
    tot_feats_sqrd <- scale(sce$total_features^2)
    total_features_norm <- scale(sce$total_features)
    prop_counts_top_500_features <- scale(sce$pct_counts_top_500_features)
    design <- model.matrix(~prop_counts_top_500_features + total_features_norm +
                               tot_feats_sqrd + log10_counts_endogenous_features,
                           data = pData(sce))
    if (npcs >= 1) {
        design <- cbind(design, redDim(sce)[, 1:npcs, drop = FALSE])
    }
    sce <- normalizeExprs(sce, method = "none", design = design,
                          exprs_values = "exprs", return_norm_as_exprs = FALSE)
    ## comparative PCA plots
    ## p0 <- plotPCA(sce, colour_by = "experiment", exprs_value = "exprs") +
    ##     ggtitle("PCA using exprs values")
    ## p2 <- plotPCA(sce, colour_by = "experiment", exprs_value = "norm_exprs") +
    ##     ggtitle("PCA using exprs residuals")
    ## plot_grid(p0, p2, labels = letters, ncol = 2)
    ## p0 <- plotPCA(sce, colour_by = "ENSG00000204531_POU5F1",
    ##               size_by = "total_counts", exprs_value = "exprs") +
    ##     ggtitle("PCA using exprs values")
    ## p1 <- plotPCA(sce, colour_by = "ENSG00000204531_POU5F1",
    ##               size_by = "total_counts", exprs_value = "norm_exprs") +
    ##     ggtitle("PCA using exprs resids")
    ## plot_grid(p0, p1, labels = letters)
    ## p0 <- plotPCA(sce, colour_by = "ENSG00000204531_POU5F1",
    ##               size_by = "total_features", exprs_value = "exprs") +
    ##     ggtitle("PCA using exprs values")
    ## p1 <- plotPCA(sce, colour_by = "ENSG00000204531_POU5F1",
    ##               size_by = "total_features", exprs_value = "norm_exprs") +
    ##     ggtitle("PCA using exprs resids")
    ## plot_grid(p0, p1, labels = letters)
    ## QC plots
    ## plotQC(sce, type = "find-pcs", variable = "total_features",
    ##        exprs_values = "norm_exprs")
    ## plotQC(sce, type = "find-pcs", variable = "total_counts",
    ##        exprs_values = "norm_exprs")
    ## plotQC(sce, type = "find-pcs", variable = "total_features",
    ##        exprs_values = "exprs")
    sce
}


#' @param sce an SCESet object
make_mean_pheno <- function(sce) {
    exprs_df <- bind_cols(data_frame(line = sce$donor_long_id),
                          as_data_frame(t(norm_exprs(sce))))
    exprs_df_long <- tidyr::gather(exprs_df, key = "gene", value = "exprs",
                                   -line) %>% group_by(gene, line)
    ave_exprs_by_line_df <- exprs_df_long %>%
        summarise(ave_exprs = mean(exprs)) %>%
        tidyr::spread(key = gene, value = ave_exprs)
    lines <- ave_exprs_by_line_df[["line"]]
    ave_exprs_by_line <- as.matrix(dplyr::select(ave_exprs_by_line_df, -line))
    rownames(ave_exprs_by_line) <- lines
    ave_exprs_by_line
}


#' @param sce an SCESet object
make_var_pheno <- function(sce) {
    exprs_df <- bind_cols(data_frame(line = sce$donor_long_id),
                          as_data_frame(t(norm_exprs(sce))))
    exprs_df_long <- tidyr::gather(exprs_df, key = "gene", value = "exprs",
                                   -line) %>% group_by(gene, line)
    ave_exprs_by_line_df <- exprs_df_long %>%
        summarise(ave_exprs = var(exprs)) %>%
        tidyr::spread(key = gene, value = ave_exprs)
    lines <- ave_exprs_by_line_df[["line"]]
    ave_exprs_by_line <- as.matrix(dplyr::select(ave_exprs_by_line_df, -line))
    rownames(ave_exprs_by_line) <- lines
    ave_exprs_by_line
}


#' @param sce an SCESet object
make_alpha_pheno <- function(sce) {
    exprs_df <- bind_cols(data_frame(line = sce$donor_long_id),
                          as_data_frame(t(norm_exprs(sce))))
    exprs_df_long <- tidyr::gather(exprs_df, key = "gene", value = "exprs",
                                   -line) %>% group_by(gene, line)
    mean_var_by_line_df <- exprs_df_long %>%
        summarise(mean_exprs = mean(exprs), var_exprs = var(exprs))
    loess_predicts <- mean_var_by_line_df %>%
        ungroup %>% group_by(line) %>%
        do(loess_fit = predict(loess(var_exprs~ mean_exprs, data = .), data = .))
    ypred_list <- loess_predicts[["loess_fit"]]
    names(ypred_list) <- loess_predicts[["line"]]
    ypred <- as_data_frame(ypred_list)
    ypred[["gene"]] <- mean_var_by_line_df %>% ungroup %>%
        group_by(line) %>% dplyr::select(gene) %>% .[["gene"]] %>% unique
    ypred_long <- tidyr::gather(ypred, key = "line",
                                value = "loess_fitted_value", -gene)
    alpha_df <- left_join(mean_var_by_line_df, ypred_long,
                          by = c("gene", "line"))
    alpha_df <- dplyr::mutate(alpha_df, alpha = var_exprs- loess_fitted_value)
    alpha_by_line_df <- alpha_df %>%
        dplyr::select(-mean_exprs, -var_exprs, -loess_fitted_value) %>%
                                tidyr::spread(key = gene, value = alpha)
    lines <- alpha_by_line_df[["line"]]
    alpha_by_line_mat <- as.matrix(dplyr::select(alpha_by_line_df, -line))
    rownames(alpha_by_line_mat) <- lines
    alpha_by_line_mat
}


#' @param pheno_mat matrix of phenotype values with rows as samples and columns
#' as genes
normalise_pheno <- function(pheno_mat, method = "standardise") {
    rnames <- rownames(pheno_mat)
    cnames <- colnames(pheno_mat)
    method = match.arg(method, c("standardise"))
    ## plot hist before normalisation
    ## to_plot <- tidyr::gather(as_data_frame(pheno_mat[, 1:16]),
    ##                          key = "gene", value = "exprs")
    ## ggplot(to_plot, aes(x = exprs)) +
    ##     geom_histogram() +
    ##     ggtitle("Distributions for first 16 genes before standardization") +
    ##     facet_wrap(~gene, ncol = 4)
    colramp_samp <- colorRampPalette(
        c("dodgerblue4","gray50", "firebrick4"))(nrow(pheno_mat))
    colramp_gene <- colorRampPalette(
        c("dodgerblue4","gray50", "firebrick4"))(1000)
    ## plot densities
    ## par(mfcol = c(1, 2))
    ## plot(density(pheno_mat[1,]), col = colramp[1], lwd = 3, ylim = c(0, 1.0))
    ## for (i in 2:nrow(pheno_mat) {
    ##     lines(density(pheno_mat[i,]),lwd=3,col=colramp[i])
    ## }
    ## plot(density(pheno_mat[,1]), col=colramp[1], lwd=3, ylim=c(0,2.0))
    ## for (i in 2:1000) {
    ##     lines(density(pheno_mat[,i]),lwd=3,col=colramp[i])
    ## }
    ## normalise phenotype matrix
    norm_pheno_mat <- switch(
        method,
        standardise = double_standardised_quantile_norm_pheno(pheno_mat)
    )
    ## plot densities
    ## par(mfcol = c(1, 2))
    ## plot(density(norm_pheno_mat[1,]), col = colramp[1], lwd = 3, ylim = c(0, 1.0))
    ## for (i in 2:nrow(norm_pheno_mat) {
    ##     lines(density(norm_pheno_mat[i,]),lwd=3,col=colramp[i])
    ## }
    ## plot(density(norm_pheno_mat[,1]), col=colramp[1], lwd=3, ylim=c(0,2.0))
    ## for (i in 2:1000) {
    ##     lines(density(norm_pheno_mat[,i]),lwd=3,col=colramp[i])
    ## }
    ## normalise phenotype matrix
    rownames(norm_pheno_mat) <- rnames
    colnames(norm_pheno_mat) <- cnames
    norm_pheno_mat
}


double_standardised_quantile_norm_pheno <- function(pheno_mat) {
    norm_edata <- t(normalize.quantiles(t(pheno_mat)))
    double_norm_edata <- normalize.quantiles.use.target(
        norm_edata, target = qnorm(seq(1 / nrow(norm_edata),
        (1 - 1 / nrow(norm_edata)), length.out = nrow(norm_edata))))
    double_norm_edata
}


make_design_for_lmm <- function(pheno_mat, sample_meta, npcs_line = 10) {
    pca <- prcomp(pheno_mat)
    pca_df <- bind_cols(data_frame(line = rownames(pheno_mat)),
                                   as_data_frame(pca$x))
    mm <- match(pca_df[["line"]], sample_meta[["name"]])
    pca_df <- bind_cols(pca_df, sample_meta[mm,])
    if (npcs_line >=1) {
        xnam <- paste0("PC", 1:min(npcs_line, ncol(pca$x)))
        fmla <- as.formula(paste("~ ", paste(xnam, collapse= "+")))
        design <- model.matrix(fmla, data = pca_df)
    } else {
        design <- model.matrix(~1, data = pca_df)
    }
    rownames(design) <- pca_df[["line"]]
    design
}


write_output <- function(pheno_mat, gene_info, design, output_prefix) {
    write.table(t(pheno_mat), file = paste0(output_prefix, "exprs.tsv"),
                sep = "\t", quote = TRUE, col.names = NA)
    write.table(design, file = paste0(output_prefix, "covs.tsv"),
                 sep = "\t", quote = TRUE, col.names = NA)
    write_tsv(gene_info, path = paste0(output_prefix, "annos.tsv"))
    write_tsv(data_frame(genotype_individual_id = rownames(pheno_mat),
                         phenotype_sample_id = rownames(pheno_mat)),
              path = paste0(output_prefix, "samples.tsv"))
    limma_fit <- limma::lmFit(object = new("EList", list(E = t(pheno_mat))),
                              design = design, method = "ls")
    pheno_resids <- limma::residuals.MArrayLM(limma_fit, t(pheno_mat))
    write.table(pheno_resids, file = paste0(output_prefix, "exprs_resids.tsv"),
                sep = "\t", quote = TRUE)
}

  
#' Main function for running the script to make QTL phenotypes
main <- function(sce_file, output_prefix, pheno_type = "mean", npcs_cell = 0,
                 npcs_line = 10) {
    cat("Read SCESet....\n")
    pheno_type <- match.arg(pheno_type, c("mean", "var", "alpha"))
    sce <- readRDS(file = sce_file)
    featureNames(sce) <- fData(sce)$ensembl_gene_id
    sample_meta <- read_tsv("/nfs/research2/hipsci/drop/hip-drop/tracked/sample_meta_data/hipsci.qc1_sample_info.20160926.tsv", col_types = cols())
    ## filter donors based on cell number
    sce <- filter_donors(sce)
    ## regress out cell-level factors
    sce <- regress_cell_factors(sce, npcs = npcs_cell)
    cat("Make phenotype matrix....\n")
    ## add switch for different pheno types
    pheno_mat <- switch(pheno_type,
                        mean = make_mean_pheno(sce),
                        var = make_var_pheno(sce),
                        alpha = make_alpha_pheno(sce)
                        )
    cat("Normalise phenotype matrix....\n")
    norm_pheno_mat <- normalise_pheno(pheno_mat)
    cat("Make design matrix....\n")
    design_eqtl <- make_design_for_lmm(norm_pheno_mat, sample_meta, npcs_line)
    gene_info <- data_frame(
        feature_id = fData(sce)$ensembl_gene_id,
        chromosome = fData(sce)$chromosome_name,
        start = fData(sce)$start_position,
        end = fData(sce)$end_position,
        ensembl_gene_id = fData(sce)$ensembl_gene_id,
        strand = fData(sce)$strand,
        gene_name = fData(sce)$hgnc_symbol
    )
    cat("Write output to file(s)....\n")
    write_output(norm_pheno_mat, gene_info, design_eqtl, output_prefix)
    cat("Done!\n")
}

## Get command line options
opt <- docopt::docopt(doc, version = "version 0.0.1\n")
if (is.null(opt$npcs_cell))
    opt$npcs_cell <- 0
if (is.null(opt$npcs_line))
    opt$npcs_line <- 10

message("input SCESet: ", opt$input_file, "\n")
message("output prefix: ", opt$output_prefix, "\n")
message("N PCs at cell level: ", opt$npcs_cell, "\n")
message("N PCs at line level: ", opt$npcs_line, "\n")

## Run main function
main(opt$input_file, opt$output_prefix, opt$pheno_type,
     as.integer(opt$npcs_cell), as.integer(opt$npcs_line))

## "Rscript ../src/R/make_qtl_phenotypes.R -i ../data_processed/merged/sceset_merged_qc_hvg_day0.rds -o day0_mean_ -p mean"
## opt <- list()
## opt$input_file <- "../data_processed/merged/sceset_merged_qc_hvg_day0.rds"
## opt$output_prefix <- "day0_mean_"
## opt$pheno_type <- "mean"


#' @param sce an SCESet object
## make_aggcount_pheno <- function(sce) {

##     day0_counts_df <- bind_cols(data_frame(line = sce$donor_long_id), as_data_frame(t(counts(sce))))
##     day0_counts_df_long <- tidyr::gather(day0_counts_df, key = "gene", value = "exprs", -line) %>% group_by(gene, line)

##     head(day0_counts_df_long)
##     counts_by_line_df <- day0_counts_df_long %>% summarise(ave_exprs = sum(exprs)) %>% tidyr::spread(key = gene, value = ave_exprs)

##     lines <- counts_by_line_df[["line"]]
##     counts_by_line_mat <- as.matrix(dplyr::select(counts_by_line_df, -line))

##     dim(counts_by_line_mat)

##     counts_by_line_mat[1:6, 1:8]

##     summary(rowSums(counts_by_line_mat) * 1e-06)

##     options(repr.plot.width=9, repr.plot.height=4)
##     data_frame(x = rowSums(counts_by_line_mat)) %>%
##         ggplot(aes(x = x * 1e-06)) +
##         geom_histogram() + xlab("Total counts (millions)")

##     library(edgeR)

##     dge <- DGEList(counts = t(counts_by_line_mat), samples = data.frame(line = lines))
##     dge <- calcNormFactors(dge)

##     hist(dge$samples$norm.factors, main = "Histogram of normalization factors")

##     norm_cpm_by_line <- t(cpm(dge, log = TRUE))

##     head(norm_cpm_by_line)
## }



