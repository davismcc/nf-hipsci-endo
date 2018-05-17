 '
Read expression values from Salmon and return .rds files with SCESets

Usage:
  read_salmon_to_scesets.R [--input_dir <DIRECTORY_IN> --output_prefix <OUT_PREFIX> --biomart <HOST> -rqhv]

Options:

  -i --input_dir <DIRECTORIES_IN>  input directory containing sub-directories with Salmon quantification results
  -o --output_prefix <OUT_PREFIX>  file path prefix for output files with merged results [default: ./test_out]
  -b --biomart <HOST>              provide the URL for the Biomart version to use; default is current Ensembl release [default: dec2015.archive.ensembl.org]
  -q --quiet                       print less text
  -h --help                        show this
  -v --version                     print version and stop

This program reads the Salmon results into an SCESet object (from the
scater package) and saves one .rds file with gene-level expression data and one
containing transcript-level expression data.

Provide the input directory containing subdirectories with Salmon
quantification results and an output file path prefix where the .rds files
containing the SCESet objects are to be written.

The R package scater is used to read the Salmon results into an SCESet object,
a convenient data format for single-cell expression data. Transcript and gene
annotation information is obtained from Biomart.

A host for Biomart can be specified to use archived Ensembl versions,
e.g. "dec2015.archive.ensembl.org". To use the current version of Ensembl use
"www.ensembl.org".

This program requires the R packages "docopt" (CRAN) and "scater" (Bioconductor).

Davis McCarthy
June 2016
' -> doc

main <- function(input_dir, output_prefix, host = NULL, verbose = TRUE) {
    ## Load scater
    library(scater)

    ## Define targets file
    input_directories <- list.dirs(input_dir, recursive = FALSE)

    ## Read kallisto results in R
    sce_tx <- readSalmonResults(
        samples = input_directories, directories = input_directories)
    if ( "exprs" %in% assayNames(sce_tx) ) {
        logcounts(sce_tx) <- assay(sce_tx, "exprs")
        assay(sce_tx, "exprs") <- NULL
    }

    ## Check results
    if ( any(apply(exprs(sce_tx), 2, function(x) {any(is.na(x))})) ) {
        warning("The following samples have NA expression values")
        ww <- which(apply(exprs(sce_tx), 2, function(x) {any(is.na(x))}))
        cat(colnames(sce_tx)[ww])
    }

    ## Get transcript feature annotations from biomaRt
    sce_tx <- getBMFeatureAnnos(
        sce_tx, filters = "ensembl_transcript_id",
        attributes = c("ensembl_transcript_id", "ensembl_gene_id",
        "hgnc_symbol", "chromosome_name", "transcript_biotype",
        "transcript_start", "transcript_end", "transcript_count"),
        feature_symbol = "hgnc_symbol", feature_id = "ensembl_gene_id",
        dataset = "hsapiens_gene_ensembl", host = host)

    ## tweak feature names to be more readable
    rownames(sce_tx) <- paste(
        rownames(sce_tx), rowData(sce_tx)$hgnc_symbol, sep = "_")

    ## Save SCESet object
    if ( verbose ) {
        cat("Saving transcript expression object to: ")
        cat(paste0(output_prefix, "_tx.rds"))
        cat("\n")
    }
    saveRDS(sce_tx, file = paste0(output_prefix, "_tx.rds"))

    ## If number of samples is large, save tx into multiple rds files
    ## each containing a max of 500 samples
    idx_list <- split(colnames(sce_tx), ceiling(seq_len(ncol(sce_tx)) / 500))
    tmp_file_list <- paste0(output_prefix, "_tmp_tx_", 1:length(idx_list), ".rds")
    for (i in seq_along(idx_list)) {
        saveRDS(sce_tx[, idx_list[[i]]], file = tmp_file_list[i])
    }
    cat("Removing transcript expression object...")
    rm(sce_tx)
        
    ## Summarise expression at gene level
    cat("Summarising transcript expression at gene level...\n")
    tmp_tx <- readRDS(tmp_file_list[1])
    sce_gene <- summariseExprsAcrossFeatures(tmp_tx)
    if (length(idx_list) >= 2){
        for (i in 2:length(idx_list)) {
            cat("File ", i, " of ", length(idx_list), "\n")
            tmp_tx <- readRDS(tmp_file_list[i])
            tmp_gene <- summariseExprsAcrossFeatures(tmp_tx)
            sce_gene <- cbind(sce_gene, tmp_gene)
        }
    }

    ## remove temporary files
    for (i in seq_along(idx_list)) {
        rm_cmd <- paste("rm", tmp_file_list[i])
        system(rm_cmd)
    }
    
    ## Get gene feature annotations from biomaRt
    sce_gene <- getBMFeatureAnnos(
        sce_gene, filters = "ensembl_gene_id",
        attributes = c("ensembl_transcript_id", "ensembl_gene_id",
        "hgnc_symbol", "chromosome_name", "start_position", "end_position",
        "strand", "gene_biotype"), feature_symbol = "hgnc_symbol",
        feature_id = "ensembl_gene_id", dataset = "hsapiens_gene_ensembl",
        host = host)

    ## tweak feature names to be more readable
    rownames(sce_gene) <- paste(
        rownames(sce_gene), rowData(sce_gene)$hgnc_symbol, sep = "_")

    ## calculate QC metrics
    ercc_genes <- grepl("^ERCC", rownames(sce_gene))
    mt_genes <- grepl("^MT-", rowData(sce_gene)$hgnc_symbol)
    sce_gene <- calculateQCMetrics(
        sce_gene, feature_controls = list(ERCC = ercc_genes, MT = mt_genes))

    ## save gene annotation data to a feather file
    ## annos <- rowData(sce_gene)
    ## annos <- as(annos, "data.frame")
    ## feather::write_feather(annos, paste0(output_prefix, ".feather"))

    if ( verbose ) {
        cat("Saving gene expression object to: ")
        cat(paste0(output_prefix, "_gene.rds"))
        cat("\n")
    }
    saveRDS(sce_gene, file = paste0(output_prefix, "_gene.rds"))

    cat("This program has completed.\n")
}

## Get command line options
opt <- docopt::docopt(doc, version = "version 0.0.1\n")

## Run main function
main(opt$input_dir, opt$output_prefix, opt$biomart, verbose = !opt$quiet)
