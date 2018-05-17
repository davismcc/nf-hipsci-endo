'
Compile and automated report from an .RData file with an existing template

Usage:
combine_ase_results_run.R --input_dir <DIR_IN> --output_file <FILE_OUT> [-hv]

Options:
-i --input_dir <DIR_IN>      input directory containing CSV output files from
                             ASEReadCounter
-o --output_file <FILE_OUT>  specify output HTML file
-h --high_thresh             should high threshold ASE results be used?
-h --help                    show this
-v --version                 print version and stop

This program reads in and combines ASE results CSV files produced by the GATK
ASEReadCounter function.

This program requires the R packages "docopt", "tibble", and "dplyr" (CRAN) and
"Biobase" and "SummarizedExperiment" (Bioconductor).

The default is to use the more permissive "low_thresh" ASE calls, but instead
"high_thresh" calls can be used if desired.

Example arguments:
--input_dir "data_raw/scrnaseq"
--output_file "data_processed/ase/ase_counts.low_thresh.aggregated.rds"

Davis McCarthy
October 2016
' -> doc

library(tibble)
library(dplyr)
library(data.table)
library(Biobase)
library(GenomicRanges)
library(SummarizedExperiment)
library(Matrix)

saveAssayData <- function(df_list, variant_tbl, output_file) {
    output_file_prefix <- gsub(".rds", "", output_file)
    samplenames <- names(df_list)
    ## initialise matrices for assayData
    assay_slots <- c("refCount", "altCount", "totalCount",
                    "lowMAPQDepth", "lowBaseQDepth", "rawDepth", "otherBases",
                    "improperPairs", "refAltRatio")

    for (ass in assay_slots) {
        cat(paste("...compiling results for", ass, "\n"))
        out_file <- paste0(output_file_prefix, "_", ass, ".rds")
        ass_mat <- matrix(0,
            nrow = nrow(variant_tbl), ncol = length(df_list),
            dimnames = list(variant_tbl[["variantID"]], samplenames))
        ## incorporate data across samples into matrices
        pb <- txtProgressBar(min = 0, max = length(df_list), style = 3)
        counter <- 0
        for (this_sample in names(df_list)) {
            counter <- counter + 1
            vars <- df_list[[this_sample]][["variantID"]]
            ass_mat[vars, this_sample] <- df_list[[this_sample]][[ass]]
            setTxtProgressBar(pb, counter)
        }
        close(pb)
        saveRDS(ass_mat, out_file)
    }
    ad <- list(refAltRatio = ass_mat)
    ad
}

readASEResults <- function(ase_files) {
    df_list <- list()
    pb <- txtProgressBar(min = 0, max = length(ase_files), style = 3)
    counter <- 0
    for (this_file in ase_files) {
        counter <- counter + 1
        df_list[[this_file]] <- readr::read_tsv(this_file,
                                                col_types = "ciccciiiiiiii")
        df_list[[this_file]][["refAltRatio"]] <-
            (df_list[[this_file]][["refCount"]] /
                df_list[[this_file]][["totalCount"]])
        setTxtProgressBar(pb, counter)
    }
    close(pb)
    df_list
}

compileVariantInfo <- function(df_list) {
    cat("Compiling variant information into one table.\n")
    variant_tbl <- df_list[[1]][, 1:5]
    variant_dt1 <- data.table::data.table(variant_tbl)
    pb <- txtProgressBar(min = 1, max = (length(df_list) - 2), style = 3)
    counter <- 0
    ### first half of the data frames
    nhalf <- floor(length(df_list) / 2)
    for (i in 2:nhalf) {
        counter <- counter + 1
        this_dt <- data.table::data.table(df_list[[i]][, 1:5])
        variant_dt1 <- merge(variant_dt1, this_dt, all = TRUE,
            by = c("contig", "position", "variantID", "refAllele", "altAllele"))
        # variant_tbl <- full_join(variant_tbl, df_list[[i]][, 1:5],
        #                          by = c("contig", "position", "variantID",
        #                                 "refAllele", "altAllele"))
        setTxtProgressBar(pb, counter)
    }
    ### second half of the data frames
    variant_dt2 <- data.table::data.table(df_list[[nhalf + 1]][, 1:5])
    for (i in (nhalf + 2):length(df_list)) {
        counter <- counter + 1
        this_dt <- data.table::data.table(df_list[[i]][, 1:5])
        variant_dt1 <- merge(variant_dt1, this_dt, all = TRUE,
            by = c("contig", "position", "variantID", "refAllele", "altAllele"))
        setTxtProgressBar(pb, counter)
    }
    ### merge the two big data tables
    variant_dt <- merge(variant_dt1, variant_dt2, all = TRUE,
        by = c("contig", "position", "variantID", "refAllele", "altAllele"))
    close(pb)
    variant_dt
}


## Define main function
main <- function(input_dir, output_file, high_thresh = FALSE) {
    cat("Reading in ASE results files in ", input_dir, "\n")
    ## define ASE results files
    runs <- dir(input_dir)
    if (high_thresh)
        ase_dirs <- file.path(input_dir, runs, "ase", "high_thresh")
    else
        ase_dirs <- file.path(input_dir, runs, "ase", "low_thresh")
    ase_files <- c()
    for (this_dir in ase_dirs) {
        these_files <- list.files(this_dir, pattern = ".ase.",
                                  full.names = TRUE, recursive = TRUE)
        ase_files <- c(ase_files, these_files)
    }
    ## read in all ASE results data
    df_list <- readASEResults(ase_files)
    ## get variant information into one table
    variant_dt <- compileVariantInfo(df_list)
    ## output information about number of samples and variants processed
    variant_tbl <- tibble::as_data_frame(variant_dt)
    oo <- order(variant_tbl[["contig"]], variant_tbl[["position"]])
    variant_tbl <- variant_tbl[oo, ]
    cat("ASE results read in from ", length(df_list), " samples.\n")
    cat("Results obtained for ", nrow(variant_tbl), " variants.\n")
    ## build SummarizedExperiment object
    samplenames <- gsub(".ase.csv", "", gsub(".*/ase/", "", names(df_list)))
    names(df_list) <- samplenames
    cat("Building SummarizedExperiment object...\n")
    ad <- saveAssayData(df_list, variant_tbl, output_file)
    gr <- GRanges(seqnames = variant_tbl[["contig"]],
                  ranges = IRanges(start = variant_tbl[["position"]],
                  width = 1),
                  names = variant_tbl[["variantID"]])
    mcols(gr) <- as.data.frame(variant_tbl)
    ase_set <- SummarizedExperiment(ad, rowRanges = gr)
    cat("Writing combined ASE counts to ", output_file, "\n")
    saveRDS(ase_set, output_file)
    ase_set
}

## Get command line options
#opt <- docopt::docopt(doc, version = "version 0.0.1\n")

## test code
opt <- list()
opt$input_dir <- "data_raw/scrnaseq"
opt$output_file <- "data_processed/ase/ase_counts.low_thresh.aggregated.rds"
ase <- main(opt$input_dir, opt$output_file)

## Run main function
#main(opt$input_dir, opt$output_file)
