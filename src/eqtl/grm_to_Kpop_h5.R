"
Read GRM produced by PLINK and save to HDF5 format

Usage:
grm_to_Kpop_h5.R --input_file <DATA_IN> --output_file <FILE_OUT> [-hv]

Options:
-i --input_file <DATA_IN>      input file
-o --output_file <FILE_OUT>    name for output file
-h --help                      show this
-v --version                   print version and stop

Convenience program to make QTL analyses easier.

This program requires the R packages 'docopt', 'tidyverse' and 'rhdf5'
(Bioconductor).

Davis McCarthy
January 2017
" -> doc

library(rhdf5)

## Define main function
main <- function(input_file, output_file) {
    grm <- readr::read_tsv(input_file, col_names = FALSE)
    donors_df <- readr::read_tsv(gsub("gz", "id", input_file),
                                 col_names = FALSE, col_types = "cc")
    donors <- donors_df[[1]]
    rhdf5::h5createFile(file = output_file)
    rhdf5::h5write(donors, file = output_file, name = 'sampleID')
    rhdf5::h5write(as.matrix(grm), file = output_file, name = 'Kpop')
    cat("Kpop saved to ", output_file)
}

## Get command line options
opt <- docopt::docopt(doc, version = "version 0.0.1\n")

## Run main function
main(opt$input_file, opt$output_file)
