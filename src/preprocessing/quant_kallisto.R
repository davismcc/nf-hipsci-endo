'
Quantify expression with kallisto and return an .RData file with SCESets

Usage:
  quant_kallisto.R [--run_on_cluster | --read_results] [--input_dir <DIRECTORY_IN> --output_dir <DIRECTORY_OUT> --kallisto_idx <IDX_FILE> --scesets_out <RDATA_FILE> --n_cores <NUM> --biomart <HOST> --single_end --no_bias_correction --kallisto_cmd <KALLISTO> -rqhv]

Options:
  --run_on_cluster                 tell the program to run kallisto jobs on the cluster
  --read_results                   tell the program to read kallisto results into an SCESet and save
  -i --input_dir <DIRECTORY_IN>    input directory [default: ./test_in]
  -o --output_dir <DIRECTORY_OUT>  specify output directory for kallisto results [default: ./test_out]
  -k --kallisto_idx <IDX_FILE>     specify kallisto index file [default: /nfs/research2/stegle/datasets/references/human/GRCh37/Homo_sapiens.GRCh37.rel75.cdna.all.ERCC92.kallisto_idx]
  -s --scesets_out <RDATA_FILE>    specify .RData file to save generated SCESets to [default: ./scesets_kallisto_out.RData]
  -n --n_cores <NUM>               number of cores to use for each kallisto run [default: 4]
  -b --biomart <HOST>              provide the URL for the Biomart version to use; default is current Ensembl release [default: feb2014.archive.ensembl.org]
  --kallisto_cmd                   provide the full path for the kallisto command if required [default: /nfs/software/stegle/kallisto_linux-v0.43.0/kallisto]
  --single_end                     tell kallisto that the the reads are single end
  --no_bias_correction             tell kallisto not to apply its bias correction
  --kallisto_cmd                   the command to call to run kallisto [default: /nfs/software/stegle/kallisto_linux-v0.43.0/kallisto]
  -r --recompute                   recompute previously computed kallisto data in the same location
  -q --quiet                       print less text
  -h --help                        show this
  -v --version                     print version and stop

This program does two distinct things:
1) if the "--run_on_cluster" option is given then it runs kallisto on the FASTQ
files in the input directory, sending the jobs off to the cluster;

2) if the "--read_results" option is given then it reads the kallisto results
into an SCESet object (from the scater package) and saves an RData file with
two objects, one containing gene-level expression data and one containing
transcript-level expression data.

Provide the input directory containing FASTQ files and the output directory
where kallisto results are to be written. kallisto will be used to quantify
transcript expression using the provided kallisto index file.

The R package scater is used to read the kallisto results into an SCESet object,
a convenient data format for single-cell expression data. Transcript and gene
annotation information is obtained from Biomart and two SCESet objects, one
containing transcript-level expression values and one containing gene-level
expression values are saved in an RData file named according to the optional
argument.

The default kallisto index used (for backwards compatibility with previous
results; Ensembl release 75, February 2015) is:
/nfs/research2/stegle/datasets/references/human/GRCh37/Homo_sapiens.GRCh37.rel75.cdna.all.ERCC92.kallisto_idx

A kallisto index produced using Ensembl version 83 (December 2015) is here:
/nfs/research2/stegle/datasets/references/human/GRCh38/Homo_sapiens.GRCh38.rel83.cdna.all.ERCC92.kallisto_idx

The kallisto command used is:
/nfs/software/stegle/kallisto_linux-v0.43.0/kallisto
/homes/davis/davis/src/kallisto_linux-v0.42.4/kallisto

A host for Biomart can be specified to use archived Ensembl versions,
e.g. "feb2014.archive.ensembl.org". To use the current version of Ensembl use
"www.ensembl.org".

This program requires the R packages "docopt" (CRAN) and "scater" (github.com/davismcc/scater).

Davis McCarthy
December 2015
' -> doc

## returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

## Define main function
main <- function(input_dir, output_dir, kallisto_idx, scesets_out, kallisto_cmd,
                 run_on_cluster = FALSE, read_results = FALSE, n_cores = 4,
                 single_end = FALSE, correct_bias = TRUE, host = NULL,
                 recompute = FALSE, verbose = TRUE) {
    ## Check that the call makes sense
    if ( !(run_on_cluster | read_results) )
        stop("Please specify either --run_on_cluster (to run kallisto jobs) or
             --read_results (to read kallisto results into SCESet objects)")

    ## Define some parameters
    #kallisto_cmd <- "/homes/davis/davis/src/kallisto_linux-v0.42.4/kallisto"

    ## Load scater
    library(scater)

    ## Define targets file
    targets_file <- file.path(input_dir, "targets_kallisto.txt")

    if ( file.exists(targets_file) ) {
        targets_df <- read.table(targets_file, stringsAsFactors = FALSE)
        if ( verbose )
            cat("Using existing targets file in input directory")
    } else {
        ## Define targets
        fastq_left <- dir(input_dir)[grep("_1.fastq$", dir(input_dir))]
        sample_names <- gsub("#", "_", gsub("_[12].fastq", "", fastq_left))
        fastq_left <-  file.path(input_dir, fastq_left)
        fastq_left <- trim(fastq_left)
        fastq_right <- gsub("_1.fastq", "_2.fastq", fastq_left)
        names_match <- identical(gsub("_[12].fastq", "", fastq_left),
                                 gsub("_[12].fastq", "", fastq_right))
        if ( !names_match )
            stop("Base names of FASTQ files do not match")
        targets_df <- data.frame(Sample = sample_names, FASTQ_1 = fastq_left,
                                 FASTQ_2 = fastq_right)
        if ( verbose ) {
            cat("Targets: \n")
            cat(targets_df)
        }
        write.table(targets_df, row.names = FALSE, file = targets_file,
                    quote = TRUE, sep = "\t")
    }

    ## Run kallisto on the fastq files
    kallisto_run <- runKallisto(targets_file,
                                 kallisto_idx, output_prefix = "out",
                                 n_cores = n_cores, single_end = single_end,
                                 correct_bias = correct_bias, verbose = verbose,
                                 dry_run = TRUE)

    ## Expected kallisto output files
    kall_files_out <- file.path(sapply(kallisto_run, function(x) x$output_dir),
                                "abundance.tsv")

    ## if argument is run_on_cluster then run cluster jobs in parallel
    if ( run_on_cluster ) {
        ## If asked to recompute, delete any existing output files in the same location
        if ( recompute ) {
            dirs_to_delete <- sapply(kallisto_run, function(x) x$output_dir)
            dirs_to_delete[file.exists(kall_files_out)]
            unlink(dirs_to_delete, recursive = TRUE)
        } else {
            if ( any(file.exists(kall_files_out)) )
                stop("Asked not to recompute results, but some results files already exist. Please clean directories or recompute with -r or --recompute flag")
        }

        ## specify cluster out folder (and create if it does not already exist)
        cluster_out_dir <- file.path(output_dir, 'cluster_out')
        if ( !dir.exists(cluster_out_dir) )
            dir.create(cluster_out_dir)
        ## build cluster command

        cluster_prefix <-  paste0('bsub -n ', n_cores)
        cluster_prefix <- paste0(cluster_prefix, ' -q research-rh6 -R "rusage[mem=8000,tmp=8000]" -M 8000 ')
        cluster_prefix <- paste0(cluster_prefix, '-o ', cluster_out_dir)
        ## Submit kallisto jobs to cluster

        for( i in seq_len(length(kallisto_run)) ) {
            cmd <- paste0(
                cluster_prefix, "_", i, "kallisto_run[[i]]$kallisto_call")
            cmd <- gsub(" kallisto ", paste("", kallisto_cmd), cmd)
            if ( verbose )
                print(cmd)
            system(cmd)
            Sys.sleep(1)
        }
        if ( verbose )
            cat("Submitted kallisto jobs to cluster")
    }

    ## Alternatively, read in kallisto results into SCESet objects and save
    if ( read_results ) {
        if( any(!file.exists(kall_files_out)) ) {
            stop("Not all of the expected output files exist. Check that the
                  correct kallisto results have been produced?")
        }

        ## Read kallisto results in R
        sce_kallisto_tx <- readKallistoResults(kallisto_run, read_h5 = TRUE)

        ## Check results
        if ( any(apply(exprs(sce_kallisto_tx), 2,
                       function(x) {any(is.na(x))})) ) {
            warning("The following samples have NA expression values")
            ww <- which(apply(exprs(sce_kallisto_tx), 2,
                              function(x) {any(is.na(x))}))
            cat(sampleNames(sce_kallisto_tx)[ww])
        }

        ## Get transcript feature annotations from biomaRt
        sce_kallisto_tx <- getBMFeatureAnnos(
            sce_kallisto_tx, filters = "ensembl_transcript_id",
            attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                           "hgnc_symbol", "chromosome_name", "transcript_biotype",
                           "transcript_start", "transcript_end", "transcript_count"),
            feature_symbol = "hgnc_symbol", feature_id = "ensembl_gene_id",
            dataset = "hsapiens_gene_ensembl", host = host)
        featureNames(sce_kallisto_tx) <- paste(
            featureNames(sce_kallisto_tx), fData(sce_kallisto_tx)$hgnc_symbol,
            sep = "_")

        ## Summarise expression at gene level
        sce_kallisto_gene <- summariseExprsAcrossFeatures(sce_kallisto_tx)

        ## Get gene feature annotations from biomaRt
        sce_kallisto_gene <- getBMFeatureAnnos(
            sce_kallisto_gene, filters = "ensembl_gene_id",
            attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                           "hgnc_symbol", "chromosome_name", "gene_biotype"),
            feature_symbol = "hgnc_symbol", feature_id = "ensembl_gene_id",
            dataset = "hsapiens_gene_ensembl", host = host)
        featureNames(sce_kallisto_gene) <- paste(
            featureNames(sce_kallisto_gene), fData(sce_kallisto_gene)$hgnc_symbol,
            sep = "_")

        ## Save SCESet objects
        if ( verbose ) {
            cat("Saving RData object to: ")
            cat(scesets_out)
        }
        save(sce_kallisto_tx, sce_kallisto_gene, file = scesets_out,
             compress = 'gzip', compression_level = 9)
        if ( verbose )
            cat(paste("Saved RData file in ", output_dir))
    }
    cat("This program has completed.")
}

## Get command line options
opt <- docopt::docopt(doc, version = "version 0.0.1\n")
n_cores <- as.numeric(opt$n_cores)
correct_bias <- !opt$no_bias_correction

## Run main function
main(opt$input_dir, opt$output_dir, opt$kallisto_idx, opt$scesets_out,
     opt$kallisto_cmd,
     run_on_cluster = opt$run_on_cluster, read_results = opt$run_on_cluster,
     n_cores = n_cores, single_end = opt$single_end,
     correct_bias = correct_bias, opt$biomart, opt$recompute,
     verbose = !opt$quiet)
