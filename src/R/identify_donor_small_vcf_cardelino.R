"
Identify donor from a cell's filtered VCF file

Usage:
identify_donor_small_vcf_cardelino.R --input_file <DATA_IN> --output_prefix <PRFX_OUT> --donor_lines <DONORS> [--donor_vcf <DONOR_VCF> -hv]

Options:
-i --input_file <DATA_IN>      input file
-o --output_prefix <PRFX_OUT>  prefix for output files; .RData, .feather and .csv files produced
-l --donor_lines <DONORS>      string providing donor line IDs separated by ';' (e.g. 'babk_2;wuye_2;podx3')
-d --donor_vcf <DONOR_VCF>     VCF file containing the genotype data for all donors [default: /nfs/research2/stegle/projects/hipsci/data/genotypes/imputed/REL-2014-11_SS/hipsci.wec.gtarray.HumanCoreExome-12_v1_0.REL-2014-11.imputed_phased.20151005.genotypes.gdid.mac1.recode.WGS.maf0.1.vcf.gz]
-h --help                      show this
-v --version                   print version and stop

This program does compares called genotypes from RNA-seq reads for a cell to
genotype data for all HipSci donors to find the best match in terms of
(weighted) average genotypic correlation across called variants.

The program returns results as a .csv file.

This program requires the R packages 'VariantAnnotation' and 'snpStats'
(Bioconductor) and 'cardelino' (github.com/davismcc/cardelino).

Davis McCarthy
March 2018
" -> doc

## Script to identify HIPSCI donor for a cell's VCF file
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(snpStats))
suppressPackageStartupMessages(library(cardelino))


init_output_df <- function(input_vcf, probs, these_donors, short_names) {
    ## define output data frame
    output_df <- data.frame(
        sample_id = gsub(".filtered.*", "", basename(input_vcf)),
        donor = as.character(names(probs)),
        nvars_called = 0,
        nvars_ovlap_hipsci = 0,
        nvars_used = 0,
        post_prob = probs,
        used_in_expt = (short_names %in% these_donors),
        stringsAsFactors = FALSE)
    output_df
}

full_output_df <- function(input_vcf, vcf_sample, vcf_hipsci, sm_sample,
                           probs, these_donors, short_names) {
    ## define output data frame
    output_df <- data.frame(
        sample_id = gsub(".filtered.*", "", basename(input_vcf)),
        donor = as.character(names(probs)),
        nvars_called = length(vcf_sample),
        nvars_ovlap_hipsci = length(vcf_hipsci),
        nvars_used = ncol(sm_sample$genotypes),
        post_prob = probs,
        used_in_expt = (short_names %in% these_donors), 
        stringsAsFactors = FALSE)
    output_df
}

read_sample_vcf <- function(input_vcf) {
    ## Read in VCF from this sample
    message("Reading sample VCF\n")
    vcf_sample <- readVcf(input_vcf, "GRCh37")
    seqlevelsStyle(vcf_sample) <- "UCSC"
    vcf_sample_filt <- vcf_sample[isSNV(vcf_sample)]
    if (length(vcf_sample_filt) > 0) {
        new_snp_names <- paste0("snp_",
                                gsub("chr", "",
                                     gsub(":", "_",
                                          gsub("_[ATCG]/[ATCG]", "",
                                               names(vcf_sample_filt)))))    
        names(vcf_sample_filt) <- new_snp_names
    }
    vcf_sample_filt
}

read_donor_vcf <- function(donor_vcf) {
    message("Reading donors VCF\n")
    vcf_donor <- readVcf(donor_vcf, "GRCh37")
    seqlevelsStyle(vcf_donor) <- "UCSC"
    isSNV_idx <- isSNV(vcf_donor)
    vcf_donor <- vcf_donor[, grep("HPS.*pf-[a-z]+$", colnames(vcf_donor))]
    colnames(vcf_donor) <- gsub("HPS.*-", "", gsub("_[0-9]+", "", colnames(vcf_donor)))
    list(vcf = vcf_donor, isSNV_idx = isSNV_idx)
}

get_snp_matrices <- function(vcf_sample, vcf_donor) {
    ## get snp matrices
    sm_donor <- genotypeToSnpMatrix(
        geno(vcf_donor, "GT"), ref = ref(vcf_donor), alt = alt(vcf_donor))
    ## filter sample VCF to those variants found in donor VCF
    message("Filtering sample VCF\n")
    vcf_sample <- sortSeqlevels(vcf_sample, X.is.sexchrom = TRUE)
    slengths_sample <- seqlengths(vcf_sample)
    vcf_donor <- sortSeqlevels(vcf_donor, X.is.sexchrom = TRUE)
    seqlengths(vcf_donor) <- slengths_sample[seqlevels(vcf_donor)]
    ovlap <- findOverlaps(vcf_sample, vcf_donor)
    if (length(ovlap) < 1L) {
        message("No common variants overlapping in sample VCF and Donor VCF\n")
        return(list(stop_program = TRUE))
    } else {
        vcf_sample2 <- vcf_sample[queryHits(ovlap)]
        vcf_donor <- vcf_donor[subjectHits(ovlap)]
        match_alleles <- unlist(ref(vcf_sample2) == ref(vcf_donor) & alt(vcf_sample2) == alt(vcf_donor))
        if (sum(match_alleles) < 1L)
            stop("No variants with matching alleles in sample and donor VCFs")
        else {
            vcf_sample2 <- vcf_sample2[match_alleles]
            vcf_donor <- vcf_donor[match_alleles]
        }
        sm_sample <- genotypeToSnpMatrix(
                geno(vcf_sample2, "GT"), ref = ref(vcf_sample2),
                alt = alt(vcf_sample2))
        sm_sample_REF <- matrix(sapply(geno(vcf_sample2, "AD"), function(x) x[[1]]), ncol = 1)
        sm_sample_ALT <- matrix(sapply(geno(vcf_sample2, "AD"), function(x) x[[2]]), ncol = 1)
        sm_sample_DEP <- matrix(sm_sample_REF + sm_sample_ALT, ncol = 1)
        na_sample <- is.na(sm_sample_DEP)
        if (sum(!na_sample) < 1L) {
            message("No common variants with non-missing genotypes overlapping in sample VCF and Donor VCF\n")
            return(list(stop_program = TRUE))
        } else {
            sm_donor <- genotypeToSnpMatrix(
                geno(vcf_donor, "GT"), ref = ref(vcf_donor),
                alt = alt(vcf_donor))
            donor_geno_mat <- matrix(as.numeric(as(sm_donor$genotypes, "numeric") > 0), nrow = nrow(sm_donor$genotypes))
            rownames(donor_geno_mat) <- rownames(sm_donor$genotypes)
            colnames(donor_geno_mat) <- colnames(sm_donor$genotypes)
            donor_geno_mat <- t(donor_geno_mat[order(rownames(donor_geno_mat)),])
            message("Doing donor assignment using ", length(vcf_sample2),
                    " variants\n")
            return(list(sm_sample_REF = sm_sample_REF, A = sm_sample_ALT, D = sm_sample_DEP, 
                        sm_sample = sm_sample, stop_program = FALSE,
                       C = donor_geno_mat))
        }
    }
}


main <- function(input_vcf, output_prefix, donor_lines, donor_vcf, fasta_idx) {
    ## define samples in donor VCF
    hdr_donor <- scanVcfHeader(donor_vcf)
    donors <- samples(hdr_donor)[grep("HPS.*pf-[a-z]+$", samples(hdr_donor))]
    donors <- gsub("HPS.*-", "", gsub("_[0-9]+", "", donors))
    ## define donors and names
    short_names <- donors
    these_donors <- gsub("HPS.*-", "", gsub("_[0-9]+", "", strsplit(donor_lines, ";")[[1]]))
    these_donors_idx <- grepl(gsub(";", "|", donor_lines), donors)
    probs <- rep(0, length(donors))
    names(probs) <- donors
    ## define output data frame
    output_df <- init_output_df(input_vcf, probs, these_donors, short_names)
    ## Read in VCF from this sample
    vcf_sample <- read_sample_vcf(input_vcf)
    if (length(vcf_sample) < 1) {
        write.csv(output_df, file = paste0(output_prefix, ".csv"),
                  row.names = FALSE)
        message("No variants in sample VCF after filtering.\n")
        return("Done.")
    }
    output_df$nvars_called <- length(vcf_sample)
    message("...read ", length(vcf_sample), " variants from sample VCF\n")
    ## Read in Donor VCF
    donor_data <- read_donor_vcf(donor_vcf)
    isSNV_idx <- donor_data$isSNV_idx
    vcf_donor <- donor_data$vcf
    message("...read ", length(vcf_donor), " variants from sample VCF\n")
    if (!any(isSNV_idx)) {
        write.csv(output_df, file = paste0(output_prefix, ".csv"),
                  row.names = FALSE)
        message("No single-nucleotide variants overlapping in sample VCF and Donor VCF\n")
        return("Done.")
    }
    if (!identical(colnames(vcf_donor), short_names))
        stop("Sample names in donor VCF do not match those read from header.")
    vcf_donor <- vcf_donor[isSNV_idx]
    ## get snp matrices
    snpmat_list <- get_snp_matrices(vcf_sample, vcf_donor)
    if (snpmat_list$stop_program) {
        write.csv(output_df, file = paste0(output_prefix, ".csv"),
                  row.names = FALSE)
        message("No common variants overlapping in sample VCF and Donor VCF\n")
        return("Done.")
    }
    assign <- cell_assign_Gibbs(
        A = snpmat_list$A, D = snpmat_list$D, C = snpmat_list$C,
        Psi = rep(1 / ncol(snpmat_list$C), ncol(snpmat_list$C)),
        model = "Binomial")
    probs <- as.vector(assign$prob)
    names(probs) <- colnames(assign$prob)
    output_df <- full_output_df(input_vcf, vcf_sample, vcf_donor, snpmat_list$sm_sample,
                                probs, these_donors, short_names)
    output_df[["n_total_reads"]] <- sum(snpmat_list$D, na.rm = TRUE)
    output_df[["n_alt_reads"]] <- sum(snpmat_list$A, na.rm = TRUE)
    ## write output to file
    message("Writing output to CSV\n")
    write.csv(output_df, file = paste0(output_prefix, ".csv"),
              row.names = FALSE)
    return("Done.")
}

## Get command line options
opt <- docopt::docopt(doc, version = "version 0.0.1\n")

message("working directory: ", getwd(), "\n")
message("input vcf: ", opt$input_file, "\n")
message("output prefix: ", opt$output_prefix, "\n")
message("donor vcf: ", opt$donor_vcf, "\n")
message("lines used: ", opt$donor_lines, "\n")

## Run main function
main(opt$input_file, opt$output_prefix, opt$donor_lines, opt$donor_vcf)

## # params for testing
## opt <-  list()
## opt[["input_file"]] <- "data_raw/scrnaseq/run_21999/vcf/21999_1#56.filtered.vcf.gz"
## opt[["output_prefix"]] <- "data_raw/scrnaseq/run_21999/donor_id/tmp.donor_id"
## opt[["donor_vcf"]] <- "data_raw/scrnaseq/run_21999/vcf/21999_1#56.filtered.hipsci.overlap.vcf.gz"
## opt[["donor_lines"]] <- "fasu_2;kegd_2;zerv_8;zoio_2;xojn_3;fuai_1;eevy_7;oaqd_3;paab_4;sita_1;toss_3;zoio_2;heth_1;jogf_2;pelm_3;vass_1;wibj_2;zapk_3"

