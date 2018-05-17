"
Conduct test for differences in ASE using a binomial GLMM

Usage:
run_ase_binomial_glmm.R --input_ase <ASE_FILE> --input_exprs <EXPRS_FILE> [-hv]

Options:
-a --input_ase <ASE_FILE>      input file containing ASE results as a
                               SummarizedExperiment object
-e --input_exprs <EXPRS_FILE>  input file containing expression values for QCed
                                cells in an SCESet object
-h --help                      show this
-v --version                   print version and stop

This function takes as input a file containing ASE results in a
SummarizedExperiment experiment object and a file containing an SCESet object
with expression values for QC'd cells.

Davis McCarthy
October 2016
" -> doc

library(SummarizedExperiment)
library(scater)
library(lme4)
library(BiocParallel)
library(tibble)
library(dplyr)
library(pryr)
library(clustermq)


fit_ase_glmm <- function(refalt_count, coldata) {
    dat4glmm <- cbind(refalt_count, as.data.frame(coldata))
    dat4glmm <- dat4glmm[!is.na(dat4glmm$refCount) &
                            !is.na(dat4glmm$altCount), ]
    out <- tryCatch({
        ## define data object for GLMM fit
        if (length(unique(dat4glmm$experiment)) < 2L) {
            ## only one environment, so don't fit random effect
            gm0 <- glm(cbind(altCount, refCount) ~ 1,
                        data = dat4glmm, family = binomial(link = "logit"))
            out <- c(0, AIC(gm0), BIC(gm0), logLik(gm0), deviance(gm0))
            ### test donor
            if (length(unique(dat4glmm$donor)) < 2L) {
                out <- c(out, rep(NA, 4))
                P_donor <- NA
            } else {
                gm1 <- glm(cbind(altCount, refCount) ~ donor,
                                data = dat4glmm, family = binomial(link = "logit"))
                out <- c(out, AIC(gm1), BIC(gm1), logLik(gm1), deviance(gm1))
                P_donor <- 1 - pchisq(gm0$deviance - gm1$deviance,
                                gm0$df.residual - gm1$df.residual)
            }
            ### test differentiation population
            if (length(unique(dat4glmm$diff_population)) < 2L) {
                out <- c(out, rep(NA, 4))
                P_donor <- NA
            } else {
                gm1 <- glm(cbind(altCount, refCount) ~ diff_population,
                                data = dat4glmm, family = binomial(link = "logit"))
                out <- c(out, AIC(gm1), BIC(gm1), logLik(gm1), deviance(gm1))
                P_diffpop <- 1 - pchisq(gm0$deviance - gm1$deviance,
                                gm0$df.residual - gm1$df.residual)
            }
            ### test ouija_pseudotime
            gm1 <- glm(cbind(altCount, refCount) ~ ouija_pseudotime,
                        data = dat4glmm, family = binomial(link = "logit"))
            out <- c(out, AIC(gm1), BIC(gm1), logLik(gm1), deviance(gm1))
            P_pt <- 1 - pchisq(gm0$deviance - gm1$deviance,
                            gm0$df.residual - gm1$df.residual)
            ### add P-values to output vector
            out <- c(P_donor, P_diffpop, P_pt, out)
        } else {
            ## fit GLMM with glmer
            gm0 <- glmer(cbind(altCount, refCount) ~ (1 | experiment),
                            data = dat4glmm, family = binomial(link = "logit"))
            out <- c(1, AIC(gm0), BIC(gm0), logLik(gm0), deviance(gm0))
            ### test donor
            if (length(unique(dat4glmm$donor)) < 2L) {
                out <- c(out, rep(NA, 4))
                P_donor <- NA
            } else {
                gm1 <- glmer(cbind(altCount, refCount) ~ donor + (1 | experiment),
                        data = dat4glmm, family = binomial(link = "logit"))
                out <- c(out, AIC(gm1), BIC(gm1), logLik(gm1), deviance(gm1))
                P_donor <- anova(gm0, gm1)$P[2]
            }
            ### test differentiation population
            if (length(unique(dat4glmm$diff_population)) < 2L) {
                out <- c(out, rep(NA, 4))
                P_donor <- NA
            } else {
                gm1 <- glmer(cbind(altCount, refCount) ~ diff_population +
                    (1 | experiment),
                    data = dat4glmm, family = binomial(link = "logit"))
                out <- c(out, AIC(gm1), BIC(gm1), logLik(gm1), deviance(gm1))
                P_diffpop <- anova(gm0, gm1)$P[2]
            }
            ### test ouija_pseudotime
            gm1 <- glmer(cbind(altCount, refCount) ~ ouija_pseudotime +
                    (1 | experiment),
                    data = dat4glmm, family = binomial(link = "logit"))
                    out <- c(out, AIC(gm1), BIC(gm1), logLik(gm1), deviance(gm1))
            P_pt <- anova(gm0, gm1)$P[2]
            ### add P-values to output vector
            out <- c(P_donor, P_diffpop, P_pt, out)
        }
    return(out)
    },
    error = function(cond) {
        return(rep(NA, 20))
    },
    finally = NULL)
    out
}



## exprs data
input_exprs <- "data_processed/merged/sceset_merged_qc.rds"
sce <- readRDS(input_exprs)
## order cells by ouija pseudotime
oo <- order(sce$ouija_pseudotime)
sce <- sce[, oo]

## load total counts to filter sites early
input_totalCount <-
    "data_processed/ase/ase_counts.low_thresh.aggregated_totalCount.rds"
totalCount <- readRDS(input_totalCount)
## drop variants that are not observed with sufficient depth in >=30 cells
## 850312 sites called in at least one cell
## 14192 sites called in at least 30 cells
keep_site <- (rowSums(totalCount > 0) >= 30)
totalCount <- totalCount[keep_site, ]
## keep cells that passed QC in the expression data
colnames_totalCount <- gsub("low_thresh/", "",
                    gsub(".ase.lowthresh.tsv", "", colnames(totalCount)))
keep_cell <- (colnames_totalCount %in% colnames(sce))
rm(totalCount)


## load SummarizedExperiment data
input_ase <- "data_processed/ase/ase_counts.low_thresh.aggregated.rds"
ase <- readRDS(input_ase)
ase <- ase[keep_site, ]
colnames(ase) <- gsub("low_thresh/", "",
                    gsub(".ase.lowthresh.tsv", "", colnames(ase)))
ase <- ase[, colnames_totalCount[keep_cell]]
## read in altCount data and add to ase object
input_alt_count <-
    "data_processed/ase/ase_counts.low_thresh.aggregated_altCount.rds"
mat <- readRDS(input_alt_count)
mat <- mat[keep_site, ]
colnames(mat) <- gsub("low_thresh/", "",
                    gsub(".ase.lowthresh.tsv", "", colnames(mat)))
mat <- mat[, colnames(ase)]
assays(ase)$altCount <- mat
rm(mat)
## read in altCount data and add to ase object
input_ref_count <-
    "data_processed/ase/ase_counts.low_thresh.aggregated_refCount.rds"
mat <- readRDS(input_ref_count)
mat <- mat[keep_site, ]
colnames(mat) <- gsub("low_thresh/", "",
                    gsub(".ase.lowthresh.tsv", "", colnames(mat)))
mat <- mat[, colnames(ase)]
assays(ase)$refCount <- mat
rm(mat)


## get the number of cells with sufficient depth to be of interest
mm <- match(colnames(ase), colnames(sce))
identical(rownames(pData(sce)[mm, ]), colnames(ase))
colData(ase) <- DataFrame(pData(sce)[mm, ])

## save intermediate file
saveRDS(ase, file =
    "data_processed/ase/ase_counts.low_thresh.aggregated.qc.rds")



## load summarised data
ase <- readRDS("data_processed/ase/ase_counts.low_thresh.aggregated.qc.rds")
pdata <- readRDS("data_processed/merged/pdata_merged_qc_for_ase.rds")
identical(colnames(ase), rownames(pdata))
ase$ouija_pseudotime <- pdata$ouija_pseudotime
head(ase$ouija_pseudotime)

# dat_list <- vector("list", nrow(ase))
# names(dat_list) <- rownames(ase)
# pb <- txtProgressBar(min = 0, max = length(dat_list), style = 3)
# counter <- 0
# for (i in seq_len(nrow(ase))) {
#     counter <- counter + 1
#     dat4glmm <- cbind(refCount = assays(ase)$refCount[i, ],
#                     altCount = assays(ase)$altCount[i, ])
#     dat_list[[i]] <- dat4glmm
#     setTxtProgressBar(pb, counter)
# }
# close(pb)
# saveRDS(dat_list, file = "data_processed/ase/data_list_for_ase_glmm.qc.rds")
dat_list <- readRDS("data_processed/ase/data_list_for_ase_glmm.qc.rds")

## try with Michael Schubert's clustermq package
Q(FUN, i=1:10, n_jobs=1)


## with bplapply
fit_all <- bplapply(dat_list, fit_ase_glmm, BPPARAM = multicoreParam, coldata = colData(ase))

## with clustermq
fit_all <- Q(fit_ase_glmm_cmq, i = seq_len(nrow(ase)),
                memory = 16000, n_jobs = 1000, job_size = 20)

fit_all_df <- as_data_frame(t(as_tibble(fit_all, validate = FALSE)))
colnames(fit_all_df) <- c("P_donor", "P_diffpop", "P_pt", "GLMM_fit",
                paste(rep(c("AIC", "BIC", "logLik", "deviance"), 4),
                rep(c("null", "donor", "diffpop", "pt")), sep = "_"))
fit_all_df <- bind_cols(fit_all_df, as.data.frame(mcols(rowRanges(ase))))
fit_all_df[["P_donor"]][is.na(fit_all_df[["P_donor"]])] <- 1
fit_all_df[["P_diffpop"]][is.na(fit_all_df[["P_diffpop"]])] <- 1
fit_all_df[["P_pt"]][is.na(fit_all_df[["P_pt"]])] <- 1

fit_all_df


saveRDS(fit_all_df, file = "data_processed/ase/ase_bin_glmm_results.rds")
saveRDS(ase_filt,
            file = "data_processed/ase/ase_counts.aggregated.filtered.rds")


############


fit_ase_glmm_cmq <- function(i) {
    library(SummarizedExperiment)
    library(scater)
    library(lme4)
    library(tibble)
    library(dplyr)
    ase <- readRDS("data_processed/ase/ase_counts.low_thresh.aggregated.qc.rds")
    pdata <- readRDS("data_processed/merged/pdata_merged_qc_for_ase.rds")
    identical(colnames(ase), rownames(pdata))
    ase$ouija_pseudotime <- pdata$ouija_pseudotime
    dat4glmm <- cbind(refCount = assays(ase)$refCount[i,],
                        altCount = assays(ase)$altCount[i,],
                        as.data.frame(colData(ase)))
    dat4glmm <- dat4glmm[!is.na(dat4glmm$refCount) &
                            !is.na(dat4glmm$altCount), ]
    all_zero <- (rowSums(dat4glmm[, 1:2]) == 0)
    dat4glmm <- dat4glmm[!all_zero,]
    out <- tryCatch({
        ## define data object for GLMM fit
        if (length(unique(dat4glmm$experiment)) < 2L) {
            ## only one environment, so don't fit random effect
            gm0 <- glm(cbind(altCount, refCount) ~ 1,
                        data = dat4glmm, family = binomial(link = "logit"))
            out <- c(0, AIC(gm0), BIC(gm0), logLik(gm0), deviance(gm0))
            ### test donor
            if (length(unique(dat4glmm$donor)) < 2L) {
                out <- c(out, rep(NA, 4))
                P_donor <- NA
            } else {
                gm1 <- glm(cbind(altCount, refCount) ~ donor,
                                data = dat4glmm, family = binomial(link = "logit"))
                out <- c(out, AIC(gm1), BIC(gm1), logLik(gm1), deviance(gm1))
                P_donor <- 1 - pchisq(gm0$deviance - gm1$deviance,
                                gm0$df.residual - gm1$df.residual)
            }
            ### test differentiation population
            if (length(unique(dat4glmm$diff_population)) < 2L) {
                out <- c(out, rep(NA, 4))
                P_donor <- NA
            } else {
                gm1 <- glm(cbind(altCount, refCount) ~ diff_population,
                                data = dat4glmm, family = binomial(link = "logit"))
                out <- c(out, AIC(gm1), BIC(gm1), logLik(gm1), deviance(gm1))
                P_diffpop <- 1 - pchisq(gm0$deviance - gm1$deviance,
                                gm0$df.residual - gm1$df.residual)
            }
            ### test ouija_pseudotime
            gm1 <- glm(cbind(altCount, refCount) ~ ouija_pseudotime,
                        data = dat4glmm, family = binomial(link = "logit"))
            out <- c(out, AIC(gm1), BIC(gm1), logLik(gm1), deviance(gm1))
            P_pt <- 1 - pchisq(gm0$deviance - gm1$deviance,
                            gm0$df.residual - gm1$df.residual)
            ### add P-values to output vector
            out <- c(P_donor, P_diffpop, P_pt, out)
        } else {
            ## fit GLMM with glmer
            gm0 <- glmer(cbind(altCount, refCount) ~ (1 | experiment),
                            data = dat4glmm, family = binomial(link = "logit"))
            out <- c(1, AIC(gm0), BIC(gm0), logLik(gm0), deviance(gm0))
            ### test donor
            if (length(unique(dat4glmm$donor)) < 2L) {
                out <- c(out, rep(NA, 4))
                P_donor <- NA
            } else {
                gm1 <- glmer(cbind(altCount, refCount) ~ donor + (1 | experiment),
                        data = dat4glmm, family = binomial(link = "logit"))
                out <- c(out, AIC(gm1), BIC(gm1), logLik(gm1), deviance(gm1))
                P_donor <- anova(gm0, gm1)$P[2]
            }
            ### test differentiation population
            if (length(unique(dat4glmm$diff_population)) < 2L) {
                out <- c(out, rep(NA, 4))
                P_donor <- NA
            } else {
                gm1 <- glmer(cbind(altCount, refCount) ~ diff_population +
                    (1 | experiment),
                    data = dat4glmm, family = binomial(link = "logit"))
                out <- c(out, AIC(gm1), BIC(gm1), logLik(gm1), deviance(gm1))
                P_diffpop <- anova(gm0, gm1)$P[2]
            }
            ### test ouija_pseudotime
            gm1 <- glmer(cbind(altCount, refCount) ~ ouija_pseudotime +
                    (1 | experiment),
                    data = dat4glmm, family = binomial(link = "logit"))
                    out <- c(out, AIC(gm1), BIC(gm1), logLik(gm1), deviance(gm1))
            P_pt <- anova(gm0, gm1)$P[2]
            ### add P-values to output vector
            out <- c(P_donor, P_diffpop, P_pt, out)
        }
    return(out)
    },
    error = function(cond) {
        return(rep(NA, 20))
    },
    finally = NULL)
    out
}
