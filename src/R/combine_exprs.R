## function to combine expression across two lanes of sequencing data

combine_exprs <- function(sceset_object, sample = sceset_object$sample) {
    ## Initialise new TPM and count matrices
    tpm_new <- counts_new <- matrix(
        NA, nrow = nrow(sceset_object),
        ncol = length(unique(sample)))
    colnames(tpm_new) <- colnames(counts_new) <- unique(sample)
    rownames(tpm_new) <- rownames(counts_new) <- rownames(sceset_object)
    if ( !is.null(sceset_object$n_obs_fragments) ) {
        nfrags <- rep(NA, length(unique(sample)))
        names(nfrags) <- unique(sample)
    }
    if ( !is.null(sceset_object$mapping_rate_pct) ) {
        maprate <- rep(NA, length(unique(sample)))
        names(maprate) <- unique(sample)
    }
    if ( !is.null(sceset_object$n_processed) ) {
        nprocessed <- rep(NA, length(unique(sample)))
        names(nprocessed) <- unique(sample)
    }
    ## average TPM and sum counts across replicates
    for (samp in unique(sample)) {
        tpm_tmp <- rowMeans(tpm(
            sceset_object)[, sample == samp])
        tpm_new[, samp] <- tpm_tmp
        count_tmp <- rowSums(counts(
            sceset_object)[, sample == samp])
        counts_new[, samp] <- count_tmp
        if ( !is.null(sceset_object$mapping_rate_pct) )
            maprate[samp] <- mean(
                sceset_object$mapping_rate_pct[sample == samp])
        if ( !is.null(sceset_object$n_obs_fragments) )
            nfrags[samp] <- sum(
                sceset_object$n_obs_fragments[sample == samp])
        if ( !is.null(sceset_object$n_processed) )
            nprocessed[samp] <- sum(
                sceset_object$n_processed[sample == samp])
        
    }
    ## Define new SCESet object and return
    pd <- pData(sceset_object)[!duplicated(sample),]
    rownames(pd) <- pd$sample
    pd$lane <- "combined"
    if ( !is.null(sceset_object$mapping_rate_pct) )
        pd$mapping_rate_pct <- maprate
    if ( !is.null(sceset_object$n_obs_fragments) )
        pd$n_obs_fragments <- nfrags
    if ( !is.null(sceset_object$n_processed) )
        pd$n_processed <- nprocessed
    pd <- new("AnnotatedDataFrame", pd)
    if ( identical(rownames(pd), colnames(tpm_new)) )
        print("Sample names match")
    else
        stop("Sample names don't match")
    fd <- new("AnnotatedDataFrame", fData(sceset_object))
    newSCESet(tpmData = tpm_new, countData = counts_new,
              phenoData = pd,
              featureData = fd)
}

