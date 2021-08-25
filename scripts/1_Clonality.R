# Analysis scripts for Cagan, Baez-Ortega et al., 2021
# Step 1: Identification of likely polyclonal samples

# Adrian Baez-Ortega, 2020-21


# TO RUN THIS SCRIPT IN THE TERMINAL
# ----------------------------------
# Run the commands below in the terminal, replacing '/path/to/CrossSpecies2021' with the 
# path to the CrossSpecies2021 directory.
#
#    cd /path/to/CrossSpecies2021
#    Rscript scripts/1_Clonality.R


# TO RUN THIS SCRIPT IN RSTUDIO
# -----------------------------
# Before starting, run the line below in RStudio, replacing '/path/to/CrossSpecies2021' with 
# the path to the CrossSpecies2021 directory.
#
#    setwd("path/to/CrossSpecies2021")


# If the paths to the input or output files differ from the ones
# used in this script, they may be updated by modifying the lines below.


INPUT = list(
    
    # Path to project info table
    SAMPLE.INFO = file.path("data", "original", "CrossSpecies_ProjectInfo.txt"),
    
    # Path to table of samples per individual
    SAMPLE.MATCH = file.path("data", "original", "CrossSpecies_SamplesPerIndividual.txt"),
    
    # Path to list of samples to 'white-list'
    WHITELIST = file.path("data", "original", "CrossSpecies_SampleWhitelist.txt"),
    
    # Generic path to allele count files
    # (obtained with the alleleCounter software; see Suppl. Methods)    
    COUNTS = file.path("data", "original", "AlleleCounts_Subs", "AlleleCounts_Subs_${SAMPLE}.txt"),
    
    # Generic path to files of filtered substitutions
    # (see Suppl. Methods for details on variant calling and filtering)    
    SUBS.FINAL = file.path("data", "original", "FinalCalls_Subs", "${SPECIES}",
                           "FinalCalls_Subs_${SAMPLE}.txt"),
    
    # Path to file of filtered mtDNA variants
    MTDNA.FINAL = file.path("data", "original", "FinalCalls_mtDNA", "FinalCalls_mtDNA_Merged.txt")

)

OUTPUT = list(
    
    # Path to PDF of VAF histograms
    HIST = file.path("output", "Clonality_VAF_Histograms.pdf"),
    
    # Path to table of filtering results
    TABLE = file.path("output", "Clonality_Results_Table.txt"),
    
    # Path to list of samples to exclude from all analyses
    EXCLUDE = file.path("data", "processed", "SamplesToExclude.txt"),
    
    # Path to list of samples to exclude from mtDNA analyses
    EXCLUDE.MT = file.path("data", "processed", "SamplesToExclude_mtDNA.txt")

)


# Threshold parameters:
# Minimum number of supporting reads
MIN.NV = 4

# Minimum number of variants per sample
# (must be greater than the maximum number of clusters)
MIN.VARS = 6

# Burden and cluster VAF thresholds for sample filtering
HIGH.BURDEN.THR = 3.0
MIN.BURDEN = 50
MIN.PROP = 0.7
MIN.VAF = 0.3

# Thresholds for mtDNA sample filtering: Samples are excluded if the number of
# variants with VAF < MT.VAF is greater than MT.FACTOR times the overall mean
MT.VAF = 0.02
MT.FACTOR = 5


# Print input and output file paths
cat("\nInput files:")
for (name in INPUT) {
    cat("\n  ", name)
}
cat("\n\nOutput files:")
for (name in OUTPUT) {
    cat("\n  ", name)
}
cat("\n\n")


cat("Loading data...\n")
whitelist = as.character(as.matrix(read.table(INPUT$WHITELIST)))
mt.vars = read.table(INPUT$MTDNA.FINAL, sep="\t", header=T, as.is=T)
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
sample.match = read.table(INPUT$SAMPLE.MATCH, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
stopifnot(!any(duplicated(sample.match$SAMPLE_NAME)))
stopifnot(all(sample.info$SAMPLE_NAME %in% sample.match$SAMPLE_NAME))
rownames(sample.match) = sample.match$SAMPLE_NAME
cat("Loaded\n")


# In order to identify samples for which the read support distribution is
# suggestive of polyclonality, a truncated binomial mixture model is fitted
# using an expectation-maximisation algorith. The EM code below is adapted
# from original code by T.H.H. COORENS (Wellcome Sanger Institute).


# Truncated binomial density function
# (by default, the minimum number of supporting reads is minx=4)
dbinomtrunc = function(x, size, prob, minx=4) {
    dbinom(x, size, prob) / pbinom(minx-0.1, size, prob, lower.tail=F)
}


# Expectation step of EM algorithm
estep = function(x, size, p.vector, prop.vector, ncomp, mode) {
    # p.vector: vector of probabilities for the individual components
    # prop.vector: vector of proportions for the individual components
    # ncomp: number of components
    # mode: form of binomial distribution ("full" or "truncated")
    p.mat.estep = matrix(0, ncol=ncomp, nrow=length(x))
    for (i in 1:ncomp){
        if (mode != "full") {
            p.mat.estep[, i] = prop.vector[i] * dbinomtrunc(x, size, prob=p.vector[i])
        }
        else {
            p.mat.estep[, i] = prop.vector[i] * dbinom(x, size, prob=p.vector[i])
        }
    }
    
    # Normalise probabilities and obtain log-likelihood
    norm = rowSums(p.mat.estep)
    p.mat.estep = p.mat.estep / norm
    LL = sum(log(norm))
    
    # Assign observations to components
    which.clust = rep(1, length(x))
    if (ncomp > 1) {
        which.clust = apply(p.mat.estep, 1, which.max)
    }
    
    list("posterior"=p.mat.estep, "LL"=LL, "which.cluster"=which.clust)
}


# Maximisation step of EM algorithm
mstep = function(x, size, e.step){
    # Estimate proportions and probabilities
    list("prop" = colMeans(e.step$posterior),
         "p" = colSums(x / size * e.step$posterior) / colSums(e.step$posterior))
}


# EM algorithm
em.algorithm = function(x, size, prop.vector.inits, p.vector.inits,
                        maxit = 1e4, tol = 1e-6, nclust, binom.mode) {
    # prop.vector.inits: initial values for the mixture proportions
    # p.vector.inits: initial values for the probabilities
    
    # Initialise EM
    flag = FALSE
    e.step = estep(x, size, p.vector=p.vector.inits, prop.vector=prop.vector.inits,
                   ncomp=nclust, mode=binom.mode)
    m.step = mstep(x, size, e.step)
    prop.cur = m.step[["prop"]]
    p.cur = m.step[["p"]]
    LL.cur = e.step[["LL"]]
    LL.vector = e.step[["LL"]]
    
    # Iterate between expectation and maximisation steps
    for (i in 2:maxit){
        e.step = estep(x, size, p.vector=p.cur, prop.vector=prop.cur,
                       ncomp=nclust, mode=binom.mode)
        m.step = mstep(x, size, e.step)
        prop.new = m.step[["prop"]]
        p.new = m.step[["p"]]
        LL.new = e.step[["LL"]]
        LL.vector = c(LL.vector, LL.new)
        LL.diff = abs((LL.cur - LL.new))
        which.clust = e.step[["which.cluster"]]
        
        # Stop iteration if the difference between the current
        # and new log-likelihoods is below a tolerance level
        if (LL.diff < tol) {
            flag = TRUE
            break
        }
        # Otherwise, continue iteration
        prop.cur = prop.new
        p.cur = p.new
        LL.cur = LL.new
    }
    
    if (!flag) warning("EM algorithm did not converge\n")
    
    BIC = 2 * log(length(x)) * nclust - 2 * LL.cur
    AIC = 4 * nclust - 2 * LL.cur
    list("LL"=LL.vector,
         "prop"=prop.cur,
         "p"=p.cur,
         "BIC"=BIC,
         "AIC"=AIC,
         "n"=nclust,
         "which.cluster"=which.clust)
}


# Binomial mixture model core function
# This function runs the EM algorithm for different numbers of components, and selects the best
# fit using the Bayesian information criterion (BIC) or the Akaike information criterion (AIC).
# Returns a list with the following information:
#   LL: log-likelihood vector of the iteration
#   prop: the optimal mixing proportion of clones
#   p: the probabilities (VAFs) associated with each clone
#   BIC: the Bayesian information criterion associated with best fit
#   AIC: the Akaike information criterion associated with best fit
#   n: the optimal number of clusters/clones
#   which.cluster: same length as input; cluster to which each mutation belongs
#   BIC.vec: the vector of BIC values for entire range
#   AIC.vec: the vector of AIC values for entire range
binom.mix = function(x, size, nrange = 1:5, criterion = "BIC",
                     maxit = 1e4, tol = 1e-6, mode = "truncated") {
    stopifnot(criterion %in% c("BIC", "AIC"))
    results = vector(mode="list", length=length(nrange))
    BIC.vec = rep(NA, length(nrange))
    AIC.vec = rep(NA, length(nrange))
    
    for (i in 1:length(nrange)) {
        # Use values from k-means clustering as EM inits
        n = nrange[i]
        init = kmeans(x/size, n)
        prop.init = init$size / length(x)
        p.init = init$centers
        # Run EM algorithm
        results[[i]] = em.algorithm(x, size, prop.init, p.init, maxit, tol, n, mode)
        BIC.vec[i] = results[[i]]$BIC
        AIC.vec[i] = results[[i]]$AIC
    }
    
    # Return best fit
    idx = switch(criterion,
                 "BIC"=which.min(BIC.vec),
                 "AIC"=which.min(AIC.vec))
    cat("Best fit: n =", results[[idx]]$n, "\n")
    results[[idx]]$BIC.vec = BIC.vec
    results[[idx]]$AIC.vec = AIC.vec
    results[[idx]]
}


# Initialise VAF histograms PDF
pdf(OUTPUT$HIST, 12, 6)
par(mar=c(4.5, 4.5, 4.5, 0.5), mgp=c(2.2, 0.75, 0))


# Run binomial mixture model for each sample
sample.filt.table = data.frame("Sample"=NA, "Species"=NA, "Mutations"=NA, "Median indiv. burden"=NA,
                               "Clusters"=NA, "Main cluster VAF"=NA, "Prop. with VAF ≥ 0.3"=NA,
                               check.names=F, stringsAsFactors=F)
binom.mix.out = structure(vector(mode="list", length=nrow(sample.info)),
                          names=sample.info$SAMPLE_NAME)

for (i in 1:nrow(sample.info)) {
    sample.id = sample.info$SAMPLE_NAME[i]
    normal.id = sample.info$NORMAL_NAME[i]
    species = sample.info$SPECIES_NAME[i]
    cat("\nProcessing sample", sample.id, "\n")
    
    # Check file paths
    paths = c(gsub("${SAMPLE}", sample.id,
                   gsub("${SPECIES}", species,
                        INPUT$SUBS.FINAL, fixed=T), fixed=T),
              gsub("${SAMPLE}", sample.id, INPUT$COUNTS, fixed=T))
    
    if (!all(file.exists(paths))) {
        cat("WARNING: File(s)", paste(paths[!file.exists(paths)], collapse=", "),
            "not found. Skipping sample.\n")
        next
    }
    
    # Load filtered variants and read counts
    variants = read.table(paths[1], sep="\t", header=T, check.names=F, as.is=T, comment.char="",
                          colClasses=c("Chr"="character", "Ref"="character", "Alt"="character"))
    counts = read.table(paths[2], header=T, comment.char="", check.names=F, as.is=T)
    cat(nrow(variants), "variants read from input file\n")
    
    if (nrow(variants) < MIN.VARS) {
        cat("WARNING: Insufficient variants in input file. Skipping sample.\n")
        sample.filt.table = rbind(sample.filt.table,
                                  data.frame("Sample"=sample.id, "Species"=species,
                                             "Mutations"=nrow(variants),
                                             "Median indiv. burden"=NA, "Clusters"=NA,
                                             "Main cluster VAF"=NA, "Prop. with VAF ≥ 0.3"=NA,
                                             check.names=F, stringsAsFactors=F))
        next
    }
    
    # Calculate median mut. burden across matched samples (from same individual)
    idx = match(c("NORMAL_NAME", "SPECIES_NAME"), colnames(sample.match))
    matched.samples = as.character(sample.match[sample.id, -idx])
    matched.samples = matched.samples[matched.samples != ""]
    median.burden = median(sapply(matched.samples, function(smp) {
        fpath = gsub("${SAMPLE}", smp,
                     gsub("${SPECIES}", species,
                          INPUT$SUBS.FINAL, fixed=T), fixed=T)
        if (file.exists(fpath)) {
            vars = read.table(fpath, sep="\t", header=T, comment.char="")
            ifelse(nrow(vars) >= MIN.BURDEN, nrow(vars), NA)
        }
        else {
            NA
        }
    }), na.rm=T)
    
    # Obtain NR (total reads) and NV (supporting reads) for filtered variants
    idx = match(paste0(variants$Chr, ":", variants$Start), paste0(counts$`#CHR`, ":", counts$POS))
    stopifnot(!any(is.na(idx)) & !any(duplicated(counts[, 1:2])))
    vars.nr = counts$Good_depth[idx]
    vars.nv = sapply(1:length(idx), function(j) {
        counts[idx[j], paste0("Count_", variants$Alt[j])]
    })
    vars.vaf = vars.nv / vars.nr
    stopifnot(length(vars.vaf) == nrow(variants) & !any(is.na(vars.vaf)) &
                  !any(vars.nr == 0) & !any(vars.nv > vars.nr))
    
    # Run truncated binomial mixture model
    if (sum(vars.nv >= MIN.NV) < MIN.VARS) {
        cat("WARNING: Insufficient variants with NV>=", MIN.NV, ". Skipping sample.\n", sep="")
        sample.filt.table = rbind(sample.filt.table,
                                  data.frame("Sample"=sample.id, "Species"=species,
                                             "Mutations"=nrow(variants),
                                             "Median indiv. burden"=NA, "Clusters"=NA,
                                             "Main cluster VAF"=NA, "Prop. with VAF ≥ 0.3"=NA,
                                             check.names=F, stringsAsFactors=F))
        next
    }
    vars.vaf = vars.vaf[vars.nv >= MIN.NV]
    vars.nr = vars.nr[vars.nv >= MIN.NV]
    vars.nv = vars.nv[vars.nv >= MIN.NV]
    out = binom.mix(vars.nv, vars.nr)
    
    # Plot VAF and fitted model
    set.seed(0xC0FFEE)
    h = hist(vars.vaf, breaks=seq(0, 1, 0.01), freq=F, col="dodgerblue4", border="white", xlab="VAF",
             main=paste0("\nBest fit for ", sample.id, "\nn = ", out$n,
                         ",   proportions = {", paste(round(out$prop, 2), collapse=", "),
                         "},   VAFs = {", paste(round(out$p, 2), collapse=", "),
                         "},   BIC = ", round(out$BIC, 2), "\n",
                         prettyNum(length(vars.vaf), big.mark=","), " variants with NV >= ", MIN.NV))
    for (j in 1:out$n) {
        nr = rpois(5000, lambda=median(vars.nr))
        nv = unlist(lapply(nr, rbinom, n=1, prob=out$p[j]))
        dens = density((nv / nr)[nv >= MIN.NV])
        lines(dens$x, out$prop[j] * dens$y, col="red2", lwd=3, lty=2)
    }
    
    # Add sample to output table
    sample.filt.table = rbind(sample.filt.table,
                              data.frame("Sample"=sample.id, "Species"=species,
                                         "Mutations"=nrow(variants),
                                         "Median indiv. burden"=median.burden, "Clusters"=out$n,
                                         "Main cluster VAF"=round(out$p[which.max(out$prop)], 3),
                                         "Prop. with VAF ≥ 0.3"=round(sum(out$prop[out$p >= MIN.VAF]), 3),
                                         check.names=F, stringsAsFactors=F))
    # Store model output
    binom.mix.out[[i]] = out
}

sample.filt.table = sample.filt.table[-1, ]
binom.mix.out = binom.mix.out[!sapply(binom.mix.out, is.null)]
invisible(dev.off())


# Build sample filter indices
# We use the following filters:
# 1) Low burden: discard samples with <50 mutations
# 2) High burden: discard samples with burden >3 * median burden across the same individual
#                 (samples with <50 mutations are not considered when calculating the median)
# 3) Low VAF: discard samples where the total proportion of variants in clusters with VAF ≥0.3
#             is <0.7, unless all the samples in the same individual have VAF <0.3
low.burden.idx = sample.filt.table$Mutations < MIN.BURDEN
high.burden.idx = (sample.filt.table$Mutations >
    sample.filt.table$`Median indiv. burden` * HIGH.BURDEN.THR) %in% TRUE

low.vaf.idx = sapply(1:nrow(sample.filt.table), function(i) {
    idx = match(c("NORMAL_NAME", "SPECIES_NAME"), colnames(sample.match))
    matched.samples = as.character(sample.match[sample.filt.table$Sample[i], -idx])
    matched.samples = matched.samples[matched.samples != ""]
    matched.vafs = sample.filt.table$`Main cluster VAF`[match(matched.samples,
                                                              sample.filt.table$Sample)]
    (sample.filt.table$`Prop. with VAF ≥ 0.3`[i] < MIN.PROP) %in% TRUE &
        any(matched.vafs >= MIN.VAF, na.rm=T)
})


# Output sample filtering table
sample.filt.table = cbind(sample.filt.table,
                          "<50 mutations"=low.burden.idx,
                          "Burden > 3x median"=high.burden.idx,
                          "<70% muts. with VAF ≥ 0.3"=low.vaf.idx)
write.table(sample.filt.table, file=OUTPUT$TABLE, sep="\t", col.names=T, row.names=F, quote=F)

# Output list of samples to exclude
exclude.samples = sample.filt.table$Sample[(low.burden.idx | high.burden.idx | low.vaf.idx) &
                                               !(sample.filt.table$Sample %in% whitelist)]
cat("\nSamples excluded from all analyses:", exclude.samples, sep="\n")
cat(exclude.samples, file=OUTPUT$EXCLUDE, sep="\n")


# For mtDNA analyses only: exclude samples whose number of variants
# with VAF < 0.02 is greater than 5 times the overall mean

# Count number of mtDNA variants with VAF < 0.02 (inc. 79 samples with zero variants)
mtdna.counts = c(table(mt.vars$Sample[mt.vars$VAF < MT.VAF]), rep(0, 79))

# Output list of samples with more than 5*mean(mtdna.counts) = 4.36 variants
exclude.samples.mt = names(mtdna.counts[mtdna.counts > MT.FACTOR * mean(mtdna.counts)])

cat("\nSamples excluded from mtDNA analyses only:", exclude.samples.mt, sep="\n")
cat(exclude.samples.mt, file=OUTPUT$EXCLUDE.MT, sep="\n")


cat("\nDone\n\n")
