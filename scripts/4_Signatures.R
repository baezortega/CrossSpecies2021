# Analysis scripts for Cagan, Baez-Ortega et al., 2022
# Step 4: Inference of mutational signatures from somatic substitutions

# Adrian Baez-Ortega, 2020-22


# TO RUN THIS SCRIPT IN THE TERMINAL
# ----------------------------------
# Run the commands below in the terminal, replacing '/path/to/CrossSpecies2021' with the 
# path to the CrossSpecies2021 directory.
#
#    cd /path/to/CrossSpecies2021
#    Rscript scripts/4_Signatures.R


# TO RUN THIS SCRIPT IN RSTUDIO
# -----------------------------
# Before starting, run the line below in RStudio, replacing '/path/to/CrossSpecies2021' with 
# the path to the CrossSpecies2021 directory.
#
#    setwd("path/to/CrossSpecies2021")


# If the paths to the input or output files differ from the ones
# used in this script, they may be updated by modifying the lines below.


# Input file paths
INPUT = list(
    
    # Path to project info table
    SAMPLE.INFO = file.path("data", "original", "CrossSpecies_ProjectInfo.txt"),
    
    # Path to list of samples to exclude
    EXCLUDE = file.path("data", "processed", "SamplesToExclude.txt"),
    
    # Path to RData file of mutational catalogues and opportunities
    COUNTS = file.path("data", "processed", "Catalogues_Opportunities.RData")
    
)

# Output file paths
OUTPUT = list(
    
    # Paths to preliminary signature extraction results (not run by default)
    # PDF.DIR.PRELIM = file.path("output", "Signature_Extraction_Preliminary"),
    # PDF.NORM.PRELIM = "Species_Signatures_HumanNormalised.pdf",
    # GOF.PRELIM = "Species_GOF.pdf",
    
    # Paths to final signature extraction results
    PDF.DIR = file.path("output", "Signature_Extraction_Definitive"),
    PDF.EXP = "Exposures_Sample_Definitive.pdf",
    PDF.REC = "Reconstructions_Sample_Definitive.pdf",
    PDF.SIGS = "Signatures_Definitive.pdf",
    PDF.SIGS.NORM = "Signatures_Definitive_Human-Normalised.pdf",
    
    # Path to PDF of cross-species SBSB spectra
    PDF.SBSB = file.path("output", "SBSB_Per_Species.pdf"),
    
    # Path to PDF of colibactin exposures
    PDF.SBS88 = file.path("output", "Colibactin_Exposure.pdf"),
    
    # Path to colibactin/APOBEC test results
    FISHER.TEST = file.path("output", "Colibactin_APOBEC_Tests.txt"),
    
    # Path to RData file of mutation burden and exposures per sample
    BURDEN.DATA = file.path("data", "processed", "Burden_Exposures.RData"),
    
    # Path to RData file of mutational signatures and exposures
    SIGS.DATA = file.path("data", "processed", "Signatures_Definitive.RData")
    
)


# sigfit parameters
ITER.1 = 15000                   # Total iterations (extraction)
ITER.2 = 10000                   # Total iterations (re-fitting)
WARMUP = 5000                    # Warmup iterations
CHAINS = 3                       # MCMC chains (re-fitting)
SEED = 0xC0FFEE                  # Random seed
CORES = parallel::detectCores()  # Number of CPUs (detected automatically)

# Final signature names
SIG.NAMES = c("SBS1", "SBSB", "SBSC")

# Index for reordering final signatures after extraction
SIG.ORDER = c(1, 3, 2)

# Function: reverse complement (for string vectors)
rev.comp = function(nucleotide.list) {
    sapply(nucleotide.list, function(nucleotides) {
        paste(
            rev(sapply(strsplit(nucleotides, "")[[1]], function(nuc) {
                if (nuc == "A") "T"
                else if (nuc == "C") "G"
                else if (nuc == "G") "C"
                else if (nuc == "T") "A"
            })),
            collapse="")
    }, USE.NAMES=F)
}

# Function: normalise vector to sum to one
normalise = function(x) {
    x / sum(x)
}


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


# Load packages and data
PACKAGES = c("scales", "sigfit", "tools")
cat("Loading packages:", paste(PACKAGES, collapse=", "), "\n")
for (package in PACKAGES) {
    suppressPackageStartupMessages(library(package, character.only=TRUE, quietly=TRUE))
}

cat("Loading data...\n")
load(INPUT$COUNTS)
data("cosmic_signatures_v3", package="sigfit")
exclude.list = as.character(as.matrix(read.table(INPUT$EXCLUDE)))
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
rownames(sample.counts) = rownames(sample.opps)
rownames(species.counts) = rownames(species.opps)
cat("Loaded\n")


# Create output directory
dir.create(OUTPUT$PDF.DIR, showWarnings=F)


# (0) Preliminary signature analysis (optional; uncomment code below to run)
# To determine the best number of signatures, signatures were initially extracted
# de novo from catalogues of mutations per species, per individual and per sample.
# Below is the code to extract 2-7 signatures from the catalogues per species,
# and plot the results (including the estimated best number of signatures).

# cat("\n\nRunning preliminary signature extraction (2-7 signatures):\n\n")
# dir.create(OUTPUT$PDF.DIR.PRELIM, showWarnings=F)
# 
# # Extract 2-7 signatures from species catalogues (using mutational opportunities)
# pdf(file.path(OUTPUT$PDF.DIR.PRELIM, OUTPUT$GOF.PRELIM), 8, 5)
# fit.species = suppressWarnings(extract_signatures(species.counts,
#                                                   nsignatures=2:7,
#                                                   opportunities=species.opps,
#                                                   iter=8000, warmup=3000, seed=SEED))
# invisible(dev.off())
# 
# cat("\nPlotting preliminary extraction results:\n")
# for (i in NSIGS) {
#     plot_all(fit.species[[i]], out_path=file.path(OUTPUT$PDF.DIR.PRELIM, i), prefix="Species")
#     plot_spectrum(convert_signatures(retrieve_pars(fit.species[[i]], "signatures"),
#                                      opportunities_to="human-genome"),
#                   pdf_path=file.path(OUTPUT$PDF.DIR.PRELIM, i, OUTPUT$PDF.NORM.PRELIM))
# }



# (1) Final signature analysis
# Optimal extraction results are achieved by extracting 3 signatures from
# species catalogues (see preliminary analysis above).
# To avoid the mixing of SBS1 and the other signatures (especially SBSB), the
# definitive set of signatures is obtained by fitting COSMIC SBS1 and extracting
# two additional signatures (SBSB, SBSC) using the Fit-Ext model in sigfit.

cat("\n\nRunning definitive signature extraction (3 signatures):\n\n")
fitext.species = suppressWarnings(
    fit_extract_signatures(species.counts,
                           convert_signatures(cosmic_signatures_v3[1, ],
                                              opportunities_from="human-genome"),
                           num_extra_sigs=2,
                           opportunities=species.opps,
                           iter=ITER.1, warmup=WARMUP, seed=SEED)
)

# Retrieve and reorder final signatures (SBS1, SBSB, SBSC)
signatures.final = retrieve_pars(fitext.species, "signatures")
for (i in 1:length(signatures.final)) {
    signatures.final[[i]] = signatures.final[[i]][SIG.ORDER, ]
    rownames(signatures.final[[i]]) = SIG.NAMES
}

# Fit inferred signatures to sample catalogues to obtain signature exposures
cat("\n\nInferring signature exposures per sample:\n\n")
fit.sample = suppressWarnings(fit_signatures(sample.counts,
                                             signatures.final,
                                             opportunities=sample.opps,
                                             iter=ITER.2, warmup=WARMUP, seed=SEED,
                                             chains=CHAINS, cores=CORES))

# Retrieve exposures and reconstructions
exposures.final = retrieve_pars(fit.sample, "exposures")
reconstructions.final = retrieve_pars(fit.sample, "reconstructions")

# Plot final signatures, exposures and reconstructions
cat("\n\nPlotting signatures, exposures and reconstructions:\n")
plot_spectrum(signatures.final, pdf_path=file.path(OUTPUT$PDF.DIR, OUTPUT$PDF.SIGS))
plot_spectrum(convert_signatures(signatures.final, opportunities_to="human-genome"),
              pdf_path=file.path(OUTPUT$PDF.DIR, OUTPUT$PDF.SIGS.NORM))
plot_exposures(fit.sample, cex_names=0.5, margin_bottom=9,
               pdf_path=file.path(OUTPUT$PDF.DIR, OUTPUT$PDF.EXP))
plot_reconstruction(fit.sample, pdf_path=file.path(OUTPUT$PDF.DIR, OUTPUT$PDF.REC))

# Build table of burden and exposures per sample
cat("\nBuilding burden and exposures table...\n")
burden.expos = data.frame("Sample" = rownames(sample.counts),
                          "Mutations" = rowSums(sample.counts),
                          "SBS1_exposure" = exposures.final$mean[, 1],
                          "SBSB_exposure" = exposures.final$mean[, 2],
                          "SBSC_exposure" = exposures.final$mean[, 3],
                          "Callable_genome_bp" = callable.size,
                          stringsAsFactors=F)

# Save signatures, exposures, reconstructions, and burden and exposures table
save(signatures.final, exposures.final, reconstructions.final, file=OUTPUT$SIGS.DATA)
save(burden.expos, file=OUTPUT$BURDEN.DATA)



# (2) Cross-species analysis of SBSB
# To explore variation in SBSB across species, the fit-extract model in sigfit is used
# to fit COSMIC SBS1, SBS18, SBS34 and extract one additional signature

# Normalise COSMIC signatures to genome-independent representation
sigs.fixed = convert_signatures(cosmic_signatures_v3[c("SBS1", "SBS18", "SBS34"), ],
                                opportunities_from="human-genome")

# For each species, analyse mutation counts per individual
cat("\nExtracting and plotting signature SBSB for each species:\n")
fit.sbsb = lapply(rownames(species.counts), function(species) {
    
    cat("\n", species, ":\n", sep="")
    species.idx = sample.info$SPECIES_NAME == species & !(sample.info$SAMPLE_NAME %in% exclude.list)
    indiv.ids = unique(sample.info$NORMAL_NAME[species.idx])
    indiv.counts = t(sapply(indiv.ids, function(id) {
        sample.ids = sample.info$SAMPLE_NAME[species.idx & sample.info$NORMAL_NAME == id]
        colSums(sample.counts[sample.ids, , drop=F])
    }))
    indiv.opps = t(sapply(indiv.ids, function(id) {
        sample.ids = sample.info$SAMPLE_NAME[species.idx & sample.info$NORMAL_NAME == id]
        normalise(colSums(sample.opps[sample.ids, , drop=F]))
    }))
    
    suppressWarnings(fit_extract_signatures(indiv.counts, sigs.fixed,
                                            num_extra_sigs=1, opportunities=indiv.opps,
                                            iter=ITER.2, warmup=WARMUP, seed=SEED, refresh = 0))
    cat(" Done\n")
})

sbsb.species = convert_signatures(rbind(as.numeric(signatures.final$mean[2, ]),
                                        t(sapply(fit.sbsb, function(f) {
                                            as.numeric(retrieve_pars(f, "signatures")$mean[4, ])
                                        }))), opportunities_to="human-genome")

sbsb.sim.sbs5 = apply(sbsb.species, 1, sigfit:::cosine_sim, cosmic_signatures_v3["SBS5", ])
sbsb.sim.sbs40 = apply(sbsb.species, 1, sigfit:::cosine_sim, cosmic_signatures_v3["SBS40", ])
rownames(sbsb.species) = paste0("SBSB as inferred from ",
                                c("all species (original)", gsub("_", " ", rownames(species.counts))),
                                " (cosine similarity to SBS5 = ", round(sbsb.sim.sbs5, 3),
                                ", SBS40 = ", round(sbsb.sim.sbs40, 3), ")")

plot_spectrum(sbsb.species, pdf_path=OUTPUT$PDF.SBSB, pdf_width=27)



# (3) Assess prevalence of colibactin/APOBEC-induced mutations in non-human crypts

# Load COSMIC v3.2 signatures (from sigfit v2.1)
data("cosmic_signatures_v3.2")
cat("\nAssessing colibactin and APOBEC prevalence...\n")

# Function: mutational signature fitting via expectation-maximization
# Based on original code by I. MARTINCORENA (Wellcome Sanger Institute)
emsig = function(counts, signatures, maxiter = 1e4, epsilon = 1e-6) {
    counts = as.numeric(counts)
    num_signatures = nrow(signatures)
    # EM algowith to estimate the signature contribution
    alpha = rep(1/num_signatures, num_signatures)  # uniform start
    for (iter in 1:maxiter) {
        contr = t(array(alpha, dim=c(num_signatures, 96))) * t(signatures)
        probs = contr / array(rowSums(contr), dim=dim(contr))
        probs = probs * counts
        old_alpha = alpha
        alpha = colSums(probs) / sum(probs)
        if (sum(abs(alpha-old_alpha)) < epsilon) {
            break
        }
    }
    # Log-likelihood
    pred_spectrum = rowSums(t(array(alpha, dim=c(num_signatures, 96))) * t(signatures))
    pred_spectrum = pred_spectrum / sum(pred_spectrum)
    LL = dmultinom(x=counts, prob=pred_spectrum, log=T)
    return(list(alpha=alpha, LL=LL))
}

# Remove human crypts from sample set
stopifnot(identical(rownames(sample.counts), rownames(sample.opps)))
human.idx = rownames(sample.counts) %in% sample.info$SAMPLE_NAME[sample.info$SPECIES_NAME == "human"]
sample.counts = sample.counts[!human.idx, ]
sample.opps = sample.opps[!human.idx, ]

# Colibactin (SBS88) test: Run EM fitting algorith for each sample, using 4 (SBS1,5,18,34)
# or 5 (+SBS88) signatures, and calculate colibactin SBS88 p-values via likelihood ratio test
sig.names = c("SBS1", "SBS5", "SBS18", "SBS34", "SBS88")
fit.colibactin = data.frame(Sample=rownames(sample.counts),
                            Mutations=rowSums(sample.counts),
                            SBS1=NA, SBS5=NA, SBS18=NA, SBS34=NA, SBS88=NA, Pval_SBS88=NA,
                            stringsAsFactors=F)

for (i in 1:nrow(sample.counts)) {
    # Convert signatures to sample-specific opportunities
    sigs = t(apply(cosmic_signatures_v3.2[sig.names, ], 1, function(s) {
        s1 = s / sigfit:::human_trinuc_freqs() * sample.opps[i, ]
        s1 / sum(s1)
    }))
    counts = sample.counts[i, ]
    out = emsig(counts, sigs)
    fit.colibactin[i, sig.names] = out$alpha  # Estimated exposures
    LL1 = out$LL                              # LogLik with all sigs
    LL0 = emsig(counts, sigs[1:4, ])$LL       # LogLik excluding SBS88
    fit.colibactin$Pval_SBS88[i] = 1 - pchisq(2*(LL1-LL0), 1)  # Likelihood ratio test
}
fit.colibactin$Qval_SBS88 = p.adjust(fit.colibactin$Pval_SBS88, method="BH")

# APOBEC (SBS2+13) test: Run EM fitting algorith for each sample, using 4 (SBS1,5,18,34)
# or 6 (+SBS2,13) signatures, and calculate APOBEC (SBS2+13) p-values via likelihood ratio test
sig.names = c("SBS1", "SBS5", "SBS18", "SBS34", "SBS2", "SBS13")
fit.apobec = data.frame(Sample=rownames(sample.counts),
                        Mutations=rowSums(sample.counts),
                        SBS2=NA, SBS13=NA, Pval_SBS2=NA, Pval_SBS13=NA,
                        stringsAsFactors=F)

for (i in 1:nrow(sample.counts)) {
    # Convert signatures to sample-specific opportunities
    sigs = t(apply(cosmic_signatures_v3.2[sig.names, ], 1, function(s) {
        s1 = s / sigfit:::human_trinuc_freqs() * sample.opps[i, ]
        s1 / sum(s1)
    }))
    counts = sample.counts[i, ]
    out2 = emsig(counts, sigs[-6, ])
    out13 = emsig(counts, sigs[-5, ])
    fit.apobec[i, "SBS2"] = out2$alpha["SBS2"]            # Estimated exposures
    fit.apobec[i, "SBS13"] = out13$alpha["SBS13"]
    LL1 = out2$LL                                         # LogLik with SBS2
    LL2 = out13$LL                                        # LogLik with SBS13
    LL0 = emsig(counts, sigs[1:4, ])$LL                   # LogLik without SBS2+13
    fit.apobec$Pval_SBS2[i] = 1 - pchisq(2*(LL1-LL0), 1)  # Likelihood ratio test
    fit.apobec$Pval_SBS13[i] = 1 - pchisq(2*(LL2-LL0), 1)
}
fit.apobec$Qval_SBS2 = p.adjust(fit.apobec$Pval_SBS2, method="BH")
fit.apobec$Qval_SBS13 = p.adjust(fit.apobec$Pval_SBS13, method="BH")

# Plot colibactin exposures for non-human samples
# (a sample is considered colibactin-positive if it has q < 0.05 for SBS88)
pdf(OUTPUT$PDF.SBS88, 12, 4)
suppressWarnings(par(mar=c(5.5, 3, 2, 0), mgp=c(1, -0.5, -1.3), tcl=-0.2))
cols = alpha(c("#228B22", "#FF8C00", "#9932CC", "grey70", "firebrick"), 0.8)
species = sample.info$SPECIES_NAME[match(fit.colibactin$Sample, sample.info$SAMPLE_NAME)]
sig.names = c("SBS1", "SBS5", "SBS18", "SBS34", "SBS88")
b = barplot(t(fit.colibactin[, sig.names]),
            ylab="Signature exposure", las=1, cex.names=1e-9, cex.lab=1.1, cex.axis=0.85,
            col=cols, border=NA, space=c(0, ifelse(species[-length(species)] == species[-1], 0, 1)))
x = c(0, b[species[-length(species)] != species[-1]], b[length(b)]) + 1
text(x=colMeans(rbind(x[-1], x[-length(x)])),
     y=par()$usr[3]-0.02*(par()$usr[4]-par()$usr[3]),
     labels=gsub("_", " ", toTitleCase(unique(species))), cex=0.9, srt=45, adj=1, xpd=TRUE)
text(b[fit.colibactin$Qval_SBS88 < 0.05], 1.02, "*", cex=1.4, font=2, xpd=NA)
legend("topright", sig.names, pch=22, pt.bg=cols, col=NA, pt.lwd=0.75, pt.cex=2,
       xpd=NA, horiz=T, bty="n", inset=c(0.03, -0.15), cex=1.0)
invisible(dev.off())

# Use Fisher's exact test to assess significance of depletion in colibactin/APOBEC relative
# to human crypts (Lee-Six et al. 2019). We test at the crypt level, given the difference
# in the numbers of crypts per individual between both datasets.
# A sample is considered colibactin-positive if it has q < 0.05 for SBS88.
# A sample is considered APOBEC-positive if it has q < 0.05 for both SBS2 and SBS13.
# Applying the EM+LRT method given above to the set of somatic SBS per crypt across 445
# human crypts sequenced by Lee-Six et al. (Nature 2019) gave the following contrasts:

#                 Colibactin- Colibactin+    APOBEC-     APOBEC+
# Human (Lee-Six)         353          92        436           9
# Non-human (CS)          179           1        179           1

# Redirect test output to text file
sink(OUTPUT$FISHER.TEST)

# Fisher's test for colibactin prevalence
cat("FISHER'S TEST FOR COLIBACTIN PREVALENCE IN HUMAN VS NON-HUMAN CRYPTS\n\n")
colibactin.counts = matrix(c(353, 92, table(fit.colibactin$Qval_SBS88 < 0.05)),
                           nrow=2, byrow=T,
                           dimnames=list(c("Human", "Non-human"), c("Colibactin-", "Colibactin+")))
print(colibactin.counts)
fisher.test(colibactin.counts)

# Fisher's test for APOBEC prevalence
cat("\n\nFISHER'S TEST FOR APOBEC PREVALENCE IN HUMAN VS NON-HUMAN CRYPTS\n\n")
apobec.counts = matrix(c(436, 9, table(fit.apobec$Qval_SBS2 < 0.05 & fit.apobec$Qval_SBS13 < 0.05)),
                       nrow=2, byrow=T,
                       dimnames=list(c("Human", "Non-human"), c("APOBEC-", "APOBEC+")))
print(apobec.counts)
fisher.test(apobec.counts)
sink()


cat("\nDone\n\n")
