# Analysis scripts for Cagan, Baez-Ortega et al., 2022
# Step 7: Calculation of somatic mutation burdens and rates per sample

# Adrian Baez-Ortega, 2020-22


# TO RUN THIS SCRIPT IN THE TERMINAL
# ----------------------------------
# Run the commands below in the terminal, replacing '/path/to/CrossSpecies2021' with the 
# path to the CrossSpecies2021 directory.
#
#    cd /path/to/CrossSpecies2021
#    Rscript scripts/7_Burdens.R


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
    
    # Path to sample data table
    SAMPLE.DATA = file.path("data", "original", "CrossSpecies_SampleData.txt"),
    
    # Path to list of samples to exclude
    EXCLUDE = file.path("data", "processed", "SamplesToExclude.txt"),
    
    # Path to list of samples to exclude from mtDNA analyses
    EXCLUDE.MT = file.path("data", "processed", "SamplesToExclude_mtDNA.txt"),
    
    # Path to table of total genome and mtDNA sizes
    GENOME.SIZE = file.path("data", "original", "GenomeSizes.txt"),
    
    # Path to table of callable mtDNA sizes
    MTDNA.SIZE = file.path("data", "original", "CallableGenome_mtDNA.txt"),
    
    # Path to table of genome and mtDNA coverage per sample
    COVERAGE = file.path("data", "original", "GenomeCoverage.txt"),
    
    # Path to file of mutation burden and exposures per sample
    BURDEN.EXPOS = file.path("data", "processed", "Burden_Exposures.RData"),
    
    # Path to table of lifespan estimates from Species 360
    LIFESPAN = file.path("data", "original", "Lifespan_Species360.txt"),
    
    # Path to table of AnAge life-history data
    ANAGE.DATA = file.path("data", "original", "AnAge_Data.txt"),
    
    # Path to file of filtered mtDNA variants
    MTDNA.FINAL = file.path("data", "original", "FinalCalls_mtDNA", "FinalCalls_mtDNA_Merged.txt"),
    
    # Generic path to files of filtered indels
    INDELS.FINAL = file.path("data", "original", "FinalCalls_Indels", "${SPECIES}",
                             "FinalCalls_Indels_${SAMPLE}.txt")
)

# Output file paths
OUTPUT = list(
    
    # Path to PDF of mutation burdens, rates and ELBs per sample
    PDF = file.path("output", "Burden_Rate_ELB.pdf"),
    
    # Path to table of mutation burdens, rates and ELBs per sample
    TABLE = file.path("output", "Burden_Rate_ELB.txt"),
    
    # Path to RData file of mutation burdens and rates
    DATA = file.path("data", "processed", "Burdens_Rates.RData")

)


# Function: multiple pattern replacement (gsub)
multi.gsub = function(patterns, replacements, x, ...) {
    stopifnot(length(patterns) == length(replacements))
    for (i in 1:length(patterns)) {
        x = gsub(patterns[i], replacements[i], x, ...)
    }
    x
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
PACKAGES = c("scales", "tools")
cat("Loading packages:", paste(PACKAGES, collapse=", "), "\n")
for (package in PACKAGES) {
    suppressWarnings(suppressPackageStartupMessages(library(package, character.only=TRUE,
                                                            quietly=TRUE)))
}

cat("Loading data...\n")
load(INPUT$BURDEN.EXPOS)
exclude.list = as.character(as.matrix(read.table(INPUT$EXCLUDE)))
exclude.list.mt = as.character(as.matrix(read.table(INPUT$EXCLUDE.MT)))
mt.vars = read.table(INPUT$MTDNA.FINAL, sep="\t", header=T, as.is=T)
genome.size = read.table(INPUT$GENOME.SIZE, sep="\t", header=T, as.is=T)
mtdna.size = read.table(INPUT$MTDNA.SIZE, sep="\t", header=T, as.is=T)
coverage = read.table(INPUT$COVERAGE, sep="\t", header=T, as.is=T)
lifespan = read.table(INPUT$LIFESPAN, sep="\t", header=T, as.is=T)
anage.data = read.table(INPUT$ANAGE.DATA, sep="\t", header=T, as.is=T, check.names=F)
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
sample.data = read.table(INPUT$SAMPLE.DATA, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
stopifnot(!any(duplicated(sample.data$Sample)))
stopifnot(!any(sample.data$Sample %in% exclude.list))
cat("Loaded\n\n")


# Integrate Lifespan_80 values into sample table (for regressions)
sample.data$Lifespan80 = lifespan$Lifespan_80[match(sample.data$Species, lifespan$Species)]
sample.data$Inverse_Lifespan80 = 1 / sample.data$Lifespan80

# Integrate body mass, litter size and metabolic rate (for regressions)
idx = match(sample.data$Species, anage.data$Name)
sample.data$Mass = anage.data$`Adult weight (g)`[idx]
sample.data$Log10_Mass = log10(sample.data$Mass)
sample.data$Litter_size = anage.data$`Litter/Clutch size`[idx]
sample.data$Metabolic_rate = anage.data$`Metabolic rate (W)`[idx]
sample.data$Log10_Metabolic_rate = log10(sample.data$Metabolic_rate)


# Integrate mutation (SBS) burden into sample table and calculate corrected burdens
# (substitutions per genome; obtained by dividing burden by the sample's callable
# genome size and multiplying by the species' total genome size)
cat("Calculating corrected mutation burdens and rates per sample...\n")
sample.data = cbind(sample.data,
                    burden.expos[match(sample.data$Sample, burden.expos$Sample), -1])
sample.data$Total_genome_bp = genome.size$Total_genome_size[match(sample.data$Species,
                                                                  genome.size$Species)]
sample.data$Mut_burden = sample.data$Mutations /
    sample.data$Callable_genome_bp * sample.data$Total_genome_bp

# Calculate corrected indel burdens
sample.data$Indels = sapply(1:nrow(sample.data), function(i) {
    nrow(read.table(multi.gsub(c("${SAMPLE}", "${SPECIES}"),
                               c(sample.data$Sample[i], sample.data$Species[i]),
                               INPUT$INDELS.FINAL, fixed=T), sep="\t", header=T))
})
sample.data$Indel_burden = sample.data$Indels /
    sample.data$Callable_genome_bp * sample.data$Total_genome_bp

# Calculate corrected mutational signature burdens
sample.data$SBS1_burden = sample.data$SBS1_exposure * sample.data$Mut_burden
sample.data$SBSB_burden = sample.data$SBSB_exposure * sample.data$Mut_burden
sample.data$SBSC_burden = sample.data$SBSC_exposure * sample.data$Mut_burden

# Calculate corrected mtDNA mutation burdens
# (raw mtDNA burden in a sample is given by the sum of mtDNA VAFs;
# 'coding' mtDNA is defined as all mtDNA except the hypervariable D-loop)
sample.data$mtDNA_burden = sapply(1:nrow(sample.data), function(i) {
    if (sample.data$Sample[i] %in% exclude.list.mt) {
        NA
    } else {
        sum(mt.vars$VAF[mt.vars$Sample == sample.data$Sample[i]]) /
            mtdna.size$Callable_mtDNA_coding[match(sample.data$Sample[i], mtdna.size$Sample)] *
            genome.size$Coding_mtDNA_size[match(sample.data$Species[i], genome.size$Species)]
    }
})

# Calculate mutation rates (mutations per genome per year)
for (mut in c("Mut", "Indel", "SBS1", "SBSB", "SBSC", "mtDNA")) {
    sample.data[paste0(mut, "_rate")] = sample.data[paste0(mut, "_burden")] / sample.data$Age_years
}

# Calculate expected end-of-lifespan burdens (ELBs; mutations per genome)
for (mut in c("Mut", "Indel", "SBS1", "SBSB", "SBSC", "mtDNA")) {
    sample.data[paste0(mut, "_ELB")] = sample.data[paste0(mut, "_rate")] * sample.data$Lifespan80
}

# Calculate mtDNA copy number
idx = match(sample.data$Sample, coverage$Sample)
sample.data$mtDNA_CN = coverage$mtDNA_coverage[idx] / coverage$Genome_coverage[idx] * 2


# Build tables of mean burdens and rates per individual and per species
cat("Calculating mean burdens and rates per individual and per species...\n")
indiv.data = sample.data[match(unique(sample.data$Individual), sample.data$Individual),
                         c("Individual", "Species", "Age_years", "Lifespan80", "Inverse_Lifespan80",
                           "Mass", "Log10_Mass", "Litter_size", "Metabolic_rate",
                           "Log10_Metabolic_rate")]
for (mut in c("Mut", "Indel", "SBS1", "SBSB", "SBSC", "mtDNA")) {
    burden = paste0(mut, "_burden")
    rate = paste0(mut, "_rate")
    indiv.data[, burden] = sapply(indiv.data$Individual, function(id) {
        mean(sample.data[sample.data$Individual == id, burden], na.rm=T)
    })
    indiv.data[, rate] = indiv.data[, burden] / indiv.data$Age_years
}

species.data = sample.data[match(unique(sample.data$Species), sample.data$Species),
                           c("Species", "Lifespan80", "Inverse_Lifespan80",
                             "Mass", "Log10_Mass", "Litter_size", "Metabolic_rate",
                             "Log10_Metabolic_rate")]
anage.idx = match(species.data$Species, anage.data$Name)
species.data$Latin_name = gsub("familiaris", "lupus",
                               paste0(anage.data$Genus, "_", anage.data$Species)[anage.idx])
for (mut in c("Mut", "Indel", "SBS1", "SBSB", "SBSC", "mtDNA")) {
    burden = paste0(mut, "_burden")
    rate = paste0(mut, "_rate")
    species.data[, burden] = sapply(species.data$Species, function(sp) {
        mean(indiv.data[indiv.data$Species == sp, burden], na.rm=T)
    })
    species.data[, rate] = sapply(species.data$Species, function(sp) {
        mean(indiv.data[indiv.data$Species == sp, rate], na.rm=T)
    })
}
rownames(sample.data) = rownames(indiv.data) = rownames(species.data) = NULL


# Round all mutation burdens (except for mtDNA)
for (mut in c("Mut", "Indel", "SBS1", "SBSB", "SBSC")) {
    species.data[paste0(mut, "_burden")] = round(species.data[paste0(mut, "_burden")])
    indiv.data[paste0(mut, "_burden")] = round(indiv.data[paste0(mut, "_burden")])
    sample.data[paste0(mut, "_burden")] = round(sample.data[paste0(mut, "_burden")])
    sample.data[paste0(mut, "_ELB")] = round(sample.data[paste0(mut, "_ELB")])
}


# Output and save burdens tables
cat("Saving and plotting mutation burdens and rates per sample...\n")
write.table(sample.data[, c(1:12, 17:18, 38, 13, 20, 14:16, 19, 21:37)],
            file=OUTPUT$TABLE, sep="\t", quote=F, row.names=F)
save(sample.data, indiv.data, species.data, file=OUTPUT$DATA)


# Plot mutation burdens, rates, ELBs and mtDNA CN per sample
cairo_pdf(OUTPUT$PDF, 14, 4, onefile=T)
par(mar=c(3, 3.0, 2, 0.75), tck=-0.02, mgp=c(1.7, 0.33, 0))
cols = c("grey35", "grey65")

# Sort species, individuals and samples by burden (for plotting)
sample.idx = sample.col = NULL
for (species in species.data$Species[order(species.data$Mut_burden)]) {
    idx = indiv.data$Species == species
    indivs = indiv.data$Individual[idx][order(indiv.data$Mut_burden[idx])]
    for (id in indivs) {
        idx = sample.data$Individual == id
        sample.col = c(sample.col, rep(cols[1], sum(idx)))
        sample.idx = c(sample.idx, which(idx)[order(sample.data$Mut_burden[idx])])
        cols = rev(cols)
    }
    sample.col = c(sample.col, NA, NA)
    sample.idx = c(sample.idx, NA, NA)
}
sample.col = sample.col[1:(length(sample.col)-2)]
sample.idx = sample.idx[1:(length(sample.idx)-2)]
species.lab = multi.gsub(c("Naked", "mole_rat", "Ring_tailed", "Harbour", "Colobus", "_"),
                         c("N", "mole-rat", "RT", "H", "BW colobus", " "),
                         toTitleCase(species.data$Species[order(species.data$Mut_burden)]))

# Plot mutation burden per sample
plot(sample.data$Mut_burden[sample.idx], type="n", xlab="", ylab="Mutations per genome",
     ylim=c(0, max(sample.data$Mut_burden) * 1.05), xlim=c(-0.5, length(sample.idx)+1.5),
     xaxt="n", yaxt="n", yaxs="i", xaxs="i", cex.lab=1.1, main="Mutation burden per sample",
     panel.first={
         abline(h=seq(1000, 4000, 1000), col="grey90", lwd=0.75);
         abline(v=colMeans(matrix(which(is.na(sample.col)), nrow=2)))
         segments(x0=seq_along(sample.idx), y0=0,
                  y1=sample.data$Mut_burden[sample.idx], col=sample.col, lwd=0.5);
     })
axis(2, cex.axis=0.8)
points(sample.data$Mut_burden[sample.idx], pch=16, col=sample.col, cex=1)
points(sample.data$Indel_burden[sample.idx], pch=18, cex=0.9, col=alpha(sample.col, 0.85))
legend("topleft", legend=c("Substitutions", "Indels"), horiz=T, inset=c(0, -0.11), xpd=T, bty="n",
       pch=c(16, 18), col=c(sample.col[1], alpha(sample.col[1], 0.85)), pt.cex=1.3, cex=0.9)
text(x=colMeans(matrix(c(0, which(is.na(sample.col)), length(sample.col)+1), nrow=2)),
     y=par()$usr[3] - 0.022 * (par()$usr[4] - par()$usr[3]),
     labels=species.lab, cex=0.85, srt=35, adj=1, xpd=TRUE)

# Plot mutation rate per sample
plot(sample.data$Mut_rate[sample.idx], type="n", xlab="", ylab="Mutations per genome per year",
     ylim=c(0, max(sample.data$Mut_rate, na.rm=T) * 1.0244), xlim=c(-0.5, length(sample.idx)+1.5),
     xaxt="n", yaxt="n", yaxs="i", xaxs="i", cex.lab=1.1, main="Mutation rate per sample",
     panel.first={
         abline(h=seq(300, 900, 300), col="grey90", lwd=0.75);
         abline(v=colMeans(matrix(which(is.na(sample.col)), nrow=2)))
         segments(x0=seq_along(sample.idx), y0=0,
                  y1=sample.data$Mut_rate[sample.idx], col=sample.col, lwd=0.5);
     })
axis(2, cex.axis=0.8, at=seq(0, 1200, 300))
points(sample.data$Mut_rate[sample.idx], pch=16, col=sample.col, cex=1.0)
points(sample.data$Indel_rate[sample.idx], pch=18, cex=0.9, col=alpha(sample.col, 0.85))
legend("topleft", legend=c("Substitutions", "Indels"), horiz=T, inset=c(0, -0.11), xpd=T, bty="n",
       pch=c(16, 18), col=c(sample.col[1], alpha(sample.col[1], 0.85)), pt.cex=1.3, cex=0.9)
text(x=colMeans(matrix(c(0, which(is.na(sample.col)), length(sample.col)+1), nrow=2)),
     y=par()$usr[3] - 0.022 * (par()$usr[4] - par()$usr[3]),
     labels=species.lab, cex=0.85, srt=35, adj=1, xpd=TRUE)

# Plot ELB per sample
plot(sample.data$Mut_ELB[sample.idx], type="n", xlab="", ylab="End-of-lifespan burden (ELB)",
     ylim=c(0, max(sample.data$Mut_ELB, na.rm=T) * 1.05), xlim=c(-0.5, length(sample.idx)+1.5),
     xaxt="n", yaxt="n", yaxs="i", xaxs="i", cex.lab=1.1, main="ELB per sample",
     panel.first={
         abline(h=seq(2000, 8000, 2000), col="grey90", lwd=0.75);
         abline(v=colMeans(matrix(which(is.na(sample.col)), nrow=2)))
         segments(x0=seq_along(sample.idx), y0=0,
                  y1=sample.data$Mut_ELB[sample.idx], col=sample.col, lwd=0.5);
     })
axis(2, cex.axis=0.8, at=seq(0, 8000, 2000))
points(sample.data$Mut_ELB[sample.idx], pch=16, col=sample.col, cex=1.0)
points(sample.data$Indel_ELB[sample.idx], pch=18, cex=0.9, col=alpha(sample.col, 0.85))
legend("topleft", legend=c("Substitutions", "Indels"), horiz=T, inset=c(0, -0.11), xpd=T, bty="n",
       pch=c(16, 18), col=c(sample.col[1], alpha(sample.col[1], 0.85)), pt.cex=1.3, cex=0.9)
text(x=colMeans(matrix(c(0, which(is.na(sample.col)), length(sample.col)+1), nrow=2)),
     y=par()$usr[3] - 0.022 * (par()$usr[4] - par()$usr[3]),
     labels=species.lab, cex=0.85, srt=35, adj=1, xpd=TRUE)

# Plot mtDNA mutation burden per sample
plot(sample.data$mtDNA_burden[sample.idx], type="n", xlab="", ylab="Mutations per mtDNA",
     ylim=c(0, max(sample.data$mtDNA_burden, na.rm=T) * 1.05), xlim=c(-0.5, length(sample.idx)+1.5),
     xaxt="n", yaxt="n", yaxs="i", xaxs="i", cex.lab=1.1, main="mtDNA mutation burden per sample",
     panel.first={
         abline(h=seq(0.5, 1.5, 0.5), col="grey90", lwd=0.75);
         abline(v=colMeans(matrix(which(is.na(sample.col)), nrow=2)))
         segments(x0=seq_along(sample.idx), y0=0,
                  y1=sample.data$mtDNA_burden[sample.idx], col=sample.col, lwd=0.5);
     })
axis(2, cex.axis=0.8, at=seq(0, 1.5, 0.5))
points(sample.data$mtDNA_burden[sample.idx], pch=16, col=sample.col, cex=1.0)
text(x=colMeans(matrix(c(0, which(is.na(sample.col)), length(sample.col)+1), nrow=2)),
     y=par()$usr[3] - 0.022 * (par()$usr[4] - par()$usr[3]),
     labels=species.lab, cex=0.85, srt=35, adj=1, xpd=TRUE)

# Plot mtDNA copy number per sample
plot(sample.data$mtDNA_CN[sample.idx], type="n", xlab="", ylab="mtDNA copy number",
     ylim=c(0, max(sample.data$mtDNA_CN, na.rm=T) * 1.05), xlim=c(-0.5, length(sample.idx)+1.5),
     xaxt="n", yaxt="n", yaxs="i", xaxs="i", cex.lab=1.1, main="mtDNA copy number per sample",
     panel.first={
         abline(h=seq(500, 2000, 500), col="grey90", lwd=0.75);
         abline(v=colMeans(matrix(which(is.na(sample.col)), nrow=2)))
         segments(x0=seq_along(sample.idx), y0=0,
                  y1=sample.data$mtDNA_CN[sample.idx], col=sample.col, lwd=0.5);
     })
axis(2, cex.axis=0.8, at=seq(0, 2000, 500))
points(sample.data$mtDNA_CN[sample.idx], pch=16, col=sample.col, cex=1.0)
text(x=colMeans(matrix(c(0, which(is.na(sample.col)), length(sample.col)+1), nrow=2)),
     y=par()$usr[3] - 0.022 * (par()$usr[4] - par()$usr[3]),
     labels=species.lab, cex=0.85, srt=35, adj=1, xpd=TRUE)

invisible(dev.off())


cat("\nDone\n\n")
