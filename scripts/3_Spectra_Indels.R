# Analysis scripts for Cagan, Baez-Ortega et al., 2022
# Step 3: Mutational spectra of insertions/deletions (indels)

# Adrian Baez-Ortega, 2020-22


# TO RUN THIS SCRIPT IN THE TERMINAL
# ----------------------------------
# Run the commands below in the terminal, replacing '/path/to/CrossSpecies2021' with the 
# path to the CrossSpecies2021 directory.
#
#    cd /path/to/CrossSpecies2021
#    Rscript scripts/3_Spectra_Indels.R


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
    
    # Path to Indelwald code
    INDELWALD = file.path("scripts", "Indelwald.R"),
    
    # Path to project info table
    SAMPLE.INFO = file.path("data", "original", "CrossSpecies_ProjectInfo.txt"),
    
    # Path to list of samples to exclude
    EXCLUDE = file.path("data", "processed", "SamplesToExclude.txt"),
    
    # Generic path to reference genome files
    GENOME = file.path("data", "original", "RefGenomes", "${SPECIES}_genome.fa.gz"),
    
    # Generic path to files of filtered indels
    INDELS.FINAL = file.path("data", "original", "FinalCalls_Indels", "${SPECIES}",
                             "FinalCalls_Indels_${SAMPLE}.txt")
)

# Output file paths
OUTPUT = list(
    
    # PDFs of spectra per sample and per species
    PDF.SAMPLE = file.path("output", "Spectra_Indels_Sample.pdf"),
    PDF.SPCS = file.path("output", "Spectra_Indels_Species.pdf")
    
)


# Function: indel spectrum plotting
plot.spectrum = function(spectra) {
    COLS = c(rep(c("#FCBD6F", "#FE8002", "#AFDC8A", "#36A02E", "#FCC9B4", "#FB896A", "#F04432",
                   "#BB191A", "#CFE0F1", "#93C3DE", "#4A97C8", "#1764AA"), each=6), "#E1E1EE",
             rep("#B5B5D7", 2), rep("#8582BC", 3), rep("#62409A", 5))
    colnames(spectra) = c(rep(c(paste0(1:5, "  "), "6+"), 2), rep(c(paste0(0:4, "  "), "5+"), 2),
                          rep(c(paste0(1:5, "  "), "6+"), 4), rep(c(paste0(0:4, "  "), "5+"), 4),
                          paste0(c(1, 1:2, 1:3, 1:4), "  "), "5+")
    for (i in 1:nrow(spectra)) {
        bars = barplot(spectra[i, ], col=COLS, mgp=c(3, 0.3, 0), border="white", yaxt="n",
                       xaxs="i", xlim=c(-0.105, 99.9), ylim=c(0, max(spectra[i, ]) * 1.095),
                       cex.names=1.5, las=2, adj=0.5)
        axis(side=2, cex.axis=1.9, lwd=2)
        mtext("Indels", side=2, cex=2.4, line=3.5)
        rect(xleft = bars[match(unique(COLS), COLS)] - 0.5, 
             xright = bars[c(match(unique(COLS), COLS)[-1] - 1, length(COLS))] + 0.5,
             ybottom = 0.94 * max(spectra[i, ]) * 1.095, ytop = max(spectra[i, ]) * 1.095,
             col = unique(COLS), border = "white")
        mtext(c("1-bp deletion", "1-bp insertion", ">1-bp deletion at repeat\n(deletion length)",
                ">1-bp insertion at repeat\n(insertion length)", "Microhomology\n(deletion length)"),
              at=c(mean(bars[c(1, 12)]), mean(bars[c(13, 24)]), mean(bars[c(25, 48)]),
                   mean(bars[c(49, 72)]), mean(bars[c(73, 83)])),
              side=3, line=2.5, cex=2.2)
        mtext(c(rep(c("C", "T"), 2), rep(c("2", "3", "4", "5+"), 3)), side=3, line=0.35, cex=2.1,
              at=rowMeans(cbind(bars[match(unique(COLS), COLS)] - 0.5,
                                bars[c(match(unique(COLS), COLS)[-1] - 1, length(COLS))] + 0.5)))
        mtext(c(rep(c("Homopolymer length", "Number of repeat units"), each=2),
                "Microhomology length"),
              side=1, line=2.8, cex=c(rep(1.95, 4), 1.85),
              at=c(mean(bars[c(1, 12)]), mean(bars[c(13, 24)]), mean(bars[c(25, 48)]),
                   mean(bars[c(49, 72)]), mean(bars[c(73, 83)]) + 0.2))
        title(paste0(rownames(spectra)[i], " (", prettyNum(sum(spectra[i, ]), big.mark=","),
                     " indels)"), line=8.5, cex.main=2.8)
        box(lwd = 2)
    }
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


# Load packages
PACKAGES = c("Biostrings", "tools")
cat("Loading packages:", paste(PACKAGES, collapse=", "), "\n")
for (package in PACKAGES) {
    suppressPackageStartupMessages(library(package, character.only=TRUE, quietly=TRUE))
}


# Load function 'indel.spectrum' from Indelwald, by MAX STAMMNITZ
# (https://github.com/MaximilianStammnitz/Indelwald)
source(INPUT$INDELWALD)


cat("Loading data...\n")
exclude.list = as.character(as.matrix(read.table(INPUT$EXCLUDE)))
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
cat("Loaded\n")


# Produce spectra for each species
species.spectra = sample.spectra = NULL
for (species in unique(sample.info$SPECIES_NAME)) {

    cat("\nProcessing species:", species, "\n")
    samples = sample.info$SAMPLE_NAME[sample.info$SPECIES_NAME == species &
                                          !(sample.info$SAMPLE_NAME %in% exclude.list)]
    cat("Loading reference genome\n")
    ref.name = sample.info$REFERENCE_GENOME[sample.info$SPECIES_NAME == species][1]
    genome = readDNAStringSet(gsub("${SPECIES}", species, INPUT$GENOME, fixed=T),
                              format="fasta", use.names=TRUE)
    names(genome) = sapply(strsplit(names(genome), " "), `[`, 1)
    
    # Produce spectra per sample
    species.vars = NULL
    spectra = t(sapply(samples, function(sample.id) {
        
        cat("Processing sample", sample.id, "\n")
        vars = read.table(gsub("${SAMPLE}", sample.id,
                               gsub("${SPECIES}", species,
                                    INPUT$INDELS.FINAL, fixed=T), fixed=T),
                          sep="\t", header=T, check.names=F, as.is=T, comment.char="",
                          col.names=c("CHROM", "POS", "REF", "ALT",	"", "", "", "", "", "", "", ""),
                          colClasses=c("CHROM"="character", "REF"="character", "ALT"="character"))
        
        if (nrow(vars) > 0) {
        
            # Add indels to species table
            species.vars <<- rbind(species.vars, vars[, c("CHROM", "POS", "REF", "ALT")])

            # Produce indel spectrum for sample
            indel.spectrum(vars[, c("CHROM", "POS", "REF", "ALT")], genome)
        
        } else {
            rep(0, 83)
        }
    }))
    
    sample.spectra = rbind(sample.spectra, spectra)
    species.spectra = rbind(species.spectra, indel.spectrum(species.vars, genome))
    rownames(species.spectra)[nrow(species.spectra)] = gsub("_", " ", toTitleCase(species))
    rm(genome); invisible(gc())
}


# Plot spectra per sample and per species
cat("\nPlotting mutational spectra...\n")
cairo_pdf(OUTPUT$PDF.SAMPLE, 24, 8, onefile=T)
par(mar = c(5.5, 7, 12.5, 2))
plot.spectrum(sample.spectra)
invisible(dev.off())

cairo_pdf(OUTPUT$PDF.SPCS, 24, 8, onefile=T)
par(mar = c(5.5, 7, 12.5, 2))
plot.spectrum(species.spectra)
invisible(dev.off())


cat("\nDone\n\n")
