# Analysis scripts for Cagan, Baez-Ortega et al., 2021
# Step 2: Mutational spectra of single-base substitutions (SBS)

# Adrian Baez-Ortega, 2020-21


# TO RUN THIS SCRIPT IN THE TERMINAL
# ----------------------------------
# Run the commands below in the terminal, replacing '/path/to/CrossSpecies2021' with the 
# path to the CrossSpecies2021 directory.
#
#    cd /path/to/CrossSpecies2021
#    Rscript scripts/2_Spectra_SBS.R
#
# NB. If the script triggers memory allocation errors, the amount of virtual memory
# assigned to R can be increased using:
#
#    env R_MAX_VSIZE=100Gb Rscript scripts/2_Spectra_SBS.R


# TO RUN THIS SCRIPT IN RSTUDIO
# -----------------------------
# Before starting, run the line below in RStudio, replacing '/path/to/CrossSpecies2021' with 
# the path to the CrossSpecies2021 directory.
#
#    setwd("path/to/CrossSpecies2021")
#
# (In the event that there is insufficient memory, see non-interactive alternative above)


# If the paths to the input or output files differ from the ones
# used in this script, they may be updated by modifying the lines below.


# Input file paths
INPUT = list(
    
    # Path to project info table
    SAMPLE.INFO = file.path("data", "original", "CrossSpecies_ProjectInfo.txt"),
    
    # Path to list of samples to exclude
    EXCLUDE = file.path("data", "processed", "SamplesToExclude.txt"),
    
    # Path to callable genome file
    CALLABLE.GNM = file.path("data", "original", "CallableGenome.RData"),
    
    # Generic path to reference genome files
    GENOME = file.path("data", "original", "RefGenomes", "${SPECIES}_genome.fa.gz"),
    
    # Generic path to files of filtered substitutions
    # (see Methods for details on variant calling and filtering)
    SUBS.FINAL = file.path("data", "original", "FinalCalls_Subs", "${SPECIES}",
                           "FinalCalls_Subs_${SAMPLE}.txt")
)

# Output file paths
OUTPUT = list(
    
    # Paths to substitution spectra per sample and per species
    PDF.SAMPLE = file.path("output", "Spectra_Subs_Sample.pdf"),
    PDF.SPECIES = file.path("output", "Spectra_Subs_Species.pdf"),
    
    # Path to RData file of mutational catalogues and opportunities
    COUNTS = file.path("data", "processed", "Catalogues_Opportunities.RData")

)


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
PACKAGES = c("Biostrings", "GenomicRanges", "sigfit", "tools")
cat("Loading packages:", paste(PACKAGES, collapse=", "), "\n")
for (package in PACKAGES) {
    suppressPackageStartupMessages(library(package, character.only=TRUE, quietly=TRUE))
}

cat("Loading data...\n")
load(INPUT$CALLABLE.GNM)
exclude.list = as.character(as.matrix(read.table(INPUT$EXCLUDE)))
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
cat("Loaded\n")


# Initialise tables of variants and mutational opportunities
species.vars = sample.vars = species.opps = sample.opps = callable.size = NULL


for (species in unique(sample.info$SPECIES_NAME)) {
    cat("\nProcessing species:", species, "\n")
    
    # Load reference genome
    cat("Loading reference genome\n")
    ref.name = sample.info$REFERENCE_GENOME[sample.info$SPECIES_NAME == species][1]
    genome = readDNAStringSet(gsub("${SPECIES}", species, INPUT$GENOME, fixed=T),
                              format="fasta", use.names=TRUE)
    names(genome) = sapply(strsplit(names(genome), " "), `[`, 1)
    
    # Process variants and mutational opportunities per sample
    samples = sample.info$SAMPLE_NAME[sample.info$SPECIES_NAME == species &
                                          !(sample.info$SAMPLE_NAME %in% exclude.list)]
    for (sample.id in samples) {
        
        cat("Processing sample", sample.id, "\n")
        var.path = gsub("${SAMPLE}", sample.id,
                       gsub("${SPECIES}", species,
                            INPUT$SUBS.FINAL, fixed=T), fixed=T)
        if (!file.exists(var.path)) {
            stop("ERROR: File", var.path, "not found.\n")
        }
        
        vars = read.table(var.path, sep="\t", header=T, check.names=F, as.is=T, comment.char="",
                          colClasses=c("Chr"="character", "Ref"="character", "Alt"="character"))
        
        # Retrieve trinucleotide contexts
        context = as.character(padAndClip(genome[vars$Chr],
                                          IRanges(vars$Start - 1, vars$Start + 1),
                                          Lpadding.letter=".", Rpadding.letter="."))
        stopifnot(identical(as.character(substr(context, 2, 2)), vars$Ref))
        
        # Append variants to variant tables
        sample.vars = rbind(sample.vars,
                            cbind(Sample = paste("Filtered variants in", sample.id),
                                  Ref = vars$Ref,
                                  Alt = vars$Alt,
                                  Context = context))
        
        species.vars = rbind(species.vars,
                             cbind(Sample = paste("Filtered variants in",
                                                  gsub("_", " ", toTitleCase(species))),
                                   Ref = vars$Ref,
                                   Alt = vars$Alt,
                                   Context = context))
        
        # Select callable genome regions
        regions = callable.genome[[sample.id]]
        genome.regions = padAndClip(genome[seqnames(regions)],
                                    IRanges(start(regions), end(regions)),
                                    Lpadding.letter=".", Rpadding.letter=".")
        callable.size = c(callable.size, sum(width(regions)))
        
        # Calculate trinucleotide frequencies and mutational opportunities
        freqs = colSums(trinucleotideFrequency(genome.regions))
        freqs.pyr = sapply(which(substr(names(freqs), 2, 2) %in% c("C", "T")), function(i) {
            rcomp = rev.comp(names(freqs)[i])
            freqs[i] + freqs[rcomp]
        })
        sample.opps = rbind(sample.opps,
                            normalise(freqs.pyr[substr(sigfit:::mut_types(), 1, 3)]))
        
        rownames(sample.opps)[nrow(sample.opps)] =
            names(callable.size)[length(callable.size)] = sample.id
    }
    
    # Combine sample-level opportunities into species-level opportunities
    species.opps = rbind(species.opps,
                         normalise(colSums(sample.opps[samples, ])))
    rownames(species.opps)[nrow(species.opps)] = species
    
    rm(genome, vars, context, regions, genome.regions)
    invisible(gc())
}


# Build mutational catalogues per sample and per species
cat("\nBuilding mutational catalogues...\n")
sample.counts = suppressWarnings(build_catalogues(sample.vars))
species.counts = suppressWarnings(build_catalogues(species.vars))
colnames(sample.opps) = colnames(species.opps) = colnames(sample.counts)
stopifnot(identical(gsub("Filtered variants in ", "", rownames(sample.counts)),
                    rownames(sample.opps)))
stopifnot(identical(tolower(gsub(" ", "_",
                                 gsub("Filtered variants in ", "", rownames(species.counts)))),
                    rownames(species.opps)))

# Plot mutational spectra
cat("Plotting mutational spectra...\n")
plot_spectrum(sample.counts, OUTPUT$PDF.SAMPLE)
plot_spectrum(species.counts, OUTPUT$PDF.SPECIES)


# Save mutational catalogues and opportunities
save(sample.counts, species.counts, sample.opps, species.opps, callable.size, file=OUTPUT$COUNTS)


cat("\nDone\n\n")
