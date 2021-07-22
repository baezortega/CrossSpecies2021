# Analysis scripts for Cagan, Baez-Ortega et al., 2021
# Step 3: Mutational spectra of insertions/deletions (indels)

# Adrian Baez-Ortega, 2020-21


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
    
    # Path to project info table
    SAMPLE.INFO = file.path("data", "original", "CrossSpecies_ProjectInfo.txt"),
    
    # Path to list of samples to exclude
    EXCLUDE = file.path("data", "processed", "SamplesToExclude.txt"),
    
    # Generic path to files of filtered indels
    INDELS.FINAL = file.path("data", "original", "FinalCalls_Indels", "${SPECIES}",
                             "FinalCalls_Indels_${SAMPLE}.txt")
)

# Output file paths
OUTPUT = list(
    
    # Path to temporary directory for SigProfiler
    TMP.DIR = "SigProfiler_tmp",
    
    # Paths to directories of spectra per sample and per species
    DIR.SAMPLE = file.path("output", "Spectra_Indels_Sample"),
    DIR.SPCS = file.path("output", "Spectra_Indels_Species"),
    
    # Generic names for PDFs of spectra per sample and per species
    PDF.SAMPLE = "Spectra_Indels_Sample_${SPECIES}.pdf",
    PDF.SPCS = "Spectra_Indels_Species_${SPECIES}.pdf"
    
)


# SigProfiler reference genome names
REF.NAMES = c("human"="GRCh37", "mouse"="mm10", "rat"="rn6", "dog"="dog")


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
PACKAGES = c("SigProfilerMatrixGeneratorR", "tools")
cat("Loading packages:", paste(PACKAGES, collapse=", "), "\n")
for (package in PACKAGES) {
    suppressPackageStartupMessages(library(package, character.only=TRUE, quietly=TRUE))
}


cat("Loading data...\n")
exclude.list = as.character(as.matrix(read.table(INPUT$EXCLUDE)))
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
cat("Loaded\n")


# Create output directories
dir.create(OUTPUT$DIR.SAMPLE, showWarnings=F)
dir.create(OUTPUT$DIR.SPCS, showWarnings=F)


for (species in names(REF.NAMES)) {
    
    cat("\nProcessing species:", species, "\n")
    species.vars = NULL
    
    # Create temporary folder for SigProfiler
    unlink(OUTPUT$TMP.DIR, recursive=T)
    dir.create(OUTPUT$TMP.DIR, showWarnings=F)
    
    samples = sample.info$SAMPLE_NAME[sample.info$SPECIES_NAME == species &
                                          !(sample.info$SAMPLE_NAME %in% exclude.list)]
    for (sample.id in samples) {
        
        cat("Processing sample", sample.id, "\n")
        vars = read.table(gsub("${SAMPLE}", sample.id,
                               gsub("${SPECIES}", species,
                                    INPUT$INDELS.FINAL, fixed=T), fixed=T),
                          sep="\t", header=T, check.names=F, as.is=T, comment.char="",
                          colClasses=c("Chrom"="character", "Ref"="character", "Alt"="character"))
        
        if (nrow(vars) > 0) {
            
            # Add indels to species table
            species.vars = rbind(species.vars, vars[, 1:4])
            
            # Output indels to VCF (for SigProfiler)
            write.table(cbind(vars[, 1:2], ".", vars[, 3:4]),
                        file=file.path(OUTPUT$TMP.DIR,
                                       paste0("Filtered_indels_in_", sample.id, ".vcf")),
                        sep="\t", quote=F, col.names=F, row.names=F)
        }
    }
    
    # Plot sample spectra
    invisible(SigProfilerMatrixGeneratorR("indels", REF.NAMES[species], OUTPUT$TMP.DIR,
                                          plot=T, exome=F, chrom_based=F, tsb_stat=F, seqInfo=F))
    
    # Copy spectra to output folder
    file.copy(file.path(OUTPUT$TMP.DIR, "output", "plots", "ID_83_plots_indels.pdf"),
              file.path(OUTPUT$DIR.SAMPLE, gsub("${SPECIES}", species, OUTPUT$PDF.SAMPLE, fixed=T)),
              overwrite=T)
    
    
    # Reset SigProfiler folder
    unlink(OUTPUT$TMP.DIR, recursive=T)
    dir.create(OUTPUT$TMP.DIR, showWarnings=F)
    
    # Output species indels to VCF
    write.table(cbind(species.vars[, 1:2], ".", species.vars[, 3:4]),
                file=file.path(OUTPUT$TMP.DIR,
                               paste0("Filtered_indels_in_", toTitleCase(species), ".vcf")),
                sep="\t", quote=F, col.names=F, row.names=F)
    
    # Plot species spectrum
    invisible(SigProfilerMatrixGeneratorR("indels", REF.NAMES[species], OUTPUT$TMP.DIR,
                                          plot=T, exome=F, chrom_based=F, tsb_stat=F, seqInfo=F))
    
    # Copy spectrum to output folder
    file.copy(file.path(OUTPUT$TMP.DIR, "output", "plots", "ID_83_plots_indels.pdf"),
              file.path(OUTPUT$DIR.SPCS, gsub("${SPECIES}", species, OUTPUT$PDF.SPCS, fixed=T)),
              overwrite=T)
    
    # Delete SigProfiler folder
    unlink(OUTPUT$TMP.DIR, recursive=T)
    invisible(gc())
}


cat("\nDone\n\n")
