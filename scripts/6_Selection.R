# Analysis scripts for Cagan, Baez-Ortega et al., 2022
# Step 6: Inference of dN/dS ratios in species with genome annotation

# Adrian Baez-Ortega, 2020-22


# TO RUN THIS SCRIPT IN THE TERMINAL
# ----------------------------------
# Run the commands below in the terminal, replacing '/path/to/CrossSpecies2021' with the 
# path to the CrossSpecies2021 directory.
#
#    cd /path/to/CrossSpecies2021
#    Rscript scripts/6_Selection.R


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
    
    # Path to callable genome file
    CALLABLE.GNM = file.path("data", "original", "CallableGenome.RData"),
    
    # Path to table of gnomAD gene metrics
    GNOMAD.LOF = file.path("data", "original", "gnomAD.v2.1.1.lof_metrics.by_gene.txt"),
    
    # Generic path to RData files of RefCDS databases
    REFCDS = file.path("data", "original", "dNdScv_RefCDS", "RefCDS_${SPECIES}.RData"),
    
    # Generic path to filtered variant (substitution) files
    # (see Methods for details on variant calling and filtering)
    SUBS.FINAL = file.path("data", "original", "FinalCalls_Subs", "${SPECIES}",
                           "FinalCalls_Subs_${SAMPLE}.txt")
)

# Output file paths
OUTPUT = list(
    
    # Path to PDF of dN/dS per species
    PDF = file.path("output", "dNdS_Species.pdf"),
    
    # Path to table of coding mutation counts
    COUNTS = file.path("output", "Coding_Mutation_Counts.txt")
    
)


# Species to consider
SPECIES = c("cat", "colobus", "cow", "dog", "ferret", "horse", "human",
            "mouse", "naked_mole_rat", "rabbit", "rat", "tiger")


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
PACKAGES = c("dndscv", "GenomicRanges")
cat("Loading packages:", paste(PACKAGES, collapse=", "), "\n")
for (package in PACKAGES) {
    suppressWarnings(suppressPackageStartupMessages(library(package, character.only=TRUE,
                                                            quietly=TRUE)))
}

cat("Loading data...\n")
load(INPUT$CALLABLE.GNM)
exclude.list = as.character(as.matrix(read.table(INPUT$EXCLUDE)))
gnomad.lof = read.table(INPUT$GNOMAD.LOF, sep="\t", header=T, as.is=T)
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
cat("Loaded\n")


# Identify haploinsufficient genes
# (LOEUF values in lower 10% in gnomAD; column 'oe_lof_upper')
loeuf.10.genes = gnomad.lof$gene[(gnomad.lof$oe_lof_upper <=
                                      quantile(gnomad.lof$oe_lof_upper, probs=0.1, na.rm=T)) %in% T]


# Run dNdScv for each species
species.dnds.all = species.dnds.hpli = annot.table = coding.bp = coding.call.bp = NULL
for (species in SPECIES) {
    cat("\n\nProcessing species:", species, "\n\n")
    refcds = ifelse(species == "human", "hg19",
                    gsub("${SPECIES}", species, INPUT$REFCDS, fixed=T))
    
    # Build mutation table
    vars.species = NULL
    species.idx = which(sample.info$SPECIES_NAME == species &
                            !(sample.info$SAMPLE_NAME %in% exclude.list))
    for (i in species.idx) {
        
        sample.id = sample.info$SAMPLE_NAME[i]
        cat("Reading variants for sample", sample.id, "\n")
        vars = read.table(gsub("${SAMPLE}", sample.id,
                               gsub("${SPECIES}", species,
                                    INPUT$SUBS.FINAL, fixed=T), fixed=T),
                          sep="\t", header=T, check.names=F, as.is=T, comment.char="",
                          colClasses=c("Chr"="character", "Ref"="character", "Alt"="character"))
        
        cat("Adding", nrow(vars), "variants to table\n")
        vars.species = rbind(vars.species,
                             data.frame(sampleID = sample.id,
                                        normalID = sample.info$NORMAL_NAME[i],
                                        chr = vars$Chr,
                                        pos = vars$Start,
                                        ref = vars$Ref,
                                        alt = vars$Alt,
                                        stringsAsFactors=F))
    }
    
    # Remove variants shared between samples from the same individual
    dup.idx = duplicated(vars.species[, -1])
    vars.species = vars.species[!dup.idx, -2]
    cat("\nDiscarding", sum(dup.idx), "duplicated variants\n")
    cat(nrow(vars.species), "variants after deduplication\n\n")
    
    # Identify genes within the callable genome of every sample
    if (species == "human") {
        data("refcds_hg19", package="dndscv")
    } else {
        load(refcds)
    }
    gene.names = gene.chrs = gene.starts = gene.ends = rep(NA, length(RefCDS))
    for (i in 1:length(RefCDS)) {
        gene = RefCDS[[i]]
        gene.names[i] = gene$gene_name
        gene.chrs[i] = gene$chr
        gene.starts[i] = gene$intervals_cds[1, 1]
        gene.ends[i] = gene$intervals_cds[nrow(gene$intervals_cds), 2]
    }
    stopifnot(all(gene.starts < gene.ends))
    gene.gr = makeGRangesFromDataFrame(data.frame(chr = gene.chrs,
                                                  start = gene.starts, end = gene.ends))
    callable.idx = rep(TRUE, length(gene.gr))
    for (sample.id in sample.info$SAMPLE_NAME[species.idx]) {
        call.gen = callable.genome[[sample.id]]
        callable.idx = callable.idx & overlapsAny(gene.gr, call.gen)
    }
    hpli.callable.idx = callable.idx & (toupper(gene.names) %in% loeuf.10.genes)
    cat(sum(callable.idx), "total genes in callable genome,",
        sum(!callable.idx), "genes excluded\n")
    cat(sum(hpli.callable.idx), "haploinsufficient genes in callable genome\n\n")
    
    
    # Run dNdScv on all genes in callable genome
    cat("Running dNdScv on all genes:\n")
    dnds.out.all = dndscv(vars.species, refdb=refcds, gene_list=gene.names[callable.idx],
                          max_muts_per_gene_per_sample=Inf, max_coding_muts_per_sample=Inf)
    cat("   ", sum(dnds.out.all$sel_cv$qallsubs_cv < 0.05 %in% TRUE),
        "genes with significant dN/dS\n\n")
    
    # Run dNdScv on haploinsufficient genes in callable genome
    cat("Running dNdScv on haploinsufficient genes:\n")
    dnds.out.hpli = dndscv(vars.species, refdb=refcds, gene_list=gene.names[hpli.callable.idx],
                           max_muts_per_gene_per_sample=Inf, max_coding_muts_per_sample=Inf, outp=1)
    
    # Collect coding genome length, callable coding genome length,
    # annotated mutation counts, and global dN/dS, and save output
    coding.bp = c(coding.bp, structure(sum(width(gr_genes)), names=species))
    coding.call.bp = c(coding.call.bp,
                       structure(sum(width(gr_genes)[gr_genes@elementMetadata$names %in%
                                                         gene.names[callable.idx]]), names=species))
    muts = table(dnds.out.all$annotmuts[, c("sampleID", "impact")])
    annot.table = rbind(annot.table,
                        muts[, c("Essential_Splice", "Missense", "Nonsense", "Synonymous"), drop=F])
    species.dnds.all = c(species.dnds.all,
                         structure(list(dnds.out.all$globaldnds), names=species))
    species.dnds.hpli = c(species.dnds.hpli,
                          structure(list(dnds.out.hpli$globaldnds), names=species))
}


# Plot global dN/dS per species, for all genes and haploinsufficient genes
cat("\nPlotting global dN/dS per species...\n")
pdf(OUTPUT$PDF, 15, 6)
par(mar=c(3, 5, 3, 1.75))
cols = c(mis="#377EB8", non="#E41A1C", spl="#4DAF4A", tru="#FF7F00", all="#984EA3", NA, NA)
ylim = c(-4, 4)

dnds.list = list("all"=species.dnds.all, "haploinsufficient"=species.dnds.hpli)

for (i in 1:length(dnds.list)) {
    dnds.mle = sapply(dnds.list[[i]], function(dnds) structure(dnds$mle, names=rownames(dnds)))
    dnds.cilow = sapply(dnds.list[[i]], function(dnds) structure(dnds$cilow, names=rownames(dnds)))
    dnds.cihigh = sapply(dnds.list[[i]], function(dnds) structure(dnds$cihigh, names=rownames(dnds)))
    plot(seq(1, length(dnds.list[[i]]) * (nrow(dnds.mle) + 2)),
         log2(as.numeric(rbind(dnds.mle, NA, NA))),
         main=paste0("Global dN/dS per species (", names(dnds.list)[i], " genes)"),
         col=cols, pch=16, xaxt="n", yaxt="n", ylab="dN/dS", xlab="", cex.lab=1.2, cex=1.1,
         ylim=ylim, xlim=c(2.75, length(dnds.list[[i]]) * (nrow(dnds.mle) + 2) - 3.75),
         panel.first=abline(h=0, col="grey"))
    segments(x0=seq(1, length(dnds.list[[i]]) * (nrow(dnds.mle)+2)),
             y0=log2(as.numeric(rbind(dnds.cilow, NA, NA))),
             y1=log2(as.numeric(rbind(dnds.cihigh, NA, NA))),
             lwd=2.1, col=cols)
    axis(side=2, at=ylim[1]:ylim[2], labels=2^(ylim[1]:ylim[2]), las=1, cex.axis=0.9)
    axis(side=1, labels=names(dnds.list[[i]]), tick=F, cex.axis=1.12,
         at=seq((nrow(dnds.mle) + 1) / 2, by=nrow(dnds.mle)+2, length=length(dnds.list[[i]])))
    abline(v=seq(nrow(dnds.mle) + 1.5, by=nrow(dnds.mle)+2,
                 length=length(dnds.list[[i]])-1), col="grey")
    legend("topleft", legend=rownames(dnds.mle),
           inset=c(0, -0.08), cex=1.1, pch=16, col=cols, horiz=T, bty="n", xpd=T)
}

invisible(dev.off())


# Output table of annotated mutation counts
species.idx = sample.info$SPECIES_NAME[match(rownames(annot.table), sample.info$SAMPLE_NAME)]
annot.table = cbind("Sample"=rownames(annot.table), as.data.frame(annot.table),
                    "Coding_Bp"=coding.bp[species.idx],
                    "Callable_Coding_Bp"=coding.call.bp[species.idx])
write.table(annot.table, file=OUTPUT$COUNTS, sep="\t", quote=F, row.names=F, col.names=T)


cat("\nDone\n\n")
