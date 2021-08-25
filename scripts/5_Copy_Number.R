# Analysis scripts for Cagan, Baez-Ortega et al., 2021
# Step 5: Inference of allele-specific copy number in chromosome-level assemblies

# Adrian Baez-Ortega, 2020-21


# TO RUN THIS SCRIPT IN THE TERMINAL
# ----------------------------------
# Run the commands below in the terminal, replacing '/path/to/CrossSpecies2021' with the 
# path to the CrossSpecies2021 directory.
#
#    cd /path/to/CrossSpecies2021
#    Rscript scripts/5_Copy_Number.R


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
    
    # Path to RData file of bin coverage data
    BIN.COV = file.path("data", "original", "BinCoverage.RData"),
    
    # Path to RData file of heterozygous SNPs per sample
    HET.SNPS = file.path("data", "original", "HetSNPs.RData")
    
)

# Output file paths
OUTPUT = list(
    
    # Path to output directory
    CN.DIR = file.path("output", "Copy_Number"),
    
    # Generic name for PDFs of copy number plots
    CN.PNG = "CN_Plots_${SAMPLE}.png",
    
    # Generic names for copy number tables
    CN.A.TABLE = "CN_Segments_Allele_${SAMPLE}.txt",
    CN.T.TABLE = "CN_Segments_Total_${SAMPLE}.txt"
    
)


# Species to consider: all species with chrom-level assemblies (longer run time)
SPECIES = c("cat", "cow", "dog", "horse", "human", "mouse", "rabbit", "rat")

# Species to consider: species with CNVs only (reduced run time)
# SPECIES = c("cow", "human", "mouse")

# Samples to exclude from CN analysis (low SNP quality)
EXCLUDE.2 = c("PD36813x15", "PD36813x16")

# Read length (bp; for coverage calculation)
READ.LEN = 150

# Minimum length (bins) for initial CN segments
MIN.SEGMENT = 5

# Minimum length (bp) for filtered segments of CN≠2 (1-1)
MIN.LEN = 1e6

# Minimum bin coverage
MIN.COV = 10

# Coverage and VAF thresholds for het SNPs
# (not used: het SNPs provided as input)
# HET.NR = 15
# HET.VAF = c(0.4, 0.6)

# Minimum VAF difference for SNP phasing within a bin
MIN.DIFF = 0.15

# Range of allele CN values to consider
N.VALUES = 0:4

# Penalty parameter for CN changes
PENALTY = 0.3

# Disable scientific notation
options(scipen=999)


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
PACKAGES = c("bbmle", "emdbook", "GenomicRanges", "MASS")
cat("Loading packages:", paste(PACKAGES, collapse=", "), "\n")
for (package in PACKAGES) {
    suppressPackageStartupMessages(library(package, character.only=TRUE, quietly=TRUE))
}

cat("Loading data...\n")
load(INPUT$BIN.COV)
load(INPUT$HET.SNPS)
exclude.list = as.character(as.matrix(read.table(INPUT$EXCLUDE)))
sample.info = read.table(INPUT$SAMPLE.INFO, sep="\t", header=T, as.is=T)
stopifnot(!any(duplicated(sample.info$SAMPLE_NAME)))
cat("Loaded\n")


# Create output directory
dir.create(OUTPUT$CN.DIR, showWarnings=F)

# Initialise list to store temporary objects per sample
sample.data = sapply(sample.info$SAMPLE_NAME, function(x) NULL, simplify=F)



# (1) Copy number calling
# Infer initial CN segments per sample
cn.segments.orig = suppressWarnings(sapply(SPECIES, function(species) {
    
    cat("\n\nProcessing species:", species, "\n")
    species.idx = sample.info$SPECIES_NAME == species &
        !(sample.info$SAMPLE_NAME %in% c(exclude.list, EXCLUDE.2))
    
    sapply(sample.info$SAMPLE_NAME[species.idx], function(sample.id) {
        
        cat("\nProcessing sample:", sample.id, "\n")
        normal.id = sample.info$NORMAL_NAME[match(sample.id, sample.info$SAMPLE_NAME)]
        
        # Calculate adjusted bin coverage for sample and matched normal
        bin.cov = bin.coverage[c(sample.id, normal.id)]
        bin.cov[[1]]$AdjCov = bin.cov[[1]]$Cov / bin.cov[[1]]$Length * READ.LEN
        bin.cov[[2]]$AdjCov = bin.cov[[2]]$Cov / bin.cov[[2]]$Length * READ.LEN
        
        bin.table = bin.cov[[1]][, c("Chrom", "Start", "End")]
        bin.table.gr = with(bin.table, GRanges(Chrom, IRanges(Start, End)))
        
        # Assign het SNPs to coverage bins
        sample.snps = snps.het[[sample.id]]
        sample.snps.gr = with(sample.snps, GRanges(Chrom, IRanges(Pos, Pos)))
        sample.snps$Bin = factor(findOverlaps(sample.snps.gr, bin.table.gr, select="first"),
                                 levels=as.character(1:length(bin.table.gr)))
        sample.snps = sample.snps[!is.na(sample.snps$Bin), ]
        cat("Assigned", nrow(sample.snps), "heterozygous SNPs to genome bins\n")
        
        
        # Infer allele-specific CN for bins with variants, total CN for bins without variants
        # (code adapted from I. MARTINCORENA, Wellcome Sanger Institute [Martincorena et al., 2015])

        # Define two sets of bins:
        #  (i) bins with coverage ≥10 in sample and matched normal, and ≥1 variant (for allele CN)
        # (ii) bins with cov ≥10 in sample and normal, but with no variants (for total CN only)
        bin.idx = 1:nrow(bin.table) %in% sample.snps$Bin
        set1.idx = bin.cov[[1]]$AdjCov >= MIN.COV & bin.cov[[2]]$AdjCov >= MIN.COV & bin.idx
        set2.idx = bin.cov[[1]]$AdjCov >= MIN.COV & bin.cov[[2]]$AdjCov >= MIN.COV & !bin.idx
        cat("Set 1: ", sum(set1.idx), " bins with coverage ≥", MIN.COV, " and ≥1 SNP\n",
            "Set 2: ", sum(set2.idx), " bins with coverage ≥", MIN.COV, " and 0 SNPs\n", sep="")
        
        
        # Beta-binomial model for bin BAF
        cat("Calculating median BAF per bin (with phasing)...\n")
        vaf.obs = sapply(which(set1.idx), function(bin) {
            idx = sample.snps$Bin == bin
            if (sum(idx) > 1) {
                # If the bin has multiple SNPs with very distant VAFs, ensure that
                # all the VAFs are on the same side of 0.5 (to improve CN1/LOH detection)
                vaf.1st = sample.snps$VAF_T[idx][1]
                chg.idx = rep(FALSE, sum(idx))
                if (vaf.1st < 0.5) {
                    chg.idx = sample.snps$VAF_T[idx] > 0.5 &
                        sample.snps$VAF_T[idx] - vaf.1st > MIN.DIFF
                }
                else if (vaf.1st > 0.5) {
                    chg.idx = sample.snps$VAF_T[idx] < 0.5 &
                        vaf.1st - sample.snps$VAF_T[idx] > MIN.DIFF
                }
                sample.snps$VAF_T[idx][chg.idx] = 1 - sample.snps$VAF_T[idx][chg.idx]
            }
            median(sample.snps$VAF_T[idx])  # Observed median bin BAF (Set 1)
        })
        
        ## VAF calculation without SNP phasing (faster)
        # cat("Calculating median BAF per bin (no phasing)...\n")
        # vaf.obs = as.numeric(tapply(sample.snps$VAF_T, sample.snps$Bin, median)[which(set1.idx)])
        
        # Calculate variant allele counts from observed VAFs
        cov.obs = bin.cov[[1]]$AdjCov[set1.idx]            # Observed bin coverage (Set 1)
        cov.obs.s2 = round(bin.cov[[1]]$AdjCov[set2.idx])  # Observed bin coverage (Set 2)
        nv.obs = round(vaf.obs * cov.obs)                  # Median no. of supporting reads (Set 1)
        cov.obs = round(cov.obs)
        
        # Fit beta-binomial to supporting read counts to estimate BAF dispersion
        # ('theta' = dispersion parameter; assuming unbiased het BAF, i.e. p=0.5)
        bb.theta = tryCatch({
            fitdistr(nv.obs, emdbook::dbetabinom, prob=0.5,
                     size=cov.obs, start=list(theta=10), method="L-BFGS-B")$estimate
        }, error = function(e) {
            fitdistr(nv.obs, emdbook::dbetabinom, prob=0.5,
                     size=cov.obs, start=list(theta=10), method="BFGS")$estimate
        })
        
        
        # Negative binomial model for bin coverage
        # The expected coverage of a bin is calculated using the matched normal coverage
        # and the median coverage ratio between sample and matched normal
        cov.exp = bin.cov[[2]]$AdjCov * median(bin.cov[[1]]$AdjCov / bin.cov[[2]]$AdjCov, na.rm=T)
        cov.exp.s2 = cov.exp[set2.idx]  # Expected bin coverage (Set 2)
        cov.exp = cov.exp[set1.idx]     # Expected bin coverage (Set 1)
        
        # Fit negative binomial to observed coverage to estimate coverage dispersion
        # ('size' = dispersion parameter [shape param of gamma mixing distribution])
        nb.size = tryCatch({
            fitdistr(c(cov.obs, cov.obs.s2), dnbinom,
                     mu=c(cov.exp, cov.exp.s2), start=list(size=1), method="L-BFGS-B")$estimate
        }, error = function(e) {
            fitdistr(c(cov.obs, cov.obs.s2), dnbinom,
                     mu=c(cov.exp, cov.exp.s2), start=list(size=1), method="BFGS")$estimate
        })

        
        # Subfunction: Log-likelihood (joint BAF and coverage)
        loglik.copynumber = function(a1, a2, exp.cov, nA, nB, rho, nb.size, bb.theta, ref.BAF) {
            # Define beta-binomial implementation to use
            dbetabinom = emdbook::dbetabinom
            
            # Expected BAF and coverage given the copy number (rho, nA and nB)
            new.baf = (1 - rho + (rho * nA)) / (2 * (1 - rho) + rho * (nA + nB))
            new.baf = new.baf * ref.BAF / (new.baf * ref.BAF + (1 - new.baf) * (1 - ref.BAF))
            new.mu = exp.cov + exp.cov * (nA + nB - 2) * rho / 2
            
            # Likelihoods
            # If a1=0 (no BAF info), return only neg binom likelihood (for total CN);
            # otherwise combine neg binom and beta-binom likelihoods (for allele CN)
            if (a1 == 0) { 
                LL = sum(dnbinom(x=a1+a2, mu=new.mu, size=nb.size, log=T))
            } else {
                ll.baf = sum(dbetabinom(x=a1, size=a1+a2, prob=new.baf, theta=bb.theta, log=T))
                ll.cov = sum(dnbinom(x=a1+a2, mu=new.mu, size=nb.size, log=T))
                LL = ll.baf + ll.cov
            }
            return(-LL)
        }
        
        
        # Subfunction: Optimise nA, nB and rho by exhaustive search of
        # nA and nB and numerical optimisation of rho
        optimise.copynumber = function(n.values, a1, a2, exp.cov,
                                       nb.size, bb.theta, penalty.matrix, ref.BAF) {
            LL.mat = array(NA, dim=c(length(n.values), length(n.values)))
            rho.mat = array(NA, dim=c(length(n.values), length(n.values)))
            rho.CIlow.mat = array(NA, dim=c(length(n.values), length(n.values)))
            rho.CIhigh.mat = array(NA, dim=c(length(n.values), length(n.values)))
            for (j in 1:length(n.values)) {
                for (h in 1:length(n.values)) {
                    nA = n.values[j]
                    nB = n.values[h]
                    # Original pars: start=list(rho=0.5), lower=c(rho=0.0001)
                    m = mle2(loglik.copynumber, start=list(rho=0.925),
                             data=list(a1=a1, a2=a2, exp.cov=exp.cov, nA=nA, nB=nB,
                                       nb.size=nb.size, bb.theta=bb.theta, ref.BAF=ref.BAF),
                             method="Brent", lower=c(rho=0.85), upper=c(rho=0.9999))
                    rho.mat[j, h] = summary(m)@coef[1]
                    LL.mat[j, h] = logLik(m)[1]
                }
            }
            
            # Penalise complexity with a prior (obtain posterior probabilities)
            PP.mat = LL.mat + penalty.matrix
            likely.cn = which(exp(PP.mat - max(PP.mat)) >= 0.10, arr.ind=T)
            
            # Report all the likely solutions in order of their posterior prob
            likely.cndf = data.frame(nA=n.values[likely.cn[, 1]], nB=n.values[likely.cn[, 2]],
                                     rhoMLE=rho.mat[likely.cn], LL=LL.mat[likely.cn],
                                     PP=PP.mat[likely.cn], rho2.5=NA, rho97.5=NA)
            likely.cndf = likely.cndf[order(likely.cndf$PP, decreasing=T), ]
            likely.cndf$rel.prob = exp(likely.cndf$PP - max(likely.cndf$PP))
            return(likely.cndf)
        }
        
        
        # Subfunction: Generate penalty matrix against complex solutions,
        # which applies a fixed penalty ('param') for every copy number step
        # away from 1-1; values closer to 1/0 imply weaker/stronger penalties
        get.penalty.matrix = function(n.values, param) {
            n = length(n.values)
            penalty.matrix = array(NA, dim=c(n, n))
            for (j in 1:n) {
                for (h in 1:n) {
                    penalty.matrix[j,h] = log(param ^ sum(abs(c(n.values[j], n.values[h]) - 1)))
                }
            }
            return(penalty.matrix)
        }
        
        
        # Search for the most likely copy number for every bin in the current sample
        cat("Estimating copy number per bin...\n")
        p.vec = rep(NA, sum(set1.idx))
        penalty.matrix = get.penalty.matrix(N.VALUES, PENALTY)
        cn.table = data.frame(CN_total = rep(NA, nrow(bin.table)),
                              CN_allele = rep(NA, nrow(bin.table)),
                              nMajor = rep(NA, nrow(bin.table)),
                              nMinor = rep(NA, nrow(bin.table)),
                              nA = rep(NA, nrow(bin.table)),
                              nB = rep(NA, nrow(bin.table)),
                              rho_MLE = rep(NA, nrow(bin.table)),
                              LL = rep(NA, nrow(bin.table)),
                              PP = rep(NA, nrow(bin.table)),
                              rho_2.5 = rep(NA, nrow(bin.table)),
                              rho_97.5 = rep(NA, nrow(bin.table)),
                              rel_prob = rep(NA, nrow(bin.table)))
        
        # Obtain most likely CN state for bins in Sets 1 and 2
        for (i in 1:sum(set1.idx)) {
            cn = optimise.copynumber(N.VALUES,
                                     a1=nv.obs[i], a2=cov.obs[i]-nv.obs[i], exp.cov=cov.exp[i],
                                     nb.size, bb.theta, penalty.matrix, ref.BAF=0.5)
            n.major = max(cn[1, 1:2])
            n.minor = min(cn[1, 1:2])
            cn.table[which(set1.idx)[i], ] = c(n.major + n.minor,
                                               paste0(n.major, "-", n.minor),
                                               n.major, n.minor, cn[1, ])
        }
        for (i in 1:sum(set2.idx)) {
            # For a1=0, optimise.copynumber estimates total CN only
            cn = optimise.copynumber(N.VALUES,
                                     a1=0, a2=cov.obs.s2[i], exp.cov=cov.exp.s2[i],
                                     nb.size, bb.theta, penalty.matrix, ref.BAF=0.5)
            cn.table[which(set2.idx)[i], ] = c(sum(cn[1, 1:2]), NA, NA, NA, cn[1, ])
        }
        
        
        # Define CN segments above minimum length
        # (segmentation is done separately for total CN and allele-specific CN)
        cat("Defining copy number segments...\n")
        cn.segs.prelim = list("total"=NULL, "allele"=NULL)
        for (name in c("total", "allele")) {
            start.idx = 1
            while (start.idx < nrow(cn.table)) {
                chrom = bin.table$Chrom[start.idx]
                chrom.end = rev(which(bin.table$Chrom == chrom))[1]
                colname = paste0("CN_", name)
                cn = cn.table[start.idx, colname]
                if (is.na(cn)) {
                    start.idx = start.idx + 1
                }
                else {
                    end.idx = start.idx - 1 + 
                        which((cn.table[-(1:start.idx), colname] != cn) |
                                  is.na(cn.table[-(1:start.idx), colname]))
                    if (length(end.idx) == 0 | end.idx[1] > chrom.end) {
                        end.idx = chrom.end
                    }
                    else {
                        end.idx = end.idx[1]
                    }
                    if (end.idx - start.idx + 1 >= MIN.SEGMENT) {
                        if (name == "total") {
                            n.major = n.minor = NA
                        }
                        else {
                            n.major = cn.table$nMajor[start.idx]
                            n.minor = cn.table$nMinor[start.idx]
                        }
                        cn.segs.prelim[[name]] = rbind(cn.segs.prelim[[name]],
                                                       data.frame(StartBin=start.idx,
                                                                  EndBin=end.idx,
                                                                  Chrom=chrom,
                                                                  Start=bin.table$Start[start.idx],
                                                                  End=bin.table$End[end.idx],
                                                                  Length=bin.table$End[end.idx] -
                                                                      bin.table$Start[start.idx] + 1,
                                                                  CN=cn,
                                                                  nMajor=n.major,
                                                                  nMinor=n.minor,
                                                                  stringsAsFactors=F))
                    }
                    start.idx = end.idx + 1
                }
            }
        }
        
        # Merge adjacent segments with same CN and separated by a gap below the min segment length
        cn.segments = list("total"=NULL, "allele"=NULL)
        for (name in c("total", "allele")) {
            i = 1
            while (i <= nrow(cn.segs.prelim[[name]])) {
                segment = cn.segs.prelim[[name]][i, ]
                j = i + 1
                while (j <= nrow(cn.segs.prelim[[name]]) &
                       cn.segs.prelim[[name]]$Chrom[j] == segment$Chrom &
                       cn.segs.prelim[[name]]$CN[j] == segment$CN &
                       cn.segs.prelim[[name]]$StartBin[j] - segment$EndBin <= MIN.SEGMENT) {
                    segment$EndBin = cn.segs.prelim[[name]]$EndBin[j]
                    segment$End = cn.segs.prelim[[name]]$End[j]
                    segment$Length = segment$End - segment$Start + 1
                    j = j + 1
                }
                cn.segments[[name]] = rbind(cn.segments[[name]], segment)
                i = j
            }
            rownames(cn.segments[[name]]) = NULL
        }
        
        # Store temporary objects and return CN segments
        sample.data[[sample.id]] <<- list("bin.table"=bin.table,
                                          "set1.idx"=set1.idx, "set2.idx"=set2.idx,
                                          "cov.obs"=cov.obs, "cov.exp"=cov.exp, "vaf.obs"=vaf.obs,
                                          "cov.obs.s2"=cov.obs.s2, "cov.exp.s2"=cov.exp.s2)
        rm(cn.segs.prelim, sample.snps, sample.snps.gr, bin.cov, bin.table,
           bin.table.gr, cov.obs, cov.exp, vaf.obs, cov.obs.s2, cov.exp.s2)
        invisible(gc())
        cn.segments
    
    }, simplify=F)
}, simplify=F))



# (2) Segment filtering
# For each species, filter copy number segments different from the default
# CN=2 (1-1) which are shared among multiple samples (as shared CNVs are not expected)

# Process CN segments for each species
cat("\n\nFiltering CN segments...\n")
cn.segments.filt = sapply(SPECIES, function(x) NULL, simplify=F)
for (species in SPECIES) {
    
    # Collect all segments without standard (1-1) CN
    cn.levels.total = unique(unlist(sapply(cn.segments.orig[[species]], function(cn) {
        cn$total$CN[cn$total$CN != 2]
    })))
    cn.levels.allele = unique(unlist(sapply(cn.segments.orig[[species]], function(cn) {
        cn$allele$CN[cn$allele$CN != "1-1"]
    })))
    
    cnvs = data.frame(Chrom="1", Start=1, End=1, stringsAsFactors=F)
    for (cn in cn.segments.orig[[species]]) {
        for (lvl in as.character(cn.levels.total)) {
            cnvs = rbind(cnvs,
                         cn$total[as.character(cn$total$CN) == lvl, c("Chrom","Start","End")])
        }
        for (lvl in cn.levels.allele) {
            cnvs = rbind(cnvs,
                         cn$allele[cn$allele$CN == lvl, c("Chrom","Start","End")])
        }
    }
    cnvs.gr = makeGRangesFromDataFrame(cnvs)
    
    # Segment filtering:
    # For each sample, discard segments of CN≠2 (1-1) and presenting
    # >2 overlaps with any segments of CN≠2 (1-1) (including itself).
    # Also discard segments of CN≠2 (1-1) with length < 1 Mb.
    cn.segments.filt[[species]] = cn.segments.orig[[species]]
    for (id in names(cn.segments.filt[[species]])) {
        cn.sample = cn.segments.filt[[species]][[id]]$total
        keep.total.idx = sapply(1:nrow(cn.sample), function(i) {
            if (cn.sample$CN[i] == 2) {
                TRUE
            } else {
                cn.sample$Length[i] >= MIN.LEN &
                    countOverlaps(makeGRangesFromDataFrame(cn.sample[i, ]), cnvs.gr) < 3
            }
        })
        cn.sample = cn.segments.filt[[species]][[id]]$allele
        keep.allele.idx = sapply(1:nrow(cn.sample), function(i) {
            if (cn.sample$CN[i] == "1-1") {
                TRUE
            } else {
                cn.sample$Length[i] >= MIN.LEN &
                    countOverlaps(makeGRangesFromDataFrame(cn.sample[i, ]), cnvs.gr) < 3
            }
        })
        cn.segments.filt[[species]][[id]]$total = 
            cn.segments.filt[[species]][[id]]$total[keep.total.idx, ]
        cn.segments.filt[[species]][[id]]$allele =
            cn.segments.filt[[species]][[id]]$allele[keep.allele.idx, ]
    }
}



# (3) Plotting and outputting
cat("Plotting and outputting CN segments...\n")
for (i in 1:length(cn.segments.filt)) {
    for (j in 1:length(cn.segments.filt[[i]])) {
        sample.id = names(cn.segments.filt[[i]])[j]
        normal.id = sample.info$NORMAL_NAME[match(sample.id, sample.info$SAMPLE_NAME)]
        bin.table = sample.data[[sample.id]]$bin.table
        set1.idx = sample.data[[sample.id]]$set1.idx
        set2.idx = sample.data[[sample.id]]$set2.idx
        cov.obs = sample.data[[sample.id]]$cov.obs
        cov.exp = sample.data[[sample.id]]$cov.exp
        vaf.obs = sample.data[[sample.id]]$vaf.obs
        cov.obs.s2 = sample.data[[sample.id]]$cov.obs.s2
        cov.exp.s2 = sample.data[[sample.id]]$cov.exp.s2
        
        chroms = unique(bin.table$Chrom)
        cutoffs = c(match(chroms, bin.table$Chrom), nrow(bin.table))
        cols = rep("black", length(chroms))
        cols[seq_along(cols) %% 2 == 0] = "darkgrey"
        cols.s2 = rep(cols, table(bin.table$Chrom)[chroms])[set2.idx]
        cols = rep(cols, table(bin.table$Chrom)[chroms])[set1.idx]
        xlm = c(1, nrow(bin.table))
        offset = 0.1; cxlab = 1.5; cxmain = 1.7
        
        png(file.path(OUTPUT$CN.DIR, gsub("${SAMPLE}", sample.id, OUTPUT$CN.PNG, fixed=T)),
            17, 8, "in", res=120)
        par(mfrow=c(3, 1), mar=c(0, 5, 3, 1.25), oma=c(2.5, 0, 3, 0))
        
        # Plot obs/exp coverage ratio
        plot(which(set1.idx), cov.obs / cov.exp,
             col=cols, pch=16, cex=0.6, las=2, xlab="", ylab="Obs/Exp coverage ratio", 
             xpd=NA, xlim=xlm, xaxs="i", xaxt="n", panel.first=abline(v=cutoffs, col="grey"),
             main="Observed/Expected bin coverage ratio", cex.lab=cxlab, cex.main=cxmain)
        points(which(set2.idx), cov.obs.s2 / cov.exp.s2, col=cols.s2, pch=16, cex=0.6)
        title(paste0(sample.id, " (normal ", normal.id, "; ",
                     prettyNum(sum(set1.idx), big.mark=","), " bins with SNPs)"),
              cex.main=2, line=3.5, xpd=NA)
        
        # Plot het BAF
        plot(which(set1.idx), vaf.obs,
             col=cols, pch=16, cex=0.6, las=2, xlab="", ylab="BAF",
             xlim=xlm, ylim=c(0, 1), xaxs="i", xaxt="n", panel.first=abline(v=cutoffs, col="grey"),
             main="Heterozygous SNP BAF", cex.lab=cxlab, cex.main=cxmain)
        
        # Plot CN segments
        plot(1, type="n", las=2, xlab="", ylab="Copy number",
             xlim=xlm, ylim=c(0, 4), xaxs="i", xaxt="n", panel.first=abline(v=cutoffs, col="grey"),
             main="Copy number segments (filtered)", cex.lab=cxlab, cex.main=cxmain)
        segments(x0=cn.segments.filt[[i]][[j]]$total$StartBin,
                 x1=cn.segments.filt[[i]][[j]]$total$EndBin,
                 y0=cn.segments.filt[[i]][[j]]$total$CN, col="forestgreen", lwd=4, lend=2)
        segments(x0=cn.segments.filt[[i]][[j]]$allele$StartBin,
                 x1=cn.segments.filt[[i]][[j]]$allele$EndBin,
                 y0=cn.segments.filt[[i]][[j]]$allele$nMajor + offset, col="red", lwd=4, lend=2)
        segments(x0=cn.segments.filt[[i]][[j]]$allele$StartBin,
                 x1=cn.segments.filt[[i]][[j]]$allele$EndBin,
                 y0=cn.segments.filt[[i]][[j]]$allele$nMinor - offset, col="blue", lwd=4, lend=2)
        
        mtext(chroms, side=1, at=(cutoffs[-length(cutoffs)] + cutoffs[-1]) / 2, line=0.7, font=2)
        invisible(dev.off())
    }
}


# Output filtered CN segments
for (i in 1:length(cn.segments.filt)) {
    for (j in 1:length(cn.segments.filt[[i]])) {
        sample.id = names(cn.segments.filt[[i]])[j]
        write.table(cn.segments.filt[[i]][[j]]$total[, 1:7], sep="\t", quote=F, row.names=F,
                    file=file.path(OUTPUT$CN.DIR,
                                   gsub("${SAMPLE}", sample.id, OUTPUT$CN.T.TABLE, fixed=T)))
        write.table(cn.segments.filt[[i]][[j]]$allele, sep="\t", quote=F, row.names=F,
                    file=file.path(OUTPUT$CN.DIR,
                                   gsub("${SAMPLE}", sample.id, OUTPUT$CN.A.TABLE, fixed=T)))
    }
}


cat("\nDone\n\n")
