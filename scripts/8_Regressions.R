# Analysis scripts for Cagan, Baez-Ortega et al., 2021
# Step 8: Regression analyses of somatic mutation burdens and rates

# Adrian Baez-Ortega, 2020-21


# TO RUN THIS SCRIPT IN THE TERMINAL
# ----------------------------------
# Run the commands below in the terminal, replacing '/path/to/CrossSpecies2021' with the 
# path to the CrossSpecies2021 directory.
#
#    cd /path/to/CrossSpecies2021
#    Rscript scripts/8_Regressions.R


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
    
    # Path to RData file of mutation burdens and rates
    BURDEN = file.path("data", "processed", "Burdens_Rates.RData"),
    
    # Path to list of samples to exclude
    EXCLUDE = file.path("data", "processed", "SamplesToExclude.txt"),
    
    # Path to table of collected species lifespan
    LIFESPAN = file.path("data", "original", "Lifespan_Collected.txt"),
    
    # Path to species tree from TimeTree
    TREE = file.path("data", "original", "SpeciesTree.newick.tree"),
    
    # Path to code for Bayesian hierarchical normal regression model
    BHN.MODEL = file.path("scripts", "BHN_Regression.stan")

)

# Output file paths
OUTPUT = list(
    
    # Path to PDF of mutation rate regressions (LME/BHN model)
    RATES.PDF = file.path("output", "Regression_Rates_LME-BHN.pdf"),
    
    # Path to PDF of allometric regressions among rate, mass and lifespan
    ALLOM.PDF = file.path("output", "Regression_Allometric.pdf"),

    # Path to PDF of mutation burden vs age regressions (simple LM)
    BURDEN.PDF = file.path("output", "Regression_Burden-Age_LM.pdf"),
    
    # Path to PDF of regression model comparison results
    COMP.PDF = file.path("output", "Regression_Model_Comparison.pdf"),
    
    # Path to PDF of bootstrapped rate vs lifespan regression
    LIFSP.PDF = file.path("output", "Regression_Lifespan_Bootstrap.pdf")

)


# Samples to exclude from mutation rate regressions: samples sharing
# the majority of variants with another sample from the same individual
EXCLUDE.2 = c("CATD0002b_lo0003", "MD6267ab_lo0003")

# Abbreviated species labels (for plots)
SP.LAB = c("cat"="CAT", "colobus"="BWC", "cow"="COW", "dog"="DOG", "ferret"="FER",
           "giraffe"="GIR", "horse"="HOR", "human"="HUM", "lion"="LION",
           "mouse"="MOU", "naked_mole_rat"="NMR", "rabbit"="RAB",
           "rat"="RAT", "ring_tailed_lemur"="RTL", "tiger"="TIG")

# MCMC/bootstrap parameters
ITER = 5000                      # Total iterations 
SEED = 0xC0FFEE                  # Random seed
CORES = parallel::detectCores()  # Number of CPUs (detected automatically)

# Function: multiple pattern replacement (gsub)
multi.gsub = function(patterns, replacements, x, ...) {
    stopifnot(length(patterns) == length(replacements))
    for (i in 1:length(patterns)) {
        x = gsub(patterns[i], replacements[i], x, ...)
    }
    x
}

# Function: Regression R-squared from `lm` model using the standard
# formula [MSS/(MSS+RSS)]. This function ignores model weights, and
# does not use a different formula for models without an intercept.
r2.lm = function(model) {
    r = model$residuals         # residuals
    f = model$fitted.values     # fitted values
    mss = sum((f - mean(f))^2)  # formula for free-intercept models
    rss = sum(r^2)
    mss / (mss + rss)
}

# Function: Regression R-squared from `lme` model
#  obj: lme object with regression model output
#    X: matrix of predictor values for prediction of fitted values
#    y: vector of observed values of the response
r2.lme = function(obj, X, y) {
    library(nlme)
    beta = fixed.effects(obj)
    stopifnot(identical(dim(X), c(length(y), length(beta))))
    f = as.numeric(X %*% beta)
    idx = !is.na(y) & !is.na(f)
    f = f[idx]                  # fitted values
    r = y[idx] - f              # residuals
    mss = sum((f - mean(f))^2)  # formula for free-intercept models
    rss = sum(r^2)
    mss / (mss + rss)
}

# Function: Regression R-squared from `stanfit` object
#  obj: stanfit object with regression model MCMC samples
#    X: matrix of predictor values for prediction of fitted values
#    y: vector of observed values of the response
r2.stan = function(obj, X, y) {
    library(rstan)
    beta = apply(extract(obj, "beta")[[1]], 2, median)
    stopifnot(identical(dim(X), c(length(y), length(beta))))
    f = as.numeric(X %*% beta)
    idx = !is.na(y) & !is.na(f)
    f = f[idx]                  # fitted values
    r = y[idx] - f              # residuals
    mss = sum((f - mean(f))^2)  # formula for free-intercept models
    rss = sum(r^2)
    mss / (mss + rss)
}

# Function: Regression 95% credible intervals from `stanfit` object
# (only for models with a single explanatory variable)
#  obj: stanfit object with regression model MCMC samples
#    X: matrix of predictor values for credible interval calculation
ci.stan = function(obj, X) {
    library(rstan)
    beta = extract(obj, "beta")[[1]]
    stopifnot(ncol(X) == ncol(beta) & ncol(X) <= 2)
    y = apply(beta, 1, function(b) X %*% b)
    t(apply(y, 1, quantile, probs=c(0.025, 0.5, 0.975)))
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
PACKAGES = c("caper", "nlme", "RColorBrewer", "rstan", "scales", "parallel", "tools")
cat("Loading packages:", paste(PACKAGES, collapse=", "), "\n")
for (package in PACKAGES) {
    suppressPackageStartupMessages(suppressWarnings(library(package, character.only=TRUE)))
}

cat("Loading data...\n")
load(INPUT$BURDEN)
species.tree = read.tree(INPUT$TREE)
exclude.list = as.character(as.matrix(read.table(INPUT$EXCLUDE)))
lifespans = read.table(INPUT$LIFESPAN, sep="\t", header=T, as.is=T)
stopifnot(!any(sample.data$Sample %in% exclude.list))
cat("Loaded\n\n")


# (1) Regress total burden and signature-specific burden on age
#     Model: simple linear model (LM) using mean burdens per individual

cairo_pdf(OUTPUT$BURDEN.PDF, 10.5, 10.5, onefile=T)
par(mfrow=c(5, 3), mar=c(2.3, 2.5, 0, 0), oma=c(1.2, 1, 1, 1), mgp=c(1.75, 0.45, 0), tck=-0.03)
cols = c("dodgerblue3", "forestgreen", "darkorange", "darkorchid")

# Discard data for individuals/species with unknown ages/lifespans
sample.data = sample.data[!is.na(sample.data$Age_years), ]
indiv.data = indiv.data[!is.na(indiv.data$Age_years), ]
species.data = species.data[!is.na(species.data$Lifespan80), ]

# (A) Plot total burden vs age for each species
cat("Using simple LM to regress mutation burdens on age...\n")
for (species in unique(sample.data$Species)) {
    
    # Retrieve ages and burdens per sample and per individual
    sample.age = sample.data$Age_years[sample.data$Species == species]
    sample.burden = sample.data$Mut_burden[sample.data$Species == species]
    indiv.age = indiv.data$Age_years[indiv.data$Species == species]
    indiv.burden = indiv.data$Mut_burden[indiv.data$Species == species]
    
    # Initialise plot
    plot(sample.age, sample.burden, type="n", cex.lab=1.2, xlab="", ylab="",
         xlim=c(0, max(sample.age, na.rm=T) * 1.05), ylim=c(0, max(sample.burden) * 1.09))
    if (match(species, unique(sample.data$Species)) %% 3 == 1) {
        mtext("Mutations", side=2, line=1.8, cex=0.8)
    }
    if (match(species, rev(unique(sample.data$Species))) %in% 1:3) {
        mtext("Age (years)", side=1, line=1.7, cex=0.8)
    }
    
    # If >1 individual with known age: regress mean burden on age
    if (sum(!is.na(indiv.age)) > 1) {
        burden.lm = lm(indiv.burden ~ indiv.age)
        ci.x = seq(-10, max(sample.age, na.rm=T) * 1.2, len=100)
        ci.y = suppressWarnings(predict(burden.lm, newdata=data.frame("indiv.age"=ci.x),
                                        interval="confidence", level=0.95))
        polygon(c(ci.x, rev(ci.x)), c(ci.y[,2], rev(ci.y[,3])), col=alpha(cols[1], 0.25), border=NA)
        abline(burden.lm, lwd=2.5, col=alpha(cols[1], 0.7))
    }
    
    # Plot points per sample, coloured by individual
    points(sample.age, sample.burden, pch=21, col=cols[1], cex=1.2, lwd=2,
           bg=brewer.pal(12, "Set3")[match(sample.data$Individual[sample.data$Species == species],
                                           indiv.data$Individual[indiv.data$Species == species])])
    title(multi.gsub(c("Ring_tailed", "mole_rat", "_"), c("Ring-tailed", "mole-rat", " "),
                     toTitleCase(species)), cex.main=1.5, line=-1.7)
}

# (B) Plot signature burdens vs age for each species
for (species in unique(sample.data$Species)) {
    
    # Retrieve ages and burdens per sample and per individual
    sig.idx = c("SBS1_burden", "SBSB_burden", "SBSC_burden")
    sample.age = sample.data$Age_years[sample.data$Species == species]
    sample.burden = sample.data[sample.data$Species == species, sig.idx]
    indiv.age = indiv.data$Age_years[indiv.data$Species == species]
    indiv.burden = indiv.data[indiv.data$Species == species, sig.idx]
    
    # Initialise plot
    plot(sample.age, sample.burden[, 1], type="n", cex.lab=1.2, xlab="", ylab="",
         xlim=c(0, max(sample.age, na.rm=T) * 1.05), ylim=c(0, max(sample.burden) * 1.09))
    if (match(species, unique(sample.data$Species)) %% 3 == 1) {
        mtext("Mutations", side=2, line=1.8, cex=0.8)
    }
    if (match(species, rev(unique(sample.data$Species))) %in% 1:3) {
        mtext("Age (years)", side=1, line=1.7, cex=0.8)
    }
    
    for (j in 1:3) {
        # If >1 individual with known age: regress mean burden on age
        if (sum(!is.na(indiv.age)) > 1) {
            burden.lm = lm(indiv.burden[, j] ~ indiv.age)
            ci.x = seq(-10, max(sample.age, na.rm=T) * 1.2, len=100)
            ci.y = suppressWarnings(predict(burden.lm, newdata=data.frame("indiv.age"=ci.x),
                                            interval="confidence", level=0.95))
            polygon(c(ci.x, rev(ci.x)), c(ci.y[,2], rev(ci.y[,3])),
                    col=alpha(cols[j+1], 0.25), border=NA)
            abline(burden.lm, lwd=2.5, col=alpha(cols[j+1], 0.7))
        }
        # Plot points per sample, coloured by signature
        points(sample.age, sample.burden[, j], pch=16, cex=1.5, col=alpha(cols[j+1], 0.9))
    }
    title(multi.gsub(c("Ring_tailed", "mole_rat", "_"), c("Ring-tailed", "mole-rat", " "),
                     toTitleCase(species)), cex.main=1.5, line=-1.7)
}

invisible(dev.off())


# (2) Regress various mutation rates/burdens on inverse lifespan (1/Lifespan)
#     Models: linear mixed-effects (LME) model, Bayesian hier. normal (BHN) model
#     (constrained intercept)

# Discard samples with a majority of variants shared
sample.data = sample.data[!(sample.data$Sample %in% EXCLUDE.2), ]

# (A) LME regression: use constrained (zero-intercept) LME model to
#     regress mutation rates on 1/Lifespan, while accounting for the
#     hierarchical data structure and heteroscedasticity among species
cat("Using LME model to regress mutation rates on 1/Lifespan...\n")
col.names = paste0(c("Mut", "Indel", "SBS1", "SBSB", "SBSC"), "_rate")
lme.rates.lifespan80 = sapply(col.names, function(name) {
    lme(fixed = as.formula(paste(name, "~ Inverse_Lifespan80 - 1")),
        random = ~ Inverse_Lifespan80 - 1 | Species/Individual,
        weights = varIdent(form = ~ 1 | Species),
        data = sample.data,
        method = "ML")
}, simplify=F)

# Use simple LM to regress mean mtDNA rate per species (high data sparsity)
lme.rates.lifespan80$mtDNA_rate = lm(mtDNA_rate ~ Inverse_Lifespan80 - 1,
                                     data = species.data)

# (B) Bayesian regression: use constrained hierarchical normal model to
#     regress mean burden per individual on 1/Lifespan, using age as an offset
#     (for model definition, see file scripts/BHN_Regression.stan)
cat("Using BHN model to regress mutation burdens on 1/Lifespan...\n")
bhn.model = stan_model(INPUT$BHN.MODEL)
col.names = paste0(c("Mut", "Indel", "SBS1", "SBSB", "SBSC", "mtDNA"), "_burden")
bhn.rates.lifespan80 = suppressWarnings(sapply(col.names, function(name) {
    sampling(bhn.model,
             iter = ITER, seed = SEED, cores = CORES, refresh = 0,
             data = list(K = 1,
                         N = nrow(indiv.data),
                         G = length(unique(indiv.data$Species)),
                         y = indiv.data[, name],
                         X = as.matrix(indiv.data$Inverse_Lifespan80),
                         t = indiv.data$Age_years,
                         group = as.integer(as.factor(indiv.data$Species))),
             init = function() {
                 # Define random initial values for beta_g within model constraints
                 list(beta_g = matrix(runif(length(unique(indiv.data$Species)), 1, 100), ncol=1))
             })
}, simplify=F))


# (3) Perform allometric regressions and linear regressions of allometric residuals,
#     for pairwise combinations of mutation rate, body mass and lifespan
#     The allometric equation is:  Y = a * X^b
#     This curve is perfectly equivalent to an exponential curve on the logarithm of X:
#           Y = a * X^b = a * exp(b*log(x))
#     Thus  log(Y) = a' + b * log(x),  which corresponds to a linear model on a log-log scale

# (A) Allometric regressions of: rate vs mass, rate vs lifespan,
#     lifespan vs mass, mass vs lifespan, metabolic rate vs mass
rownames(species.data) = species.data$Species
species.data$Log10_Mut_rate = log10(species.data$Mut_rate)
species.data$Log10_Lifespan80 = log10(species.data$Lifespan80)
lm.allom.rate.mass = lm(Log10_Mut_rate ~ Log10_Mass, species.data)
lm.allom.rate.lifespan = lm(Log10_Mut_rate ~ Log10_Lifespan80, species.data)
lm.allom.lifespan.mass = lm(Log10_Lifespan80 ~ Log10_Mass, species.data)
lm.allom.mass.lifespan = lm(Log10_Mass ~ Log10_Lifespan80, species.data)
lm.allom.bmr.mass = lm(Log10_Metabolic_rate ~ Log10_Mass, species.data)

# (B) Allometric regression residuals
stopifnot(identical(species.data$Species, names(residuals(lm.allom.rate.mass))) &
              identical(species.data$Species, names(residuals(lm.allom.rate.lifespan))) &
              identical(species.data$Species, names(residuals(lm.allom.lifespan.mass))) &
              identical(species.data$Species, names(residuals(lm.allom.mass.lifespan))))
species.data$Residuals_Rate_Mass = residuals(lm.allom.rate.mass)
species.data$Residuals_Rate_Lifespan = residuals(lm.allom.rate.lifespan)
species.data$Residuals_Lifespan_Mass = residuals(lm.allom.lifespan.mass)
species.data$Residuals_Mass_Lifespan = residuals(lm.allom.mass.lifespan)
species.data$Residuals_BMR_Mass =
    residuals(lm.allom.bmr.mass)[match(species.data$Species, names(residuals(lm.allom.bmr.mass)))]
sample.data$Residuals_BMR_Mass =
    species.data$Residuals_BMR_Mass[match(sample.data$Species, species.data$Species)]

# (C) Residual regressions for: rate vs lifespan (controlled for mass),
#                               rate vs mass (controlled for lifespan)
lm.resid.rate.lifespan = lm(Residuals_Rate_Mass ~ Residuals_Lifespan_Mass, species.data)
lm.resid.rate.mass = lm(Residuals_Rate_Lifespan ~ Residuals_Mass_Lifespan, species.data)


# (4) Regress substitution rate/burden on life-history variables,
#     both individually and in combination with 1/Lifespan
#     Models: LME model, BHN model (free intercept)

# (A) LME regression: use free-intercept LME model to regress mutation 
#     rate on life-history traits, while accounting for the hierarchical
#     structure of the data and heteroscedasticity among species
cat("Using LME model to regress mutation rate on life-history traits...\n")
col.names = c("Inverse_Lifespan80", "Log10_Mass", "Litter_size",
              "Log10_Metabolic_rate", "Residuals_BMR_Mass")

# Regression on individual variables
lme.rate.anage.single = sapply(col.names, function(name) {
    lme(fixed = as.formula(paste("Mut_rate ~", name)),
        random = as.formula(paste("~", name, "- 1 | Species/Individual")),
        weights = varIdent(form = ~ 1 | Species),
        data = sample.data[!is.na(sample.data[, name]), ],
        method = "ML")
}, simplify=F)

# Regression on variables paired with 1/Lifespan
lme.rate.anage.paired = sapply(col.names[-1], function(name) {
    lme(fixed = as.formula(paste("Mut_rate ~ Inverse_Lifespan80 +", name)),
        random = ~ Inverse_Lifespan80 - 1 | Species/Individual,
        weights = varIdent(form = ~ 1 | Species),
        data = sample.data[!is.na(sample.data[, name]), ],
        method = "ML")
}, simplify=F)

# Obtain nonparametric 95% bootstrap intervals (resampling individuals)
cat("Calculating nonparametric 95% bootstrap intervals for LME regressions...\n")
lme.rate.anage.single.bs = sapply(col.names, function(name) {
    # Make parallel cluster and export variables and packages
    cl = makeCluster(CORES)
    clusterExport(cl, c("sample.data", "name"))
    clusterEvalQ(cl, library(nlme))
    set.seed(SEED)
    # Run parallel sapply
    out = t(parSapply(cl, 1:ITER, function(i) {
        # In each bootstrap replicate, resample individuals with replacement
        indivs = sample(unique(sample.data$Individual),
                        length(unique(sample.data$Individual)), replace=T)
        idx = sample.data$Individual %in% indivs & !is.na(sample.data[, name])
        tryCatch({
            # Return fixed effects of fitted LME model
            fixed.effects(lme(fixed = as.formula(paste("Mut_rate ~", name)),
                              random = as.formula(paste("~", name, "- 1 | Species/Individual")),
                              weights = varIdent(form = ~ 1 | Species),
                              data = sample.data[idx, ],
                              method="ML"))
        }, error = function(e) {
            rep(NA, 2)
        })
    }))
    stopCluster(cl)
    out
}, simplify=F)


# (B) Bayesian regression: use free-intercept hierarchical normal model to regress
#     mean burden per individual on life-history traits, using age as an offset
#     (for model definition, see file scripts/BHN_Regression.stan)
cat("Using BHN model to regress mutation burden on life-history traits...\n")

# Regression on individual variables
bhn.rate.anage.single = suppressWarnings(sapply(col.names, function(name) {
    idx = !is.na(indiv.data[, name])
    species = as.factor(indiv.data$Species[idx])
    # Build multiplier matrix to ensure that initial values for beta_g are >0
    sign.matrix = cbind(1, ifelse(species.data[match(levels(species), species.data$Species),
                                               name] < 0, -1, 1))
    sampling(bhn.model,
             iter = ITER, seed = SEED, cores = CORES, refresh = 0,
             data = list(K = 2,
                         N = sum(idx),
                         G = length(levels(species)),
                         y = indiv.data$Mut_burden[idx],
                         X = cbind(1, indiv.data[idx, name]),
                         t = indiv.data$Age_years[idx],
                         group = as.integer(species)),
             init = function() {
                 # Define random initial values for beta_g within model constraints
                 list(beta_g = sign.matrix * matrix(runif(2 * length(levels(species)), 1, 100),
                                                    ncol=2))
             })
}, simplify=F))

# Regression on variables paired with 1/Lifespan
bhn.rate.anage.paired = suppressWarnings(sapply(col.names[-1], function(name) {
    idx = !is.na(indiv.data[, name])
    species = as.factor(indiv.data$Species[idx])
    # Build multiplier matrix to ensure that initial values for beta_g are >0
    sign.matrix = cbind(1, 1, ifelse(species.data[match(levels(species), species.data$Species),
                                                  name] < 0, -1, 1))
    sampling(bhn.model,
             iter = ITER, seed = SEED, cores = CORES, refresh = 0,
             data = list(K = 3,
                         N = sum(idx),
                         G = length(levels(species)),
                         y = indiv.data$Mut_burden[idx],
                         X = cbind(1, indiv.data$Inverse_Lifespan80[idx], indiv.data[idx, name]),
                         t = indiv.data$Age_years[idx],
                         group = as.integer(species)),
             init = function() {
                 # Define random initial values for beta_g within model constraints
                 list(beta_g = sign.matrix * matrix(runif(3 * length(levels(species)), 1, 100),
                                                    ncol=3))
             })
}, simplify=F))


# (5) Regress mean substitution rate per species on inverse lifespan
#     Models: linear model (LM), phylogenetic least-squares (PGLS) model
#     (constrained intercept)

# (A) LM regression: use constrained (zero-intercept) simple LM to
#     regress mean mutation rates per species on 1/Lifespan and log10(Mass)
#     (two versions: [1] including all species; [2] excluding mouse and rat)
cat("Using LM and PGLS models to regress mean mutation rate on life-history traits...\n")
idx = !(species.data$Species %in% c("mouse", "rat"))
lm.rate.lifespan80.1 = lm(Mut_rate ~ Inverse_Lifespan80 - 1, species.data)
lm.rate.lifespan80.2 = lm(Mut_rate ~ Inverse_Lifespan80 - 1, species.data[idx, ])
lm.rate.mass.1 = lm(Mut_rate ~ Log10_Mass, species.data)
lm.rate.mass.2 = lm(Mut_rate ~ Log10_Mass, species.data[idx, ])

# (B) PGLS regression: use constrained phylogenetic generalised least-squares
#     model to regress mean mutation rates per species on 1/Lifespan and
#     log10(Mass), while accounting for phyl. relationships among species
#     (two versions: [1] including all species; [2] excluding mouse and rat)
pgls.rate.lifespan80.1 = pgls(Mut_rate ~ Inverse_Lifespan80 - 1,
                              comparative.data(species.tree, species.data,
                                               "Latin_name", vcv=T, na.omit=F))
pgls.rate.lifespan80.2 = pgls(Mut_rate ~ Inverse_Lifespan80 - 1,
                              comparative.data(species.tree, species.data[idx, ],
                                               "Latin_name", vcv=T, na.omit=F))
pgls.rate.mass.1 = pgls(Mut_rate ~ Log10_Mass,
                        comparative.data(species.tree, species.data,
                                           "Latin_name", vcv=T, na.omit=F))
pgls.rate.mass.2 = pgls(Mut_rate ~ Log10_Mass,
                        comparative.data(species.tree, species.data[idx, ],
                                           "Latin_name", vcv=T, na.omit=F))


# (6) Bootstrapped regression of substitution rate on 1/Lifespan,
#     using random combinations of published lifespan estimates
#     Model: linear mixed-effects model (constrained intercept)

cat("Performing bootstrapped regression of mutation rate on published lifespan estimates...\n")
lme.rate.lifespan.boot = matrix(NA, ITER, 2, dimnames=list(NULL, c("k", "FVE")))
set.seed(SEED)

for (i in 1:ITER) {
    
    # In each bootstrap replicate, sample a random lifespan estimate per species
    boot.lifespan = sapply(species.data$Species, function(species) {
        idx = which(lifespans$Species == species & lifespans$Measure == "Max_longevity")
        lifespans$Lifespan[sample(idx, 1)]
    })
    sample.data$BS_inverse_lifespan = 1 / boot.lifespan[sample.data$Species]
    
    # Perform LME regression and save fixed effects and FVE
    boot.lme = lme(fixed = Mut_rate ~ BS_inverse_lifespan - 1,
                   random = ~ BS_inverse_lifespan - 1 | Species/Individual,
                   weights = varIdent(form = ~ 1 | Species),
                   data = sample.data,
                   method = "ML")
    
    lme.rate.lifespan.boot[i, ] = c(boot.lme$coefficients$fixed,
                                    r2.lme(boot.lme, as.matrix(1 / boot.lifespan),
                                           species.data$Mut_rate))
}


# (7) Plot mutation rate regression results
cat("\nPlotting regression results...\n")

# (A) Plot fitted LME/BHN models for mutation rates/burdens
cairo_pdf(OUTPUT$RATES.PDF, 12, 6, onefile=T)
par(mfrow=c(1, 2), oma=c(0, 1, 0, 0))
col = "dodgerblue3"

# Plot rate vs 1/lifespan regressions
for (i in 1:length(lme.rates.lifespan80)) {
    
    # Obtain k and FVE (fraction of inter-species variance explained = R-squared)
    name = names(lme.rates.lifespan80)[i]
    if (name == "mtDNA_rate") {
        k.est = round(c(confint(lme.rates.lifespan80[[i]])[1], coef(lme.rates.lifespan80[[i]]),
                        confint(lme.rates.lifespan80[[i]])[2]), 2)
        fve = round(r2.lm(lme.rates.lifespan80[[i]]), 2)
    } else {
        k.est = round(intervals(lme.rates.lifespan80[[i]], which="fixed")$fixed, 2)
        fve = round(r2.lme(lme.rates.lifespan80[[i]],
                           model.matrix(~ Inverse_Lifespan80 - 1, species.data),
                           species.data[, name]), 2)
    }
    
    # Plot LME regression
    ci.x = seq(0, 95, len=200)
    rates = species.data[, name]
    plot(species.data$Lifespan80, rates,
         xlim=c(0, 90), ylim=c(0, 1.1 * max(rates, na.rm=T)), xaxs="i", yaxs="i",
         pch=16, type="n", cex.main=1.1,
         xlab="Lifespan (years)", ylab="Mutations per genome per year",
         main=ifelse(name == "mtDNA_rate", "\n\nSimple linear model",
                     "\n\nLinear mixed-effects model"))
    polygon(c(ci.x, rev(ci.x)), c(2*k.est[2]/ci.x, rev(0.5*k.est[2]/ci.x)),
            col=alpha(col, 0.15), border=NA)
    polygon(c(ci.x, rev(ci.x)), c(k.est[1]/ci.x, rev(k.est[3]/ci.x)),
            col=alpha(col, 0.25), border=NA)
    lines(ci.x, k.est[2]/ci.x, col=col, lwd=3)
    points(species.data$Lifespan80, rates, pch=16, cex=1.2)
    text(species.data$Lifespan80, rates + max(rates, na.rm=T) * 0.025,
         labels=SP.LAB[species.data$Species], cex=0.5)
    mtext(c(paste(ifelse(name == "mtDNA_rate", "LM: ", "LME: "), name, "~ 1/Lifespan80 – 1"),
            paste0("k = ", k.est[2], " (95% CI: ", k.est[1], "–", k.est[3], ")"),
            paste("FVE =", fve)),
          side=3, at=25, line=-(3:5), adj=0, cex=0.9)
    
    # Obtain k and FVE for BHN model
    k.est = round(quantile(extract(bhn.rates.lifespan80[[i]], "beta")[[1]],
                           c(0.025, 0.5, 0.975)), 2)
    fve = round(r2.stan(bhn.rates.lifespan80[[i]],
                        model.matrix(~ Inverse_Lifespan80 - 1, species.data),
                        species.data[, name]), 2)
    
    # Plot BHN regression
    ci.x = seq(0, 95, len=200)
    plot(species.data$Lifespan80, rates,
         xlim=c(0, 90), ylim=c(0, 1.1 * max(rates, na.rm=T)), xaxs="i", yaxs="i",
         pch=16, type="n", xlab="Lifespan (years)", ylab=paste("Mutations per genome per year"),
         main="\n\nBayesian hier. normal model (individual means)", cex.main=1.1)
    polygon(c(ci.x, rev(ci.x)), c(2*k.est[2]/ci.x, rev(0.5*k.est[2]/ci.x)),
            col=alpha(col, 0.15), border=NA)
    polygon(c(ci.x, rev(ci.x)), c(k.est[1]/ci.x, rev(k.est[3]/ci.x)),
            col=alpha(col, 0.25), border=NA)
    lines(ci.x, k.est[2]/ci.x, col=col, lwd=3)
    points(species.data$Lifespan80, rates, pch=16, cex=1.2)
    text(species.data$Lifespan80, rates + max(rates, na.rm=T) * 0.025,
         labels=SP.LAB[species.data$Species], cex=0.5)
    mtext(c(paste("BHN: ", name, "~ 1/Lifespan80 – 1"),
            paste0("k = ", k.est[2], " (95% CI: ", k.est[1], "–", k.est[3], ")"),
            paste("FVE =", fve)),
          side=3, at=25, line=-(3:5), adj=0, cex=0.9)
    title(paste("Regression of", name, "on 1/Lifespan80 (constrained)"), outer=T, line=-1.8)
}

# Plot rate vs AnAge regressions (single)
for (i in 1:length(lme.rate.anage.single)) {
    
    # Obtain k, FVE and 95% bootstrap intervals for LME
    name = names(lme.rate.anage.single)[i]
    predictor = species.data[, name]
    ci.x = seq(min(0, 1.3 * min(predictor, na.rm=T)), 1.2 * max(predictor, na.rm=T), len=80)
    y.samples = apply(lme.rate.anage.single.bs[[i]], 1, function(k) {
        k[1] + k[2] * ci.x
    })
    ci.y = t(apply(y.samples, 1, quantile, probs=c(0.025, 0.975), na.rm=T))
    k.est = round(intervals(lme.rate.anage.single[[i]], which="fixed")$fixed, 2)
    fve = round(r2.lme(lme.rate.anage.single[[i]],
                       model.matrix(~ predictor), species.data$Mut_rate[!is.na(predictor)]), 2)
    
    # Plot LME regression
    plot(predictor, species.data$Mut_rate,
         xlim=c(min(0, 1.2 * min(predictor, na.rm=T)), 1.1 * max(predictor, na.rm=T)),
         ylim=c(0, 1.1 * max(species.data$Mut_rate, na.rm=T)),
         xaxs="i", yaxs="i", pch=16, type="n",
         xlab=name, ylab=paste("Mutations per genome per year"),
         main="\n\nLinear mixed-effects model", cex.main=1.1)
    polygon(c(ci.x, rev(ci.x)), c(ci.y[, 1], rev(ci.y[, 2])), col=alpha(col, 0.25), border=NA)
    abline(k.est[, 2], col=col, lwd=3)
    points(predictor, species.data$Mut_rate, pch=16, cex=1.2)
    text(predictor, species.data$Mut_rate + max(species.data$Mut_rate) * 0.025,
         labels=SP.LAB[species.data$Species], cex=0.5)
    mtext(c(paste("LME:  Mut_rate ~", name),
            paste0("Intercept = ", k.est[1, 2], " (95% CI: ", k.est[1, 1], "–", k.est[1, 3], ")"),
            paste0("Slope = ", k.est[2, 2], " (95% CI: ", k.est[2, 1], "–", k.est[2, 3], ")"),
            paste("FVE =", fve)),
          side=3, at=0.12 * max(predictor, na.rm=T), line=-(3:6), adj=0, cex=0.9)
    
    # Obtain k, FVE and 95% credible intervals for BHN
    ci.x = seq(min(0, 1.3 * min(predictor, na.rm=T)), 1.2 * max(predictor, na.rm=T), len=40)
    k.samples = extract(bhn.rate.anage.single[[i]], "beta")[[1]]
    y.samples = apply(k.samples, 1, function(k) {
        k[1] + k[2] * ci.x
    })
    ci.y = t(apply(y.samples, 1, quantile, probs=c(0.025, 0.975)))
    k.est = round(t(apply(k.samples, 2, quantile, c(0.025, 0.5, 0.975))), 2)
    fve = round(r2.stan(bhn.rate.anage.single[[i]],
                        model.matrix(~ predictor), species.data$Mut_rate[!is.na(predictor)]), 2)
    
    # Plot BHN regression
    plot(predictor, species.data$Mut_rate,
         xlim=c(min(0, 1.2 * min(predictor, na.rm=T)), 1.1 * max(predictor, na.rm=T)),
         ylim=c(0, 1.1 * max(species.data$Mut_rate, na.rm=T)),
         xaxs="i", yaxs="i", pch=16, type="n",
         xlab=name, ylab=paste("Mutations per genome per year"),
         main="\n\nBayesian hier. normal model (individual means)", cex.main=1.1)
    polygon(c(ci.x, rev(ci.x)), c(ci.y[, 1], rev(ci.y[, 2])), col=alpha(col, 0.25), border=NA)
    abline(k.est[, 2], col=col, lwd=3)
    points(predictor, species.data$Mut_rate, pch=16, cex=1.2)
    text(predictor, species.data$Mut_rate + max(species.data$Mut_rate) * 0.025,
         labels=SP.LAB[species.data$Species], cex=0.5)
    mtext(c(paste("BHN:  Mut_rate ~", names(bhn.rate.anage.single)[i]),
            paste0("Intercept = ", k.est[1, 2], " (95% CI: ", k.est[1, 1], "–", k.est[1, 3], ")"),
            paste0("Slope = ", k.est[2, 2], " (95% CI: ", k.est[2, 1], "–", k.est[2, 3], ")"),
            paste("FVE =", fve)),
          side=3, at=0.12 * max(predictor, na.rm=T), line=-(3:6), adj=0, cex=0.9)
    title(paste("Regression of Mut_rate on", names(bhn.rate.anage.single)[i]), outer=T, line=-1.8)
}

# Plot comparative of FVE values for life-history regressions
# Obtain LME and BHN FVEs
lme.rate.anage.single.fve = sapply(names(lme.rate.anage.single), function(name) {
    r2.lme(lme.rate.anage.single[[name]], cbind(1, species.data[, name]), species.data$Mut_rate)
})
bhn.rate.anage.single.fve = sapply(names(bhn.rate.anage.single), function(name) {
    r2.stan(bhn.rate.anage.single[[name]], cbind(1, species.data[, name]), species.data$Mut_rate)
})
lme.rate.anage.paired.fve = sapply(names(lme.rate.anage.paired), function(name) {
    r2.lme(lme.rate.anage.paired[[name]],
           cbind(1, species.data$Inverse_Lifespan80, species.data[, name]),
           species.data$Mut_rate)
})
bhn.rate.anage.paired.fve = sapply(names(bhn.rate.anage.paired), function(name) {
    r2.stan(bhn.rate.anage.paired[[name]],
            cbind(1, species.data$Inverse_Lifespan80, species.data[, name]),
            species.data$Mut_rate)
})

# For each pair of variables, plot FVE for individual and paired regressions
cols = c("firebrick", "dodgerblue4")
par(mfrow=c(1, 1), mar=c(3.5, 5, 4, 2.5), mgp=c(3.4, 0.8, 0))
idx = c("1/Lifespan"="Inverse_Lifespan80", "Litter size"="Litter_size",
        "log10(Adult mass)"="Log10_Mass", "log10(BMR)"="Log10_Metabolic_rate",
        "BMR residuals"="Residuals_BMR_Mass")
b = barplot(rbind(lme.rate.anage.single.fve[idx], lme.rate.anage.paired.fve[idx]),
            beside=T, col=cols, border="white", ylim=c(0, 1), las=1,
            cex.names=1e-10, cex.lab=1.25, cex.axis=1.1, cex.main=1.3, xpd=F,
            ylab="Fraction of inter-species variance explained (FVE)",
            main="Comparison of FVEs for different explanatory variables (life-history traits)",
            panel.first=abline(h=seq(0, 1, 0.2), col="grey85", lty=1))
axis(1, c(b[1, 1], colMeans(b[, -1])), line=-1, tick=F,
     cex.axis=1.2, mgp=c(0, 2, 0), labels=names(idx))
legend("topleft", legend=c("Individually", "Combined with 1/Lifespan"), bty="n", pch=15,
       cex=1.15, pt.cex=2, col=cols, inset=c(0.028, 0), xpd=NA)
box()

invisible(dev.off())


# (B) Plot allometric regressions among mutation rate, mass and lifespan
cairo_pdf(OUTPUT$ALLOM.PDF, 7.5, 7)
par(mfrow=c(2, 2), mar=c(4,4,1,1), mgp=c(2.3,0.7,0))
col = "dodgerblue3"

plot.regr = function(x, y, model, ci.x, newdata, xlab, ylab) {
    ci.y = predict(model, newdata=newdata, interval="confidence", level=0.95)
    plot(y ~ x, pch=16, type="n", xlab=xlab, ylab=ylab, cex.axis=0.9,
         xlim=c(min(x, na.rm=T) - max(x, na.rm=T) * 0.1, max(x, na.rm=T) * 1.1),
         ylim=c(min(y, na.rm=T) - max(y, na.rm=T) * 0.1, max(y, na.rm=T) * 1.1))
    polygon(c(ci.x, rev(ci.x)), c(ci.y[,2], rev(ci.y[,3])), col=alpha(col, 0.25), border=NA)
    abline(model, lwd=2.5, col=alpha(col, 0.7))
    points(y ~ x, pch=16, cex=1.4)
    pval = summary(model)$coefficients[2, 4]
    legend("topright", bty="n", inset=c(0.035, 0.03), cex=1.1,
           c(paste("FVE =", round(r2.lm(model), 2)),
             paste("P =", ifelse(pval < 1e-5, paste0(round(pval / 1e-6, 1), "e-6"),
                                 round(pval, 3)))))
    #print(model$call); cat("P =", summary(model)$coefficients[2, 4], "\n")
}

# log10(rate) vs log10(mass)
y = species.data$Log10_Mut_rate
x = species.data$Log10_Mass
ylab = "log10(Mutation rate)"
xlab = "log10(Adult mass)"
ci.x = seq(min(x, na.rm=T) - 1, max(x, na.rm=T) + 1, len=100)
newdata = data.frame("Log10_Mass"=ci.x)
model = lm.allom.rate.mass
plot.regr(x, y, model, ci.x, newdata, xlab, ylab)

# log10(rate) vs log10(lifespan)
y = species.data$Log10_Mut_rate
x = species.data$Log10_Lifespan80
ylab = "log10(Mutation rate)"
xlab = "log10(Lifespan)"
model = lm.allom.rate.lifespan
ci.x = seq(min(x, na.rm=T) - 1, max(x, na.rm=T) + 1, len=100)
newdata = data.frame("Log10_Lifespan80"=ci.x)
plot.regr(x, y, model, ci.x, newdata, xlab, ylab)

# rate-mass residuals vs lifespan-mass residuals
y = species.data$Residuals_Rate_Mass
x = species.data$Residuals_Lifespan_Mass
ylab = "Mutation rate residual [mass-adjusted rate]"
xlab = "Lifespan residual [mass-adjusted lifespan]"
model = lm.resid.rate.lifespan
ci.x = seq(min(x, na.rm=T) - 1, max(x, na.rm=T) + 1, len=100)
newdata = data.frame("Residuals_Lifespan_Mass"=ci.x)
plot.regr(x, y, model, ci.x, newdata, xlab, ylab)

# rate-lifespan residuals vs mass-lifespan residuals
y = species.data$Residuals_Rate_Lifespan
x = species.data$Residuals_Mass_Lifespan
ylab = "Mutation rate residual [lifespan-adjusted rate]"
xlab = "Adult mass residual [lifespan-adjusted mass]"
model = lm.resid.rate.mass
ci.x = seq(min(x, na.rm=T) - 1, max(x, na.rm=T) + 1, len=100)
newdata = data.frame("Residuals_Mass_Lifespan"=ci.x)
plot.regr(x, y, model, ci.x, newdata, xlab, ylab)

invisible(dev.off())


# (C) Plot comparative of fitted models
cairo_pdf(OUTPUT$COMP.PDF, 12, 6)
par(mfrow=c(1, 2), mar=c(4.5, 4.5, 4, 1.5), mgp=c(2.5, 0.7, 0))
cols = c("darkblue", "cornflowerblue", "darkorange1", "darkgoldenrod1", "red", "darkgreen")

# Plot fitted models for mutation rate vs 1/Lifespan
k.est = c(coef(lm.rate.lifespan80.1), coef(lm.rate.lifespan80.2),
          coef(pgls.rate.lifespan80.1), coef(pgls.rate.lifespan80.2),
          fixef(lme.rates.lifespan80$Mut_rate),
          median(extract(bhn.rates.lifespan80$Mut_burden, "beta")[[1]]))
names(k.est) = c("Simple linear model (species means, all species)",
                 "Simple linear model (species means, excl. mouse & rat)",
                 "Phylogenetic GLS model (species means, all species)",
                 "Phylogenetic GLS model (species means, excl. mouse & rat)",
                 "Linear mixed-effects model (all samples)",
                 "Bayesian hier. normal model (all samples)")
plot(species.data$Inverse_Lifespan80, species.data$Mut_rate,
     xlim=c(0, 1.1 * max(species.data$Inverse_Lifespan80)),
     ylim=c(0, 1.1 * max(species.data$Mut_rate)),
     xaxs="i", yaxs="i", pch=16, type="n", cex.main=1.2,
     xlab="1 / Lifespan", ylab="Mutations per genome per year",
     main="Regression of mutation rate on 1/Lifespan")
for (i in 1:length(k.est)) {
    abline(0, k.est[i], col=cols[i], lwd=3)
    text(0.1, seq(175, len=length(cols), by=-31)[i],
         col=cols[i], cex=0.75, font=2, adj=0, labels=names(k.est)[i])
    text(0.025, seq(800, len=length(cols), by=-43)[i],
         col=cols[i], font=2, cex=1.1, adj=0, labels=paste("k =", round(k.est[i], 2)))
}
points(species.data$Inverse_Lifespan80, species.data$Mut_rate, pch=16, cex=1.2)
text(species.data$Inverse_Lifespan80, species.data$Mut_rate + 0.025*max(species.data$Mut_rate),
     labels=SP.LAB[species.data$Species], cex=0.5)

# Plot fitted models for mutation rate vs log10(Mass)
k.est = rbind(coef(lm.rate.mass.1), coef(lm.rate.mass.2),
              coef(pgls.rate.mass.1), coef(pgls.rate.mass.2),
              fixef(lme.rate.anage.single$Log10_Mass),
              apply(extract(bhn.rate.anage.single$Log10_Mass, "beta")[[1]], 2, median))
plot(species.data$Log10_Mass, species.data$Mut_rate,
     xlim=c(0, 1.1 * max(species.data$Log10_Mass)),
     ylim=c(0, 1.1 * max(species.data$Mut_rate)),
     xaxs="i", yaxs="i", pch=16, type="n", cex.main=1.2,
     xlab="log10(Adult mass)", ylab="Mutations per genome per year",
     main="Regression of mutation rate on log10(Adult mass)")
for (i in 1:nrow(k.est)) {
    abline(k.est[i, ], col=cols[i], lwd=3)
    text(3.5, seq(800, len=length(cols), by=-43)[i],
         col=cols[i], font=2, cex=1.1, adj=0,
         labels=paste0("α = ", round(k.est[i, 1], 2), ",  β = ", round(k.est[i, 2], 2)))
}
points(species.data$Log10_Mass, species.data$Mut_rate, pch=16, cex=1.2)
text(species.data$Log10_Mass, species.data$Mut_rate + 0.025*max(species.data$Mut_rate),
     labels=SP.LAB[species.data$Species], cex=0.5)

invisible(dev.off())


# (D) Plot bootstrapped lifespan regression results
cairo_pdf(OUTPUT$LIFSP.PDF, 8, 6, onefile=T)
par(mar=c(4.5, 4.5, 4, 1.5), mgp=c(2.5, 0.7, 0))
x = seq(0, 150, len=200)
col = "dodgerblue3"

lifespan.mean = sapply(species.data$Species, function(sp) {
    mean(lifespans$Lifespan[lifespans$Species == sp & lifespans$Measure == "Max_longevity"])
})
lifespan.min = sapply(species.data$Species, function(sp) {
    min(lifespans$Lifespan[lifespans$Species == sp & lifespans$Measure == "Max_longevity"])
})
lifespan.max = sapply(species.data$Species, function(sp) {
    max(lifespans$Lifespan[lifespans$Species == sp & lifespans$Measure == "Max_longevity"])
})

plot(lifespan.mean, species.data$Mut_rate, pch=16,
     xlim=c(0, max(lifespan.max) * 1.05), ylim=c(0, max(species.data$Mut_rate) * 1.1),
     ylab="Mutations per genome per year", xlab="Maximum lifespan (years)",
     main="Bootstrapped regression of mutation rate on published lifespan estimates",
     type="n", xaxs="i", yaxs="i", cex.lab=1.25)
k.ci = round(summary(lme.rate.lifespan.boot[, 1])[c(1, 3, 6)], 2)
fve = round(summary(lme.rate.lifespan.boot[, 2])[c(1, 3, 6)], 2)
polygon(c(x, rev(x)), c(k.ci[1]/x, rev(k.ci[3]/x)), col=alpha(col, 0.3), border=NA)
lines(x, k.ci[2]/x, col=col, lwd=4)
segments(x0=lifespan.min, x1=lifespan.max, y0=species.data$Mut_rate, lwd=2.5)
points(lifespan.mean, species.data$Mut_rate, pch=16, cex=1.9)
text(lifespan.mean, species.data$Mut_rate + 26, SP.LAB[names(lifespan.mean)], cex=0.6)
mtext(c(paste(ITER, "bootstrap samples"),
        paste0("k:  Min = ", k.ci[1], ",  Median = ", k.ci[2], ",  Max = ", k.ci[3]),
        paste0("FVE:  Min = ", fve[1], ",  Median = ", fve[2], ",  Max = ", fve[3]),
        "Error bars: ranges of published lifespan estimates",
        "Blue area: range of k estimates",
        "Blue line: median k"),
      side=3, at=25, line=-c(3,4.5,6,9,10.5,12), adj=0, font=1, cex=1.25,
      col=rep(c("black", "grey50"), each=3))

hist(lme.rate.lifespan.boot[, 2], breaks=100, col="dodgerblue4", border="white",
     xlim=c(0.5,1), xaxs="i", yaxs="i",  cex.lab=1.25,
     xlab="Fraction of inter-species variance explained (FVE)",
     main=paste0("Distribution of regression FVEs (", ITER, " bootstrap samples)"))
abline(h=0)

invisible(dev.off())


cat("\nDone\n\n")
