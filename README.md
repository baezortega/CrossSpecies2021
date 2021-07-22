## Analysis scripts accompanying publication Cagan, Baez-Ortega _et al._, 2021

_Adrian Baez-Ortega  
Wellcome Sanger Institute  
2020–21_

This repository contains custom R and bash scripts which can be used to replicate the analyses presented in the article "UNPUBLISHED" (Cagan, Baez-Ortega _et al._, 2021).

The set-up and usage of the data and scripts are explained below.

If you use some of these methods for your own research, please use the following citation:

---

#### A. Cagan, A. Baez-Ortega _et al_. [UNPUBLISHED]

---

### Set-up and requirements

In order to run the scripts, the repository must first be cloned into a directory in your local machine (replace the destination path `~/Desktop/CrossSpecies2021` below with your own). This can be done from the terminal using the command below.

```
git clone https://github.com/baezortega/CrossSpecies2021.git ~/Desktop/CrossSpecies2021
```

The scripts are written to run on [**R**](https://www.r-project.org/) version **3.6.2** or later. You should be able to check your current version of R by running one of the commands below (depending on your installation):

```
R --version
/Library/Frameworks/R.framework/Resources/bin/R --version
```

The current R version is also shown when opening RStudio or the R Console.

The scripts require the following R packages:  [**`bbmle`**](https://cran.r-project.org/web/packages/bbmle/index.html),  [**`Biostrings`**](https://bioconductor.org/packages/release/bioc/html/Biostrings.html),  [**`caper`**](https://cran.r-project.org/web/packages/caper/index.html), [**`dNdScv`**](https://github.com/im3sanger/dndscv),  [**`emdbook`**](https://cran.r-project.org/web/packages/emdbook/index.html), 
[**`GenomicRanges`**](https://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html),  [**`MASS`**](https://cran.r-project.org/web/packages/MASS/index.html),  [**`nlme`**](https://cran.r-project.org/web/packages/nlme/index.html),  [**`RColorBrewer`**](https://cran.r-project.org/web/packages/RColorBrewer/index.html),  [**`scales`**](https://cran.r-project.org/web/packages/scales/index.html),  [**`sigfit` (>=2.1)**](https://github.com/kgori/sigfit), 
[**`SigProfilerMatrixGeneratorR`**](https://github.com/AlexandrovLab/SigProfilerMatrixGeneratorR).

In addition, note that [**Python 3**](https://www.python.org) and the Python package [**`SigProfilerMatrixGenerator`**](https://osf.io/s93d5/wiki/home/) are required before installing the `SigProfilerMatrixGeneratorR` R package, which is needed for one of the analysis steps below (Step 3). In addition, human, dog, mouse and rat reference genomes for `SigProfilerMatrixGenerator` need to be [installed manually](https://osf.io/s93d5/wiki/1.%20Installation%20-%20Python/). However, note that this particular step is not required for any subsequent data analyses, and so **can be omitted** if you encounter difficulties in installing Python 3 or `SigProfilerMatrixGenerator`. You should be able to check your current version of Python by running the command `python --version`.

Although care has been taken to make the code distribution-independent, some of the scripts may only work on Unix/MacOS systems, and may need to be modified in order to run on Windows systems.

Please note that some of the scripts may take considerable time to run, and some of the intermediate files generated will occupy up to a few gigabytes.

**Finally, before running the steps indicated below, it is necessary to open the terminal and navigate to the cloned directory, `CrossSpecies2021`, using the command below (replace the path below with your own).**

```
cd ~/Desktop/CrossSpecies2021
```

---

### Step 0: Setting up project data and directories

Before commencing the analyses, this script initialise the directory structure, and downloads the project data files from [Zenodo](http://zenodo.org/record/XXXXXXXX) and the required reference genomes from Ensembl/NCBI.

* This step uses the `curl` command and requires Internet access.
* The output of this step is required for subsequent steps.
* The approximate run time of this step is **30 minutes**.
* The output files produced include all the necessary data files (`data/original` folder).

This step is performed by the R script `0_Setup.sh`, which is located in the [`scripts`](scripts) directory and can be run from the terminal using the `bash` command as follows.

```
bash scripts/0_Setup.sh
```

---

### Step 1: Identification of likely polyclonal samples

This step identifies samples that are likely polyclonal (based on mutation burden and VAF), and produces a list of samples to be excluded from subsequent analyses.

* This step requires data files produced in Step 0.
* The output of this step is required for subsequent steps.
* The approximate run time of this step is **1 hour**.
* The output files produced include VAF histograms (`output/Clonality_VAF_Histograms.pdf`), a table of clonality-based sample filtering results (`output/Clonality_Results_Table.txt`), a list of samples to be excluded from all analyses (`data/processed/SamplesToExclude.txt`), and a list of samples to be excluded from mtDNA analyses only (`data/processed/SamplesToExclude_mtDNA.txt`).

This step is performed by the R script `1_Clonality.R`, which is located in the [`scripts`](scripts) directory and can be run either from RStudio (following the instructions at the beginning of the script), or from the terminal using the `Rscript` command as follows.

```
Rscript scripts/1_Clonality.R
```

 **NB.** Depending on your installation, you may have to replace `Rscript` with `/Library/Frameworks/R.framework/Resources/bin/Rscript`.

---

### Step 2: Mutational spectra of single-base substitutions (SBS)

This step produces mutational spectra of the somatic substitutions in each sample, and in all samples from each species.

* This step requires R packages `Biostrings`, `GenomicRanges` and `sigfit`.
* This step requires data files produced in Steps 0 and 1.
* The output of this step is required for subsequent steps.
* The approximate run time of this step is **1.5 hours**.
* The output files produced include plots of SBS spectra per sample and per species (`output/Spectra_Subs_Sample.pdf`, `output/Spectra_Subs_Species.pdf`), and an RData file containing mutational catalogues and opportunities (`data/processed/Catalogues_Opportunities.RData`).

This step is performed by the R script `2_Spectra_SBS.R`, which is located in the [`scripts`](scripts) directory and can be run either from RStudio (following the instructions at the beginning of the script), or from the terminal using the `Rscript` command as follows.

```
Rscript scripts/2_Spectra_SBS.R
```

**NB.** If running the script triggers memory allocation errors, the amount of virtual memory assigned to R can be increased using:

```
env R_MAX_VSIZE=100Gb Rscript scripts/2_Spectra_SBS.R
```

---

### Step 3: Mutational spectra of insertions/deletions (indels)

This step produces mutational spectra of the somatic indels in each sample, and in all samples from each species.

* This step requires the `SigProfilerMatrixGenerator` Python 3 package, and the `SigProfilerMatrixGeneratorR` R package.
* This step requires data files produced in Step 1.
* The output of this step is **not** required for subsequent steps.
* The approximate run time of this step is **10 minutes**.
* The output files produced include plots of indel spectra per sample and per species (folders `output/Spectra_Indels_Sample` and `output/Spectra_Indels_Species`).

This step is performed by the R script `3_Spectra_Indels.R`, which is located in the [`scripts`](scripts) directory and can be run either from RStudio (following the instructions at the beginning of the script), or from the terminal using the `Rscript` command as follows.

```
Rscript scripts/3_Spectra_Indels.R
```

---

### Step 4: Inference of mutational signatures from somatic substitutions


This step infers mutational signatures from the catalogues of somatic substitutions in each species; performs a cross-species analysis of signature SBSB; and examines the prevalence of colibactin and APOBEC mutagenesis in non-human samples.

* This step requires the `sigfit` R package (v2.1 or higher).
* This step requires data files produced in Steps 1 and 2.
* The output of this step is required for subsequent steps.
* The approximate run time of this step is **1.5 hours** (on 4 CPUs).
* Parts of this step are run in parallel; the number of available CPUs is detected automatically.
* The output files produced include plots of mutational signatures and exposures (folder `output/Signature_Extraction_Definitive`), plots of signature SBSB as inferred from the mutations in each species (`output/SBSB_Per_Species.pdf`), results from the analysis of colibactin and APOBEC prevalence (`output/Colibactin_Exposure.pdf`, `output/Colibactin_APOBEC_Tests.txt`), and RData files containing mutational signatures and signature exposures and mutation burdens per sample (`data/processed/Burden_Exposures.RData`, `data/processed/Signatures_Definitive.RData`).

This step is performed by the R script `4_Signatures.R`, which is located in the [`scripts`](scripts) directory and can be run either from RStudio (following the instructions at the beginning of the script), or from the terminal using the `Rscript` command as follows.

```
Rscript scripts/4_Signatures.R
```

---

### Step 5: Inference of allele-specific copy number in chromosome-level assemblies

This step infers segments of total and allele-specific copy number for samples in species with chromosome-level genome assemblies.

Note that this step takes a very long time to run. However, its output is not required for subsequent analyses, and so it **can be omitted** if necessary.

* This step requires R packages `bbmle`, `emdbook`, `GenomicRanges` and `MASS`.
* This step requires data files produced in Step 1.
* The output of this step is **not** required for subsequent steps.
* The approximate run time of this step is **250 hours**.
* The output files produced include plots and tables of copy number segments per sample (folder `output/Copy_Number`).

This step is performed by the R script `5_Copy_Number.R`, which is located in the [`scripts`](scripts) directory and can be run either from RStudio (following the instructions at the beginning of the script), or from the terminal using the `Rscript` command as follows.

```
Rscript scripts/5_Copy_Number.R
```

**NB.** The run time of this step can be reduced substantially by limiting the analysis to species with copy number changes (see line 65 in script `5_CopyNumber.R`).

---

### Step 6: Inference of dN/dS ratios in species with genome annotation

This step calculate the ratio between non-synonymous and synonymous mutation rates (dN/dS) from somatic substitutions in each species.

* This step requires R packages `dndscv` and `GenomicRanges`.
* This step requires data files produced in Step 1.
* The output of this step is **not** required for subsequent steps.
* The approximate run time of this step is **10 minutes**.
* The output files produced include a plot of exome-wide dN/dS per species (`output/dNdS_Species.pdf`) and a table of coding mutation counts per sample (`output/Coding_Mutation_Counts.txt`).

This step is performed by the R script `6_Selection.R`, which is located in the [`scripts`](scripts) directory and can be run either from RStudio (following the instructions at the beginning of the script), or from the terminal using the `Rscript` command as follows.

```
Rscript scripts/6_Selection.R
```

---

### Step 7: Calculation of somatic mutation burdens and rates per sample

This step calculates corrected mutation burdens, rates for somatic substitutions, indels, mtDNA mutations and signature-specific mutations.

* This step requires the `scales` R package.
* This step requires data files produced in Steps 1 and 2.
* The output of this step is required for subsequent steps.
* The approximate run time of this step is **1 minute**.
* The output files produced include plots, tables and RData files of mutation burdens, rates and end-of-lifespan burdens per sample (`output/Burden_Rate_ELB.pdf`, `output/Burden_Rate_ELB.txt`, `data/processed/Burdens_Rates.RData`).

This step is performed by the R script `7_Burdens.R`, which is located in the [`scripts`](scripts) directory and can be run either from RStudio (following the instructions at the beginning of the script), or from the terminal using the `Rscript` command as follows.

```
Rscript scripts/7_Burdens.R
```

---

### Step 8: Regression analyses of somatic mutation burdens and rates

This step applies a range of regression models (simple linear, linear mixed-effects (LME), hierarchical Bayesian, allometric, bootstrapped LME, and phylogenetic generalised least-squares) to quantify the associations between somatic mutation burdens/rates and several biological variables.

* This step requires R packages `caper`, `nlme`, `RColorBrewer`, `scales` and `rstan` (installed with `sigfit`).
* This step requires data files produced in Steps 1, 7.
* The approximate run time of this step is **2 hours** (on 4 CPUs).
* Parts of this step are run in parallel; the number of available CPUs is detected automatically.
* The output files produced include plots of the fitted models for each regression analysis (`output/Regression_Burden-Age_LM.pdf`, `output/Regression_Rates_LME-BHN.pdf`, `output/Regression_Allometric.pdf`, `output/Regression_Lifespan_Bootstrap.pdf`, `output/Regression_Model_Comparison.pdf`)

This step is performed by the R script `8_Regressions.R`, which is located in the [`scripts`](scripts) directory and can be run either from RStudio (following the instructions at the beginning of the script), or from the terminal using the `Rscript` command as follows.

```
Rscript scripts/8_Regressions.R
```
