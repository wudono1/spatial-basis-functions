
## CORGI (Change Of Resolution in GLM Inference)

This repository contains an implementation in R and C++ of spatial-statistical downscaling using basis-function models for biodiversity modelling. The code in this repository is based on 

Zheng, X., Cressie, N., Clarke, D. A., McGeoch, M. A., and Zammit-Mangion, A. (2025). 
"Spatial-statistical downscaling with uncertainty quantification in biodiversity modelling". To appear in *Methods in Ecology and Evolution*.
<!-- NIASRA Working Paper Series 06-24. URL: [https://www.uow.edu.au/niasra/publications/](https://www.uow.edu.au/niasra/publications/) -->

and its corresponding repository [here](https://github.com/xzheng42/corgi-examples).

The basis-function implementation draws from

Cressie, N., Sainsbury-Dale, M., and Zammit-Mangion, A. (2022). "Basis-Function Models in Spatial Statistics"


### A two-stage protocol for biodiversity inference

CORGI represents the implementation of a two-stage protocol for biodiversity inference. The first stage 
of CORGI is used to downscale coarse-resolution environmental covariates to fine resolutions. The second 
stage of CORGI involves a GLM-based biodiversity model that uses the downscaled covariates and incorporates the downscaling uncertainties.

### Installing and using the **spds** package


We use [**spds**](https://github.com/xzheng42/corgi-examples/tree/main/spds) (currently developer's version) for 
spatial-statistical downscaling with uncertainty quantification.

You can install the package with **devtools**
```
devtools::install_github("xzheng42/corgi-examples", subdir = "spds")
library(spds)
```

The main functions of the package are `gpbau` and `predict.gpbau`:

- `gpbau` fits a GPB (Gaussian process over Basic Areal Units) via Markov chain Monte Carlo (MCMC).
- `predict.gpbau` (or simply `predict`) generates posterior predictive samples of the downscaled covariate given a gpbau object.

Notes:

- The current version was tested on macOS 10.15.7 under R version 4.2.2 and on Fedora Linux 38 under R versions 4.3.3.

- Instead of installing the package, you can load the package in R via devtools::load_all("~/path/to/spds/").

### Workflow to reproduce numerical results

R scripts to reproduce results in the manuscript are available in
[*data-examples/*](https://github.com/xzheng42/corgi-examples/tree/main/data-examples).

Notes: 

- The master script (`run_all_rscripts.R`) has instructions for running the code. In each R script, users need to specify the input and/or output directories and may need to install additional packages.

- Data can be assessed in [*data/*](https://github.com/xzheng42/corgi-examples/tree/main/data). 
  The folder contains Bunger Hills' shape files (bh), landscape data (dem), climate reanalysis (ear5-land), and moss data (moss).

# Reporting bugs and issues

If you find a bug or issue, please use the Github Issues (i.e., click the `Issues` tab, click the green `New issue` button and create a new issue). It is much appreciated if you could include the following in the issue: 

- A description of the bug or issue;
- A minimal reproducible example (e.g., source code and data files) with your comments.
