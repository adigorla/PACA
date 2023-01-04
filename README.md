
# PACA

<!-- badges: start -->
![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)
![release: v0.1.0](https://img.shields.io/badge/release-v0.1.0-green)
![coverage: 100%](https://img.shields.io/badge/coverage-80%25-brightgreen)
![docs: in-progress](https://img.shields.io/badge/docs-in--progress-yellow)
<!-- badges: end -->

_**Documenation development in-progress**_

Phenotype Aware Components Analysis (**PACA**) is a contrastive learning approach leveraging canonical correlation analysis to robustly capture weak sources of subphenotypic variation. PACA can be used to define *de novo* subtypes that are more likely to reflect molecular heterogeneity, especially in challenging cases where the phenotypic heterogeneity may be masked by a myriad of strong unrelated effects in the data.

## Installation

**PACA** is implemented as a **R** packages which depends on the following :

* Rcpp
* RcppArmadillo
* rsvd
* stats
* gtools


You can install **PACA** using *devtools*:

``` r
devtools::install_github("Adigorla/PACA")
```

Please see troubleshooting at the bottom for compilation issues.

## Usage

### autoPACA

The `autoPACA` algorithm runs the basic `PACA` algorithm after automatically estimating the number of shared dimensions `k` to be removed and returns the unique components in the cases. It chooses `k` which maximizes the variation unique to cases in a given case/control dataset. The 

``` r
# load package
library(PACA)

# load data
X <- read.table("case_data1.txt")
Y <- read.table("control_data1.txt")
```
The input data, `X` & `Y` needs to be of form samples-by-features (NxM), where M >> N. The data also need to be standardized along the feature axis, e.g. gene quantile normalization for RNAseq data. 
NOTE: for all examples we assume the the number of samples, N, is the same for cases and controls for simplicity. In reality, the number if cases and controls can be different. `PACA` only requires the number of features, M, to be the same and alinged in the case/control data.

``` r
# run autoPACA
resPACA <- autoPACA(X, Y)

# number of the top shared components removed
print(resPACA$k)

# the dimension of the returned unique components of the cases
print(dim(resPACA$x))

```
Please refer to the [autoPACA man page](man/autoPACA.Rd) for more detailed usage information.

### rPACA

`rPACA` is a randomized version of the basic `PACA` algorithm. `rPACA` allows us to apply `PACA` in regimes where M << N, i.e., in cases where the number of samples is greater than the number of features.

``` r
# load package
library(PACA)

# load data
X <- read.table("case_data1.txt")
Y <- read.table("control_data1.txt")
```
The input data, `X` & `Y` needs to be of form samples-by-features (NxM), where M << N. The data also need to be standardized along the feature axis. Also note that the number of shared dimensions `k` to be removed needs to be user specified. This can be chosen by performing a grid search over a range of `k` and picking a value that maximizes some application specific metric or based on the user's domain specific knowledge.

``` r
# run rPACA
resPACA <- rPACA(X, Y, k = 10, nIter = 20, stepSize = 600, pcRank = 4)

# the dimension of the returned unique components of the cases
print(dim(resPACA$x))

```
`nIter`, `pcRank` and `stepSize` are optional params. However, the users needs to make sure to set `stepSize` to `stepSize < min({M, N}` and `k < stepSize-1`. Increasing `nIter` and/or `pcRank` empirically seems to increase the estimation accuracy of the randomized alogoritim, at the expense of increase runtime. 
Please refer to the [rPACA man page](man/rPACA.Rd) for more detailed usage information.

### PACA
One also has the option to run the `PACA` algorithm and set a user defined `k`. This would return the unique components in the cases at the user definde `k`.

``` r
# load package
library(PACA)

# load data
X <- read.table("case_data1.txt")
Y <- read.table("control_data1.txt")

# transpose input data NxM --> MxN and scale along sample axis
# DO NOT SKIP
inputDat <- transformCCAinput(X, Y, .center = TRUE, .scale = TRUE)
```
The input data, `X` & `Y` needs to be of form samples-by-features (NxM), where M >>N. The data also need to be standardized along the feature axis. Note that here, the user needs to explicitly run the `transformCCAinput` function, which transpose input data and scales along sample axis. This is a mandatory step to transform the data into the appropriate form for the core `PACA` algorithm. This method is called internally in `autoPACA` and `rPACA`.

``` r
# run PACA
resPACA <- PACA(X, Y, k = 10)

# the dimension of the returned unique components of the cases
print(dim(resPACA$x))

```
Please refer to the [PACA man page](man/PACA.Rd) and [transformCCAinput man page](man/transformCCAinput.Rd) and for more detailed usage information.

### Null Testing

The `nullEvalPACA` algorithm allows users to test for the statistical significance of the presence of subphenotypic variation unique to the cases, for a given fixed `k`. This procedure should be able to reject the null (no subphenotypic variation) when there is sufficently strong variation unique to the cases.

``` r
# load package
library(PACA)

# load data
X <- read.table("case_data1.txt")
Y <- read.table("control_data1.txt")

resNulltest <- nullEvalPACA(X, Y, k, nPerm = 100)

# p-value of rejecting the null
print(resNulltest$pval)
```
The input data, `X` & `Y` needs to be of form samples-by-features (NxM), where M >> N. The data also need to be standardized along the feature axis. Increasing `nPerm` increases the precision of the `pval` estimate.
Please refer to the [nullEvalPACA man page](man/nullEvalPACA.Rd) for more detailed usage information.

## Troubleshooting

If you are using a mac and having installation issues, try installing homebrew or xcode then reinstalling **Rcpp** and **RcppArmadillo**. 

#### R >= 4.0+ on M1/2 Macs
If you are having issues compiling R/Rcpp code on the newer ARM (M1/2) Mac hardware, make you have `gcc(11+)` and `llvm` installed using homebrew.
``` bash 
brew install gcc && brew install llvm 
```

Then update the `Makevars` file in the `~/.R/` directory to the following:
```
# custom G++ makevars 
# adapeted from here: https://stackoverflow.com/questions/65860439/installing-data-table-on-macos

GCC_LOC = /opt/homebrew/Cellar/gcc/11.2.0_3                    # UPATDTE & CHECK  path is valid
FLIBS=-L$(GCC_LOC)/lib/gcc/11 -L$(GCC_LOC)/lib -lgfortran -lm
CXX1X=$(GCC_LOC)/bin/g++-11
CXX98=$(GCC_LOC)/bin/g++-11
CXX11=$(GCC_LOC)/bin/g++-11
CXX14=$(GCC_LOC)/bin/g++-11
CXX17=$(GCC_LOC)/bin/g++-11

LLVM_LOC = /opt/homebrew/opt/llvm                              # UPATDTE & CHECK path is valid
CC=$(GCC_LOC)/bin/gcc-11 -fopenmp
CXX=$(GCC_LOC)/bin/g++-11 -fopenmp
CFLAGS=-g -O3 -Wall -pedantic -std=gnu99 -mtune=native -pipe
CXXFLAGS=-g -O3 -Wall -pedantic -std=c++11 -mtune=native -pipe
LDFLAGS=-L$(LLVM_LOC)/lib -Wl,-rpath,$(LLVM_LOC)/lib
RARM_LOC = /opt/R/arm64                                        # UPATDTE & CHECK path is valid
BREW_LOC = /opt/homebrew                                       # UPATDTE & CHECK path is valid
CPPFLAGS=-I$(LLVM_LOC)/include -I$(BREW_LOC)/include -I$(RARM_LOC)/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include
```
Make sure that the four "UPATDTE & CHECK path is valid" lines point to valid location on your machine. 

For all older versions of R and Intel Mac installation issues, please refer to the detailed instructions on the [The Coatless Professor](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/) website.

### Credits

Our CCA algorithm is entirely based on that of R package [**CONFINED**](https://github.com/cozygene/CONFINED) by Mike Thompson, plese refer to their GitHub repo for more information on the CCA implementation. 

