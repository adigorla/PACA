# PACA

<!-- badges: start -->
![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)
![release: v0.2.0](https://img.shields.io/badge/release-v0.2.0-green)
![coverage: 100%](https://img.shields.io/badge/coverage-80%25-brightgreen)
![docs: in-progress](https://img.shields.io/badge/docs-in--progress-yellow)
<!-- badges: end -->

>[!WARNING]  
>___Documenation development in-progress___

Phenotype Aware Components Analysis (**PACA**) is a contrastive learning approach leveraging canonical correlation analysis to robustly capture weak sources of subphenotypic variation. PACA can be used to define *de novo* subtypes that are more likely to reflect molecular heterogeneity, especially in challenging cases where the phenotypic heterogeneity may be masked by a myriad of strong unrelated effects in the data.

## Installation

**PACA** is implemented as a **R** packages which depends on the following :

* Rcpp
* RcppEigen
* rsvd
* stats

You can install **PACA** using *devtools*:

``` r
devtools::install_github("Adigorla/PACA")
```

Please see troubleshooting at the bottom for compilation issues.

## Usage

### PACA

The `paca` commmand runs the basic **PACA** algorithm after automatically estimating the number of shared dimensions `k` to be removed and returns the unique components in the cases. It chooses `k` which maximizes the variation unique to cases in a given case/control dataset. The 
``` r
# load package
library(PACA)

# load data
X <- read.table("case_data1.txt")
Y <- read.table("control_data1.txt")
# IF: input data NxM, transpose both matrices to MxN
#       
#       Y <- t(Y)
```
>[!NOTE]
>all PACA functions require the input matrics to be in from features-by-samples (MxN). So if input data is NxM, transpose both matrices to MxN
> ```r
> X <- t(X)
> Y <- t(Y)
> ```
The input data, `X` & `Y` needs to be of the form samples-by-features (NxM), where M >> N. The data also need to be standardized along the feature axis, e.g. quantile normalization for RNAseq data. 
NOTE: for all examples we assume the the number of samples, N, is the same for cases and controls for simplicity. In reality, the number if cases and controls can be different. `paca` only requires the number of features, M, to be the same and alinged in the case/control data.

``` r
# run PACA (and infer k)
PACA.res <- paca(X, Y)

# the dimension of the returned unique components of the cases
print(dim(PACA.res$x)) # Nxk

# the dimension of the corrected case data
print(dim(PACA.res$xtil)) # MxN

```

One also has the option to run the `paca` algorithm and set a user defined `k`. This would return the unique components in the cases at the user definde `k`. The input data, `X` & `Y` needs to be of form samples-by-features (NxM), where M >> N.  
``` r
# run PACA w/ fixed k
resPACA <- PACA(X, Y, k = 10)

# the dimension of the returned unique components of the cases
print(dim(PACA.res$x)) # Nx10

# the dimension of the corrected case data
print(dim(PACA.res$xtil)) # MxN
```

Please refer to the [PACA man page](man/PACA.Rd) for more detailed usage information.

### rPACA

`paca_r` is a randomized version of the basic `paca` algorithm. `paca_r` allows us to apply **PACA** in regimes where M << N, i.e., in cases where the number of samples is greater than the number of features.

``` r
# load package
library(PACA)

# load data
X <- read.table("case_data1.txt")
Y <- read.table("control_data1.txt")
```
The input data, `X` & `Y` needs to be of form samples-by-features (NxM), where M << N. The data also need to be standardized along the feature axis. Also note that the number of shared dimensions `k` to be removed needs to be user specified. This can be chosen by performing a grid search over a range of `k` and picking a value that maximizes some application specific metric or based on the user's domain specific knowledge.

``` r
# run paca_r
resPACA <- paca_r(X, Y, k = 10, niter = 20, batch = 600, rank = 4)

# the dimension of the returned unique components of the cases
print(dim(resPACA$x))

```
`niter`, `rank` and `batch` are optional params. However, the users needs to make sure to set `batch` to `batch < min({M, N}` and `k < batch-1`. Increasing `niter` and/or `rank` empirically seems to increase the estimation accuracy of the randomized alogoritim, at the expense of increase runtime. 
Please refer to the [rPACA man page](man/rPACA.Rd) for more detailed usage information.

### Null Testing

The `paca_null` algorithm allows users to test for the statistical significance of the presence of subphenotypic variation unique to the cases, for a given fixed `k`. This procedure should be able to reject the null (no subphenotypic variation) when there is sufficently strong variation unique to the cases.

``` r
# load package
library(PACA)

# load data
X <- read.table("case_data1.txt")
Y <- read.table("control_data1.txt")

resNulltest <- paca_null(X, Y, k, nperm = 100)

# p-value of rejecting the null
print(resNulltest$pval)
```
The input data, `X` & `Y` needs to be of form samples-by-features (NxM), where M >> N. The data also need to be standardized along the feature axis. Increasing `nperm` increases the precision of the `pval` estimate.
Please refer to the [PACA man page](man/ACA.Rd) for more detailed usage information.

## Troubleshooting

<details>
<summary>Instructions</summary>

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
</details>

## License and Disclaimer

PACA is publicly released under the GPL-3.0 license (full license text found [here](LICENSE.Rd)). Note however that the programs it calls may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running this script.

