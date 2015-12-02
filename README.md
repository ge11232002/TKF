# TKF

[![Travis-CI Build Status](https://travis-ci.org/ge11232002/TKF.svg?branch=master)](https://travis-ci.org/ge11232002/TKF)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/TKF)](http://cran.r-project.org/package=TKF)

Pairwise Distance Estimation with TKF91 and TKF92 Model

## Installation of the stable version of `TKF` from CRAN

```R
install.packages("TKF")
```

## Installation of the development version of `TKF` from github
1. Make sure you have a working development environment.
    * **Windows**: Install [Rtools](http://cran.r-project.org/bin/windows/Rtools/).
    * **Mac**: Install Xcode from the Mac App Store.
    * **Linux**: Install a compiler and various development libraries (details vary across different flavors of Linux).
2. Install the GNU gsl.
    * **Windows**: It comes with Rtools.
    * **Mac**: `brew install gsl`
    * **Linux**: In Ubuntu, `apt-get install libgsl0-dev`
3. Install within R

* **Mac and Linux**:

  ```R
  devtools::install_github("ge11232002/TKF")
  ```
