

[![Build Status](https://travis-ci.org/stewid/SimInf.svg)](https://travis-ci.org/stewid/SimInf)
[![Build status](https://ci.appveyor.com/api/projects/status/pe68xiu1anxvet2n?svg=true)](https://ci.appveyor.com/project/stewid/SimInf)
[![CRAN status](http://www.r-pkg.org/badges/version/SimInf)](http://cran.r-project.org/web/packages/SimInf/index.html)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/last-month/SimInf)](http://cran.r-project.org/web/packages/SimInf/index.html)
[![Coverage Status](https://coveralls.io/repos/stewid/SimInf/badge.svg?branch=master&service=github)](https://coveralls.io/github/stewid/SimInf?branch=master)

# SimInf

A flexible and efficient framework for stochastic disease spread
simulations in animal populations.

## Installation

To install the development version of `SimInf`, it's easiest to use
the devtools package:


```r
# install.packages("devtools")
library(devtools)
install_github("stewid/SimInf"")
```

Another alternative is to use `git` and `make`


```r
$ git clone https://github.com/stewid/SimInf.git
$ cd SimInf
$ make install
```

### Dependencies

#### On a Debian based system

Before installing `SimInf` you need to install
[The GNU Scientific Library](http://www.gnu.org/software/gsl/) (GSL).


```r
apt-get install libgsl0-dev
```

You may need to do:

```r
sudo apt-get install libgsl0-dev
```

#### On a Windows machine

* Download and install
  [Rtools](http://cran.r-project.org/bin/windows/Rtools/)

# Authors

* Pavol Bauer (Uppsala University, Sweden)
* Stefan Engblom (Uppsala University, Sweden)
* Stefan Widgren (National Veterinary Institute, Sweden) **(Maintainer)**

# Acknowledgment

This work was financially supported by the Swedish Research Council
within the UPMARC Linnaeus center of Excellence (Pavol Bauer and
Stefan Engblom).

# License

The `SimInf` package is licensed under the GPLv3.

- [LICENSE](LICENSE) - `SimInf` package license

Any suggestions, bug reports, forks and pull requests are
appreciated. Get in touch.
