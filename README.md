[![Build Status](https://travis-ci.org/stewid/siminf.svg)](https://travis-ci.org/stewid/siminf)
[![Build status](https://ci.appveyor.com/api/projects/status/pe68xiu1anxvet2n?svg=true)](https://ci.appveyor.com/project/stewid/siminf)
[![Coverage Status](https://coveralls.io/repos/stewid/siminf/badge.svg?branch=master)](https://coveralls.io/r/stewid/siminf?branch=master)

# siminf

A flexible and efficient framework for stochastic disease spread
simulations in animal populations.

## Installation

To install the development version of `siminf`, it's easiest to use
the devtools package:

```r
# install.packages("devtools")
library(devtools)
install_github("stewid/siminf"")
```

Another alternative is to use `git` and `make`

```r
$ git clone https://github.com/stewid/siminf.git
$ cd siminf
$ make install
```

### Dependencies

#### On a Debian based system

Before installing `siminf` you need to install
[The GNU Scientific Library](http://www.gnu.org/software/gsl/) (GSL).

```
apt-get install libgsl0-dev
```

You may need to do:
```
sudo apt-get install libgsl0-dev
```

#### On a Windows machine

* Download and install
  [Rtools](http://cran.r-project.org/bin/windows/Rtools/)

# Authors

* Pavol Bauer (Uppsala University, Sweden)
* Stefan Engblom (Uppsala University, Sweden)
* Stefan Widgren (National Veterinary Institute, Sweden) **(Maintainer)**

Any suggestions, bug reports, forks and pull requests are
appreciated. Get in touch.

# License

The `siminf` package is licensed under the GPLv3.

- [LICENSE](LICENSE) - `siminf` package license
