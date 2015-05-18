[![Build Status](https://travis-ci.org/stewid/siminf.svg)](https://travis-ci.org/stewid/siminf)
[![Build status](https://ci.appveyor.com/api/projects/status/pe68xiu1anxvet2n?svg=true)](https://ci.appveyor.com/project/stewid/siminf)

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

Before installing `siminf` you need to install
[The GNU Scientific Library](http://www.gnu.org/software/gsl/) (GSL).

#### On a Debian based system

```
apt-get install libgsl0-dev
```

#### On a Windows machine

* Download and install
  [Rtools](http://cran.r-project.org/bin/windows/Rtools/)

* Download `local320.zip` from
  http://www.stats.ox.ac.uk/pub/Rtools/libs.html and unpack to a
  folder

* Change `LOCAL_SOFT` in `c:/Program/R/R-3.1.2/etc/i386/Makeconf` and
  `c:/Program/R/R-3.1.2/etc/x64/Makeconf` to the folder of the
  unpacked `local320.zip`

# Authors

* Pavol Bauer (Uppsala University, Sweden)
* Stefan Engblom (Uppsala University, Sweden)
* Stefan Widgren (National Veterinary Institute, Sweden) **(Maintainer)**

Any suggestions, bug reports, forks and pull requests are
appreciated. Get in touch.

# License

The `siminf` package is licensed under the GPLv3.

- [LICENSE](LICENSE) - `siminf` package license
