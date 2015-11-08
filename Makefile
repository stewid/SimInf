# Determine package name and version from DESCRIPTION file
PKG_VERSION=$(shell grep -i ^version DESCRIPTION | cut -d : -d \  -f 2)
PKG_NAME=$(shell grep -i ^package DESCRIPTION | cut -d : -d \  -f 2)

# Roxygen version to check before generating documentation
ROXYGEN_VERSION=5.0.0

# Name of built package
PKG_TAR=$(PKG_NAME)_$(PKG_VERSION).tar.gz

# Install package
install:
	cd .. && R CMD INSTALL $(PKG_NAME)

# Build documentation with roxygen
# 1) Check version of roxygen2 before building documentation
# 2) Remove old doc
# 3) Generate documentation
roxygen:
	Rscript -e "library(roxygen2); stopifnot(packageVersion('roxygen2') == '$(ROXYGEN_VERSION)')"
	rm -f man/*.Rd
	cd .. && Rscript -e "library(roxygen2); roxygenize('$(PKG_NAME)')"

# Generate PDF output from the Rd sources
# 1) Rebuild documentation with roxygen
# 2) Generate pdf, overwrites output file if it exists
pdf: roxygen
	cd .. && R CMD Rd2pdf --force $(PKG_NAME)

# Generate README
README.md: README.Rmd
	Rscript -e "library(knitr); knit('README.Rmd')"

# Build and check package
check: clean
	cd .. && R CMD build --no-build-vignettes $(PKG_NAME)
	cd .. && _R_CHECK_CRAN_INCOMING_=FALSE R CMD check --as-cran \
	--no-manual --no-vignettes --no-build-vignettes $(PKG_TAR)

# Build and check package with gctorture
check_gctorture: clean
	cd .. && R CMD build --no-build-vignettes $(PKG_NAME)
	cd .. && R CMD check --no-manual --no-vignettes --no-build-vignettes --use-gct $(PKG_TAR)

# Build and check package with valgrind
check_valgrind: clean
	cd .. && R CMD build --no-build-vignettes $(PKG_NAME)
	cd .. && _R_CHECK_CRAN_INCOMING_=FALSE R CMD check --as-cran \
	--no-manual --no-vignettes --no-build-vignettes --use-valgrind $(PKG_TAR)

# Run all tests with valgrind
test_objects = $(wildcard tests/*.R)
valgrind:
	$(foreach var,$(test_objects),R -d "valgrind --tool=memcheck --leak-check=full" --vanilla < $(var);)

configure: configure.ac
	autoconf ./configure.ac > ./configure
	chmod +x ./configure

clean:
	-rm -f config.log
	-rm -f config.status
	-rm -rf autom4te.cache
	-rm -f local320.zip
	-rm -f gsl.zip
	-rm -f src/Makevars
	-rm -f src/*.o
	-rm -f src/core/*.o
	-rm -f src/events/*.o
	-rm -f src/misc/*.o
	-rm -f src/models/*.o
	-rm -f src/*.so
	-rm -rf src-x64
	-rm -rf src-i386
	-rm -rf gsl
	-rm -rf src/gsl

.PHONY: install roxygen pdf check check_gctorture check_valgrind clean
