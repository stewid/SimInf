# Determine package name and version from DESCRIPTION file
PKG_VERSION=$(shell grep -i ^version DESCRIPTION | cut -d : -d \  -f 2)
PKG_NAME=$(shell grep -i ^package DESCRIPTION | cut -d : -d \  -f 2)

# Name of built package
PKG_TAR=$(PKG_NAME)_$(PKG_VERSION).tar.gz

# Install package
install:
	cd .. && R CMD INSTALL $(PKG_NAME)

# Build documentation with roxygen
# 1) Remove old doc
# 2) Generate documentation
roxygen:
	rm -f man/*.Rd
	Rscript -e "library(methods)" \
                -e "library(devtools)" \
                -e "devtools::document()"

# Generate PDF output from the Rd sources
# 1) Rebuild documentation with roxygen
# 2) Generate pdf, overwrites output file if it exists
pdf: roxygen
	cd .. && R CMD Rd2pdf --force $(PKG_NAME)

# Generate README
README.md: README.Rmd
	Rscript -e "library(knitr); knit('README.Rmd')"

# Generate vignette
vignette:
	cd vignettes && Rscript -e "library('methods')" \
                                -e "Sweave('SimInf')" \
                                -e "tools::texi2pdf('SimInf.tex')"

# Build package
build: clean
	cd .. && R CMD build --compact-vignettes=both $(PKG_NAME)

# Check package
check: build
	cd .. && _R_CHECK_CRAN_INCOMING_=FALSE R CMD check \
        --no-stop-on-test-error --as-cran $(PKG_TAR)

# Check package (without manual and vignettes)
check_quick: clean
	cd .. && R CMD build --no-build-vignettes --no-manual $(PKG_NAME)
	cd .. && _R_CHECK_CRAN_INCOMING_=FALSE R CMD check --no-stop-on-test-error \
        --no-vignettes --no-manual --as-cran $(PKG_TAR)

# Build and check package with gctorture
check_gctorture:
	cd .. && R CMD build --no-build-vignettes $(PKG_NAME)
	cd .. && R CMD check --no-manual --no-vignettes --no-build-vignettes --use-gct $(PKG_TAR)

# Build and check package with valgrind
check_valgrind:
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
	./cleanup

.PHONY: install roxygen pdf build check check_quick check_gctorture check_valgrind clean vignette
