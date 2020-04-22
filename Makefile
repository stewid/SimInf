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
	cd .. && Rscript -e "roxygen2::roxygenize('$(PKG_NAME)')"

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
	cd .. && OMP_THREAD_LIMIT=2 _R_CHECK_CRAN_INCOMING_=FALSE R CMD check \
        --no-stop-on-test-error --as-cran --run-dontrun $(PKG_TAR)

# Check package (without manual and vignettes)
check_quick: clean
	cd .. && R CMD build --no-build-vignettes --no-manual $(PKG_NAME)
	cd .. && OMP_THREAD_LIMIT=2 _R_CHECK_CRAN_INCOMING_=FALSE R CMD check \
        --no-stop-on-test-error --no-vignettes --no-manual --as-cran $(PKG_TAR)

# Build and check package with gctorture
check_gctorture:
	cd .. && R CMD build --no-build-vignettes $(PKG_NAME)
	cd .. && R CMD check --no-manual --no-vignettes --no-build-vignettes --use-gct $(PKG_TAR)

# Build and check package with valgrind
check_valgrind:
	cd .. && R CMD build --no-build-vignettes $(PKG_NAME)
	cd .. && _R_CHECK_CRAN_INCOMING_=FALSE R CMD check --as-cran \
	--no-manual --no-vignettes --no-build-vignettes --use-valgrind $(PKG_TAR)

# Check to create a package with 'package_skeleton' and then run 'R
# CMD check' on that package.
check_pkg_skeleton:
	cd .. && rm -rf pkggdata
	cd .. && Rscript \
            -e "library('SimInf')" \
            -e "model <- mparse(transitions = c('S->b*S*I->I', 'I->g*I->R')," \
            -e "    compartments = c('S', 'I', 'R')," \
            -e "    gdata = c(b = 0.16, g = 0.077)," \
            -e "    u0 = data.frame(S = 100, I = 1, R = 0)," \
            -e "    tspan = 1:100)" \
            -e "package_skeleton(model = model, name = 'pkggdata')" \
        && R CMD build pkggdata \
        && R CMD check pkggdata_1.0.tar.gz \
        && R CMD INSTALL pkggdata_1.0.tar.gz
	cd .. && rm -rf pkgldata
	cd .. && Rscript \
            -e "library('SimInf')" \
            -e "model <- mparse(transitions = c('S->b*S*I->I', 'I->g*I->R')," \
            -e "    compartments = c('S', 'I', 'R')," \
            -e "    ldata = data.frame(b = 0.16, g = 0.077)," \
            -e "    u0 = data.frame(S = 100, I = 1, R = 0)," \
            -e "    tspan = 1:100)" \
            -e "package_skeleton(model = model, name = 'pkgldata')" \
        && R CMD build pkgldata \
        && R CMD check pkgldata_1.0.tar.gz \
        && R CMD INSTALL pkgldata_1.0.tar.gz
	cd .. && rm -rf pkgv0
	cd .. && Rscript \
            -e "library('SimInf')" \
            -e "model <- mparse(transitions = c('S->b*S*I->I', 'I->g*I->R')," \
            -e "    compartments = c('S', 'I', 'R')," \
            -e "    v0 = data.frame(b = 0.16, g = 0.077)," \
            -e "    u0 = data.frame(S = 100, I = 1, R = 0)," \
            -e "    tspan = 1:100)" \
            -e "package_skeleton(model = model, name = 'pkgv0')" \
        && R CMD build pkgv0 \
        && R CMD check pkgv0_1.0.tar.gz \
        && R CMD INSTALL pkgv0_1.0.tar.gz
	cd .. && Rscript \
            -e "library('pkggdata')" \
            -e "library('pkgldata')" \
            -e "library('pkgv0')" \
            -e "u0 <- data.frame(S = 100, I = 1, R = 0)" \
            -e "model_gdata <- pkggdata(u0 = u0, gdata = c(b = 0.16, g = 0.077), tspan = 1:100)" \
            -e "set.seed(22)" \
            -e "result_gdata <- prevalence(run(model_gdata), I ~ .)" \
            -e "model_ldata <- pkgldata(u0 = u0, ldata = data.frame(b = 0.16, g = 0.077), tspan = 1:100)" \
            -e "set.seed(22)" \
            -e "result_ldata <- prevalence(run(model_ldata), I ~ .)" \
            -e "model_v0 <- pkgv0(u0 = u0, v0 = data.frame(b = 0.16, g = 0.077), tspan = 1:100)" \
            -e "set.seed(22)" \
            -e "result_v0 <- prevalence(run(model_v0), I ~ .)" \
            -e "stopifnot(identical(result_gdata, result_ldata))" \
            -e "stopifnot(identical(result_gdata, result_v0))"
	R CMD REMOVE pkggdata pkgldata pkgv0

# Build and check package on R-hub
rhub: clean check
	cd .. && Rscript -e "rhub::check(path='$(PKG_TAR)', rhub::platforms()[['name']], show_status = FALSE)"

# Build and use 'rchk' on package on R-hub
rchk: clean check
	cd .. && Rscript -e "rhub::check(path='$(PKG_TAR)', 'ubuntu-rchk', show_status = FALSE)"

# Run all tests with valgrind
test_objects = $(wildcard tests/*.R)
valgrind:
	$(foreach var,$(test_objects),R -d "valgrind --tool=memcheck --leak-check=full" --vanilla < $(var);)

# Run static code analysis on the C code.
# https://github.com/danmar/cppcheck/
.PHONY: cppcheck
cppcheck:
	cppcheck src

configure: configure.ac
	autoconf ./configure.ac > ./configure
	chmod +x ./configure

clean:
	./cleanup

.PHONY: install roxygen pdf build check check_quick check_gctorture check_valgrind check_pkg_skeleton rhub clean vignette
