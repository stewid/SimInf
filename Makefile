# Determine package name and version from DESCRIPTION file
PKG_VERSION=$(shell grep -i ^version DESCRIPTION | cut -d : -d \  -f 2)
PKG_NAME=$(shell grep -i ^package DESCRIPTION | cut -d : -d \  -f 2)

# Name of built package
PKG_TAR=$(PKG_NAME)_$(PKG_VERSION).tar.gz

# Install package
.PHONY: install
install:
	cd .. && R CMD INSTALL $(PKG_NAME)

# Build documentation with roxygen (first delete previous roxygen files)
.PHONY: roxygen
roxygen:
	Rscript -e "roxygen2::roxygenize(clean = TRUE)"

# Generate PDF output from the Rd sources
# 1) Rebuild documentation with roxygen
# 2) Generate pdf, overwrites output file if it exists
.PHONY: pdf
pdf: roxygen
	cd .. && R CMD Rd2pdf --force $(PKG_NAME)

# Generate README
README.md: README.Rmd
	Rscript -e "library(knitr); knit('README.Rmd')"

# Generate vignette
.PHONY: vignette
vignette:
	cd vignettes && Rscript -e "library('methods')" \
                                -e "Sweave('SimInf')" \
                                -e "tools::texi2pdf('SimInf.tex')"

# Build package
.PHONY: build
build: clean
	cd .. && R CMD build --compact-vignettes=both $(PKG_NAME)

# Check package
.PHONY: check
check: build
	cd .. && OMP_THREAD_LIMIT=2 _R_CHECK_CRAN_INCOMING_=FALSE R CMD check \
        --no-stop-on-test-error --as-cran --run-dontrun $(PKG_TAR)

# Check package (without manual and vignettes)
.PHONY: check_quick
check_quick: clean
	cd .. && R CMD build --no-build-vignettes --no-manual $(PKG_NAME)
	cd .. && \
        OMP_THREAD_LIMIT=2 \
        _R_CHECK_CRAN_INCOMING_=FALSE \
        _R_CHECK_SYSTEM_CLOCK_=0 \
        R CMD check \
        --no-stop-on-test-error --ignore-vignettes --no-manual --as-cran $(PKG_TAR)

# Build and check package with gctorture
.PHONY: check_gctorture
check_gctorture:
	cd .. && R CMD build --no-build-vignettes $(PKG_NAME)
	cd .. && R CMD check --no-manual --no-vignettes --no-build-vignettes --use-gct $(PKG_TAR)

# Build and check package with valgrind
.PHONY: check_valgrind
check_valgrind:
	cd .. && R CMD build --no-build-vignettes $(PKG_NAME)
	cd .. && _R_CHECK_CRAN_INCOMING_=FALSE R CMD check --as-cran \
	--no-manual --no-vignettes --no-build-vignettes --use-valgrind $(PKG_TAR)

# Check to create a package with 'package_skeleton' and then run 'R
# CMD check' on that package.
.PHONY: check_pkg_skeleton
check_pkg_skeleton:
	cd .. && rm -rf pkg.gdata
	cd .. && Rscript \
            -e "library('SimInf')" \
            -e "model <- mparse(transitions = c('S->b*S*I->I', 'I->g*I->R')," \
            -e "    compartments = c('S', 'I', 'R')," \
            -e "    gdata = c(b = 0.16, g = 0.077)," \
            -e "    u0 = data.frame(S = 100, I = 1, R = 0)," \
            -e "    tspan = 1:100)" \
            -e "package_skeleton(model = model, name = 'pkg.gdata')" \
        && R CMD build pkg.gdata \
        && R CMD check pkg.gdata_1.0.tar.gz \
        && R CMD INSTALL pkg.gdata_1.0.tar.gz
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
            -e "library('pkg.gdata')" \
            -e "library('pkgldata')" \
            -e "library('pkgv0')" \
            -e "u0 <- data.frame(S = 100, I = 1, R = 0)" \
            -e "model_gdata <- pkg.gdata(u0 = u0, gdata = c(b = 0.16, g = 0.077), tspan = 1:100)" \
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
	R CMD REMOVE pkg.gdata pkgldata pkgv0

# Check reverse dependencies
#
# 1) Install packages (in ../revdep/lib) to check the reverse dependencies.
# 2) Check the reverse dependencies using 'R CMD check'.
# 3) Collect results from the '00check.log' files.
.PHONY: revdep
revdep: revdep_install revdep_check revdep_results

# Install packages to check reverse dependencies
.PHONY: revdep_install
revdep_install: clean
	mkdir -p ../revdep/lib
	cd .. && R CMD INSTALL --library=revdep/lib $(PKG_NAME)
	R_LIBS_USER=../revdep/lib Rscript --vanilla \
          -e "options(repos = c(CRAN='https://cran.r-project.org'))" \
          -e "pkg <- tools::package_dependencies('$(PKG_NAME)', which = 'all', reverse = TRUE)" \
          -e "pkg <- as.character(unlist(pkg))" \
          -e "dep <- sapply(pkg, tools::package_dependencies, which = 'all')" \
          -e "dep <- as.character(unlist(dep))" \
          -e "if ('BiocInstaller' %in% dep) {" \
          -e "    source('https://bioconductor.org/biocLite.R')" \
          -e "    biocLite('BiocInstaller')" \
          -e "}" \
          -e "install.packages(pkg, dependencies = TRUE)" \
          -e "download.packages(pkg, destdir = '../revdep')"

# Check reverse dependencies with 'R CMD check'
.PHONY: revdep_check
revdep_check:
	$(foreach var,$(wildcard ../revdep/*.tar.gz),R_LIBS_USER=../revdep/lib \
          _R_CHECK_CRAN_INCOMING_=FALSE R --vanilla CMD check --as-cran \
          --no-stop-on-test-error --output=../revdep $(var) \
          | tee --append ../revdep/00revdep.log;)

# Collect results from checking reverse dependencies
.PHONY: revdep_results
revdep_results:
	Rscript --vanilla \
          -e "options(repos = c(CRAN='https://cran.r-project.org'))" \
          -e "pkg <- tools::package_dependencies('$(PKG_NAME)', which = 'all', reverse = TRUE)" \
          -e "pkg <- as.character(unlist(pkg))" \
          -e "results <- do.call('rbind', lapply(pkg, function(x) {" \
          -e "    filename <- paste0('../revdep/', x, '.Rcheck/00check.log')" \
          -e "    if (file.exists(filename)) {" \
          -e "        lines <- readLines(filename)" \
          -e "        status <- sub('^Status: ', '', lines[grep('^Status: ', lines)])" \
          -e "    } else {" \
          -e "        status <- 'missing'" \
          -e "    }" \
          -e "    data.frame(Package = x, Status = status)" \
          -e "}))" \
          -e "results <- results[order(results[, 'Status']), ]" \
          -e "rownames(results) <- NULL" \
          -e "cat('\n\n*** Results ***\n\n')" \
          -e "results" \
          -e "cat('\n\n')"

# Build and check package on R-hub
.PHONY: rhub
rhub: clean check
	cd .. && Rscript -e "rhub::check(path='$(PKG_TAR)', rhub::platforms()[['name']], show_status = FALSE)"

# Build and use 'rchk' on package on R-hub
.PHONY: rchk
rchk: clean check
	cd .. && Rscript -e "rhub::check(path='$(PKG_TAR)', 'ubuntu-rchk', show_status = FALSE)"

# Build and check package on https://win-builder.r-project.org/
.PHONY: winbuilder
winbuilder: clean check
	cd .. && curl -T $(PKG_TAR) ftp://win-builder.r-project.org/R-oldrelease/
	cd .. && curl -T $(PKG_TAR) ftp://win-builder.r-project.org/R-release/
	cd .. && curl -T $(PKG_TAR) ftp://win-builder.r-project.org/R-devel/

# Run covr
.PHONY: covr
covr:
	Rscript -e "covr::report(file = 'covr.html')"
	xdg-open covr.html

# Run static code analysis
.PHONY: lintr
lintr:
	Rscript -e "lintr::lint_package()"

# Run all tests with valgrind
test_objects = $(wildcard tests/*.R)
.PHONY: valgrind
valgrind:
	$(foreach var,$(test_objects),R -d "valgrind --tool=memcheck --leak-check=full" --vanilla < $(var);)

# Run static code analysis on the C code.
# https://github.com/danmar/cppcheck/
.PHONY: cppcheck
cppcheck:
	cppcheck -I ./inst/include/ src

configure: configure.ac
	autoconf ./configure.ac > ./configure
	chmod +x ./configure

.PHONY: clean
clean:
	./cleanup
	-rm -rf ../revdep
