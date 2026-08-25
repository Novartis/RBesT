# makefile written using https://yuukidach.github.io/p/makefile-for-projects-with-subdirectories/ as template 

TARGET = r-source

OUTDIR = ./build

# includes all src dirs excluding R/
SRCDIR = ./demo ./inst/stan ./inst/stan/include ./man-roxygen ./vignettes
##DIR_OBJ = ./obj

OUTDIR_ABS=$(abspath $(OUTDIR))
PROJROOT_ABS=$(abspath .)

RPKG=$(patsubst ‘%’, %, $(word 2, $(shell grep ^Package: DESCRIPTION)))
INCS = 
R_PKG_SRCS = $(wildcard R/*.R inst/examples/*R)
R_SRCS = $(wildcard *.R $(foreach fd, $(SRCDIR), $(fd)/*.R))
R_TEST_SRCS = $(wildcard tests/testthat/test*.R)
R_TEST_HELPER_SRCS = $(wildcard tests/testthat/helper*.R)
R_TEST_OBJS = $(R_TEST_SRCS:.R=.Rtest)
R_TESTFAST_OBJS = $(R_TEST_SRCS:.R=.Rtestfast)
FIXTURE_SRCS = $(wildcard tests/testthat/fixtures-mcmc-src/*_fixture.R)
FIXTURE_OBJS = $(patsubst tests/testthat/fixtures-mcmc-src/%_fixture.R,tests/testthat/fixtures-mcmc/%.rds,$(FIXTURE_SRCS))
COMPACT_FIXTURE_SRCS = $(wildcard tests/testthat/fixtures-compact/*_spec.R)
COMPACT_FIXTURE_RECIPE_SRCS = $(wildcard tests/testthat/fixtures-compact-src/*_fixture.R)
COMPACT_FIXTURE_OBJS = $(patsubst tests/testthat/fixtures-compact-src/%_fixture.R,tests/testthat/fixtures-compact/%_spec.R,$(COMPACT_FIXTURE_RECIPE_SRCS))
COMPACT_FIXTURE_REPORT_OBJS = $(patsubst tests/testthat/fixtures-compact-src/%_fixture.R,tests/testthat/fixtures-compact/%.report,$(COMPACT_FIXTURE_RECIPE_SRCS))
RMD_SRCS = $(wildcard *.Rmd $(foreach fd, $(SRCDIR), $(fd)/x*.Rmd))
STAN_SRCS = $(wildcard *.stan $(foreach fd, $(SRCDIR), $(fd)/*.stan))
SRCS = $(R_PKG_SRCS) $(R_SRCS) $(RMD_SRCS) $(STAN_SRCS)
NODIR_SRC = $(notdir $(SRCS))
BIN_OBJS = src/package-binary R/sysdata.rda
MANUAL_PDF = inst/doc/$(RPKG).pdf
DOC_OBJS = man/package-doc $(MANUAL_PDF)
# RCMD ?= R_PROFILE_USER="$(PROJROOT_ABS)/.Rprofile" "${R_HOME}/bin/R" -q
RCMD ?= "${R_HOME}/bin/R" -q
FIXTURE_FORCE ?= false

R_HOME ?= $(shell R RHOME)
PKG_VERSION ?= $(patsubst ‘%’, %, $(word 2, $(shell grep ^Version DESCRIPTION)))
GIT_TAG ?= v$(PKG_VERSION)

## Refresh the git export whenever HEAD moves. Wrapped in $(wildcard ...)
## because the branch ref may be packed away or HEAD may be detached.
GIT_DIR := $(shell git rev-parse --git-dir 2>/dev/null)
GIT_EXPORT_DEPS := $(wildcard $(GIT_DIR)/HEAD $(GIT_DIR)/packed-refs \
  $(GIT_DIR)/$(shell git symbolic-ref --quiet HEAD 2>/dev/null))
GIT_EXPORT := build/$(RPKG)-$(GIT_TAG).tar.gz

MD5 ?= md5sum
TMPDIR := $(realpath $(shell mktemp -d))

# When rendering vignettes/articles the recipes cd into the source
# directory, where R no longer picks up the repo-root .Renviron that puts
# the dev-installed RBesT on the library path. Point R_ENVIRON_USER at it,
# but only if the caller has not already set it and the file exists.
R_ENVIRON_PREFIX =
ifndef R_ENVIRON_USER
ifneq ($(wildcard $(CURDIR)/.Renviron),)
R_ENVIRON_PREFIX = R_ENVIRON_USER=$(CURDIR)/.Renviron
endif
endif

all : $(TARGET)

ifneq ($(filter true TRUE 1 yes YES,$(FIXTURE_FORCE)),)
FIXTURE_FORCE_PREREQ = FORCE
endif

.PHONY: FORCE
FORCE:

# tell makefile how to turn a Rmd into an md file
%.md : %.Rmd
	cd $(@D); echo running $(RCMD) -e "rmarkdown::render('$(<F)', output_format=rmarkdown::md_document(variant='markdown'))"
	cd $(@D); $(R_ENVIRON_PREFIX) $(RCMD) -e "rmarkdown::render('$(<F)', output_format=rmarkdown::md_document(variant='markdown'))"

%.md : %.R
	cd $(@D); echo running $(RCMD) -e "rmarkdown::render('$(<F)', output_format=rmarkdown::md_document(variant='markdown'))"
	cd $(@D); $(R_ENVIRON_PREFIX) $(RCMD) -q -e "rmarkdown::render('$(<F)', output_format=rmarkdown::md_document(variant='markdown'))"

# render an html via the respective md file
%.html : %.md
	cd $(@D); echo running $(RCMD) -e "rmarkdown::render('$(<F)', output_format=rmarkdown::html_document(self_contained=TRUE))"
	cd $(@D); $(R_ENVIRON_PREFIX) $(RCMD) -e "rmarkdown::render('$(<F)', output_format=rmarkdown::html_document(self_contained=TRUE))"

tests/testthat/fixtures-mcmc/%.rds : tests/testthat/fixtures-mcmc-src/%_fixture.R tools/build-test-fixture.R NAMESPACE $(BIN_OBJS) $(FIXTURE_FORCE_PREREQ)
	install -d $(@D)
	NOT_CRAN=true $(RCMD) --slave --file=tools/build-test-fixture.R --args $< $@

tests/testthat/fixtures-compact/%_spec.R : tests/testthat/fixtures-compact-src/%_fixture.R tools/build-compact-gmap-fixture.R tools/compact-gmap-fixture-utils.R $(R_TEST_HELPER_SRCS) NAMESPACE $(BIN_OBJS) $(FIXTURE_FORCE_PREREQ)
	install -d $(@D)
	NOT_CRAN=true $(RCMD) --slave --file=tools/build-compact-gmap-fixture.R --args $< $(@D) $*

tests/testthat/fixtures-compact/%.report : tests/testthat/fixtures-compact-src/%_fixture.R tests/testthat/fixtures-compact/%_spec.R tools/report-compact-gmap-fixtures.R tools/compact-gmap-fixture-utils.R $(R_TEST_HELPER_SRCS) NAMESPACE $(BIN_OBJS)
	install -d $(@D)
	@status=0; NOT_CRAN=true $(RCMD) --slave --file=tools/report-compact-gmap-fixtures.R --args $< $@ > $@.log 2>&1 || status=$$?; \
	cat $@.log; \
	if [ $$status -eq 0 ]; then rm -f $@.log; fi; \
	exit $$status

tests/%.Rtest : tests/%.R $(R_TEST_HELPER_SRCS) $(COMPACT_FIXTURE_SRCS) $(R_PKG_SRCS) NAMESPACE tools/run-test-file.R $(BIN_OBJS) $(FIXTURE_OBJS)
	@status=0; NOT_CRAN=true $(RCMD) --slave --file=tools/run-test-file.R --args $< > $@ 2>&1 || status=$$?; \
	printf "Test summary for $(<F): "; \
	grep '^\[' $@ | tail -n 1 || true; \
	exit $$status

# Fast/CRAN-like tests intentionally omit $(FIXTURE_OBJS); fixture-backed tests
# should skip cleanly when the local cache is unavailable.
tests/%.Rtestfast : tests/%.R $(R_TEST_HELPER_SRCS) $(COMPACT_FIXTURE_SRCS) $(R_PKG_SRCS) NAMESPACE tools/run-test-file.R $(BIN_OBJS)
	@status=0; NOT_CRAN=false $(RCMD) --slave --file=tools/run-test-file.R --args $< > $@ 2>&1 || status=$$?; \
	printf "Test summary for $(<F): "; \
	grep '^\[' $@ | tail -n 1 || true; \
	exit $$status


R/stanmodels.R: $(STAN_SRCS)
	## ensure that NAMESPACE contains load directive
	echo "# Generated by roxygen2: do not edit by hand" > NAMESPACE
	echo "import(Rcpp)" >> NAMESPACE
	echo "import(methods)" >> NAMESPACE
	echo "importFrom(rstan, sampling)" >> NAMESPACE
	echo "useDynLib($(RPKG), .registration = TRUE)" >> NAMESPACE
	install -d src
	"${R_HOME}/bin/Rscript" -e "rstantools::rstan_config()"
	touch R/stanmodels.R

src/package-binary: R/stanmodels.R
	## ensure that NAMESPACE contains load directive
	echo "# Generated by roxygen2: do not edit by hand" > NAMESPACE
	echo "import(Rcpp)" >> NAMESPACE
	echo "import(methods)" >> NAMESPACE
	echo "importFrom(rstan, sampling)" >> NAMESPACE
	echo "useDynLib($(RPKG), .registration = TRUE)" >> NAMESPACE
	install -d src
	"${R_HOME}/bin/Rscript" -e 'pkgbuild::compile_dll(debug=FALSE)'
	touch src/package-binary

man/package-doc: $(R_PKG_SRCS) $(BIN_OBJS)
	## NOTE: On a clean tree (after `make clean` removes man/*.Rd) roxygen2
	## 8.x cannot resolve intra-package [topic()] links because it reads the
	## on-disk man/*.Rd topic database, which does not exist yet. This emits
	## "Could not resolve link to topic ..." warnings for valid links such as
	## [gMAP()] or [mixfit()]. The warnings are benign: the generated Rd is
	## identical and a subsequent roxygenize (e.g. the next incremental build)
	## resolves all links with no warnings.
	"${R_HOME}/bin/Rscript" -e 'roxygen2::roxygenize()'
	touch man/package-doc

inst/sbc/sbc_report.html : inst/sbc/sbc_report.R inst/sbc/calibration.rds
	cd $(@D); echo running $(RCMD) -e "rmarkdown::render('$(<F)', output_format=rmarkdown::html_vignette(self_contained=TRUE))"
	cd $(@D); $(RCMD) -e "rmarkdown::render('$(<F)', output_format=rmarkdown::html_vignette(self_contained=TRUE))"


inst/sbc/calibration.rds :
	echo "Please run inst/sbc/make_reference_rankhist.R"
	exit 1

R/sysdata.rda: inst/sbc/calibration.rds
	"${R_HOME}/bin/R" --slave --file=tools/make-ds.R

# The reference manual cites the central bibliography via the Rdpack
# \insertRef{key}{RBesT} macro. At render time Rdpack resolves the key against
# the *installed* RBesT (system.file("REFERENCES.bib", package = "RBesT")), so
# the manual must be built against a dev-install that ships inst/REFERENCES.bib
# -- otherwise the \references sections render empty. Depend on the dev-install
# target and point R_LIBS_USER at it so \insertRef finds the bibliography.
inst/doc/$(RPKG).pdf : man/package-doc build/installed/$(RPKG)/DESCRIPTION
	install -d inst/doc
	R_LIBS_USER="$(CURDIR)/build/installed" "${R_HOME}/bin/R" CMD Rd2pdf --batch --no-preview --force --output=inst/doc/$(RPKG).pdf .
	"${R_HOME}/bin/R" --vanilla --slave -e 'library(tools); tools::compactPDF("inst/doc/$(RPKG).pdf")'


NAMESPACE: man/package-doc


# The canonical bibliography lives in inst/REFERENCES.bib (read by Rdpack for
# the Rd documentation via \insertRef). CRAN's "Writing R Extensions" requires
# the vignette BibTeX file to reside inside the vignette source directory, so
# generate vignettes/REFERENCES.bib from the canonical file. It is a derived
# artifact (git-ignored, like man/*.Rd and NAMESPACE): edit inst/REFERENCES.bib
# only and let the build keep the vignette copy in sync.
vignettes/REFERENCES.bib : inst/REFERENCES.bib
	@printf '%s\n' \
	  '% DO NOT EDIT -- this file is auto-generated from inst/REFERENCES.bib.' \
	  '% Edit inst/REFERENCES.bib (the canonical bibliography) instead and run' \
	  '% make to regenerate this vignette copy. Any manual changes here will be' \
	  '% overwritten by the build.' \
	  '' > $@
	cat $< >> $@


PHONY := $(TARGET)
$(TARGET): build/r-source-fast

## Standalone target so that the git export is one of the first artifacts to
## appear in build/, ahead of the expensive documentation and vignette
## prerequisites. Both source targets below consume it rather than rolling
## their own copy.
$(GIT_EXPORT) : $(GIT_EXPORT_DEPS)
	install -d build
	git archive --format=tar.gz --prefix $(RPKG)-$(GIT_TAG)/ HEAD > $@

build/r-source-fast : $(GIT_EXPORT) $(BIN_OBJS) man/package-doc $(SRCS) DESCRIPTION vignettes/REFERENCES.bib
	install -d build
	rm -rf build/$(RPKG)-$(GIT_TAG)
	cd build; tar x -C $(TMPDIR) -f $(RPKG)-$(GIT_TAG).tar.gz
	cp -v NAMESPACE $(TMPDIR)/$(RPKG)-$(GIT_TAG)
	## Overlay the working-tree DESCRIPTION so the built tarball reflects the
	## current (possibly uncommitted) version -- PKG_VERSION is derived from the
	## working-tree DESCRIPTION, so without this the tarball name built from the
	## git-archived HEAD DESCRIPTION would diverge and the mv below would fail.
	cp -v DESCRIPTION $(TMPDIR)/$(RPKG)-$(GIT_TAG)
	install -d $(TMPDIR)/$(RPKG)-$(GIT_TAG)/man
	cp -v man/*.Rd $(TMPDIR)/$(RPKG)-$(GIT_TAG)/man
	install -d $(TMPDIR)/$(RPKG)-$(GIT_TAG)/inst
	cp -v inst/REFERENCES.bib $(TMPDIR)/$(RPKG)-$(GIT_TAG)/inst
	install -d $(TMPDIR)/$(RPKG)-$(GIT_TAG)/vignettes
	cp -v vignettes/REFERENCES.bib $(TMPDIR)/$(RPKG)-$(GIT_TAG)/vignettes
	cd $(TMPDIR)/$(RPKG)-$(GIT_TAG); "${R_HOME}/bin/R" --slave --file=tools/make-ds.R
	cd $(TMPDIR); NOT_CRAN=false "${R_HOME}/bin/R" CMD build $(RPKG)-$(GIT_TAG) --no-build-vignettes --no-manual
	rm -rf $(TMPDIR)/$(RPKG)-$(GIT_TAG)
	mv $(TMPDIR)/$(RPKG)_$(PKG_VERSION).tar.gz build/$(RPKG)-source.tar.gz
	touch build/r-source-fast

build/r-source-release : $(GIT_EXPORT) $(BIN_OBJS) $(DOC_OBJS) $(SRCS) DESCRIPTION vignettes/REFERENCES.bib inst/sbc/sbc_report.html
	install -d build
	rm -rf build/$(RPKG)-$(GIT_TAG)
	## -P disables tar's delayed-symlink handling, which extracts relative
	## symlinks such as vignettes/articles/REFERENCES.bib via a mode 0000
	## placeholder file. That placeholder cannot be reopened on shared
	## (virtiofs/9p) working directories, making extraction fail there.
	cd build; tar -Pxf $(RPKG)-$(GIT_TAG).tar.gz
	cp -v NAMESPACE build/$(RPKG)-$(GIT_TAG)
	install -d build/$(RPKG)-$(GIT_TAG)/inst/doc
	cp -v inst/doc/$(RPKG).pdf build/$(RPKG)-$(GIT_TAG)/inst/doc
	cp -v inst/REFERENCES.bib build/$(RPKG)-$(GIT_TAG)/inst
	cp -v inst/sbc/sbc_report.html build/$(RPKG)-$(GIT_TAG)/inst/sbc/sbc_report.html
	install -d build/$(RPKG)-$(GIT_TAG)/vignettes
	cp -v vignettes/REFERENCES.bib build/$(RPKG)-$(GIT_TAG)/vignettes
	cd build/$(RPKG)-$(GIT_TAG); "${R_HOME}/bin/R" --slave --file=tools/make-ds.R
	install -d build/$(RPKG)-$(GIT_TAG)/man
	cp -v man/*.Rd build/$(RPKG)-$(GIT_TAG)/man
	# set NOT_CRAN=true to get vignettes render with full sampling
	cd build; NOT_CRAN=true $(RCMD) CMD build --compact-vignettes=both $(RPKG)-$(GIT_TAG)
	#cd build; NOT_CRAN=false "${R_HOME}/bin/R" CMD build $(RPKG)-$(GIT_TAG) --no-build-vignettes --no-manual
	rm -rf build/$(RPKG)-$(GIT_TAG)
	cd build; $(MD5) $(RPKG)-$(GIT_TAG).tar.gz > $(RPKG)-$(GIT_TAG).md5
	cd build; $(MD5) $(RPKG)_$(PKG_VERSION).tar.gz > $(RPKG)_$(PKG_VERSION).md5
	touch build/r-source-release

PHONY += r-source-release
r-source-release : build/r-source-release

PHONY += binary
binary : NAMESPACE src/package-binary

PHONY += derived
derived : NAMESPACE $(BIN_OBJS) $(DOC_OBJS)

PHONY += r-source-check
r-source-check : r-source
	cd build; tar xvzf $(RPKG)-source.tar.gz
	cd build; NOT_CRAN=true $(RCMD) CMD check $(RPKG)

PHONY += r-source-release-check
r-source-release-check : r-source-release
	cd build; tar xvzf $(RPKG)_$(PKG_VERSION).tar.gz
	cd build; NOT_CRAN=true $(RCMD) CMD check $(RPKG)

# Reverse dependency checks (CRAN readiness).
# REVDEP_WORKERS: reverse deps checked in parallel (each in its own subprocess).
# REVDEP_JOBS: C++ compile threads per package build (-jN). Keep
# REVDEP_WORKERS * REVDEP_JOBS near your core count and watch RAM: Stan/rstan
# compiles are memory-hungry (~1-2 GB each). Start at 6 x 1, tune from there.
REVDEP_WORKERS ?= 6
REVDEP_JOBS ?= 1

# Local run against the current working tree (includes uncommitted changes).
# Prerequisites make the source tree buildable/installable by revdepcheck:
# NAMESPACE + man/*.Rd (via NAMESPACE), compiled Stan src/ and R/sysdata.rda
# (via BIN_OBJS). The PDF manual/vignettes are intentionally omitted:
# revdepcheck builds with --no-manual --no-build-vignettes.
PHONY += revdepcheck
revdepcheck: NAMESPACE $(BIN_OBJS)
	$(RCMD) -e 'if (!requireNamespace("revdepcheck", quietly=TRUE)) pak::pak("r-lib/revdepcheck")'
	MAKEFLAGS="-j$(REVDEP_JOBS)" $(RCMD) -e 'revdepcheck::revdep_check(num_workers = $(REVDEP_WORKERS))'

# Run on GitHub's standard runners via workflow_dispatch. Uses the *pushed* tip
# of the current branch, NOT your local tree: commit and push before running.
PHONY += revdepcheck-ci
revdepcheck-ci:
	gh workflow run revdepcheck.yaml --ref "$$(git rev-parse --abbrev-ref HEAD)"
	sleep 5
	gh run watch --exit-status
	gh run download -n revdep-results || true

build/installed/$(RPKG)/DESCRIPTION : build/r-source-fast
	rm -rf build/installed
	install -d build/installed
	cd build; $(RCMD) CMD INSTALL --library=./installed --no-docs --no-multiarch --no-test-load --no-clean-on-error $(RPKG)-source.tar.gz

docs/index.html : doc $(SRCS) vignettes/REFERENCES.bib
	## Render the reference examples with the package default (full) sampling
	## so the embedded website plots and output look appropriate. The shipped
	## slim example template is temporarily replaced by the full-sampling
	## variant, the Rd files are regenerated, the site is built, and finally
	## the slim template and Rd are restored -- even if the build fails.
	cp -f man-roxygen/example-start.R man-roxygen/example-start.R.slim.bak; \
	cp -f man-roxygen/example-start-full.R man-roxygen/example-start.R; \
	"${R_HOME}/bin/Rscript" -e 'roxygen2::roxygenize(roclets="rd")'; \
	status=0; NOT_CRAN=true $(RCMD) -e 'pkgdown::build_site()' || status=$$?; \
	mv -f man-roxygen/example-start.R.slim.bak man-roxygen/example-start.R; \
	"${R_HOME}/bin/Rscript" -e 'roxygen2::roxygenize(roclets="rd")'; \
	exit $$status

PHONY += pkgdown
pkgdown: docs/index.html

PHONY += dev-install
dev-install: build/installed/$(RPKG)/DESCRIPTION

PHONY += test-all
test-all : $(R_TEST_OBJS)

PHONY += testfast-all
testfast-all : $(R_TESTFAST_OBJS)

PHONY += test-fixtures
test-fixtures : $(FIXTURE_OBJS)

PHONY += compact-fixtures
compact-fixtures : $(COMPACT_FIXTURE_OBJS)

PHONY += compact-fixture-report
compact-fixture-report : $(COMPACT_FIXTURE_REPORT_OBJS)
	@cat $(COMPACT_FIXTURE_REPORT_OBJS)

PHONY += clean-fixtures
# NOTE: the compact fixture *.report files are committed test-evidence and are
# intentionally NOT removed here (nor by `clean`, which depends on this target).
# Only `clean-all` removes them. The *.report.log files are transient build logs
# (git-ignored) and are safe to remove on every clean.
clean-fixtures:
	rm -f $(FIXTURE_OBJS)
	rm -f $(COMPACT_FIXTURE_REPORT_OBJS:%=%.log)

PHONY += clean-test-fixtures
clean-test-fixtures: clean-fixtures

PHONY += retestfast-all
retestfast-all : clean-test $(R_TESTFAST_OBJS)

PHONY += retest-all
retest-all : clean-test $(R_TEST_OBJS)

PHONY += check-winbuilder-devel
check-winbuilder-devel : r-source-release
	cd build; $(RCMD) -e 'target <- tempdir()' \
			  -e 'untar("$(RPKG)_$(PKG_VERSION).tar.gz", exdir=target)' \
			  -e 'devtools::check_win_devel(pkg=file.path(target, "$(RPKG)"))'

PHONY += check-winbuilder-release
check-winbuilder-release : r-source-release
	cd build; $(RCMD) -e 'target <- tempdir()' \
			  -e 'untar("$(RPKG)_$(PKG_VERSION).tar.gz", exdir=target)' \
			  -e 'devtools::check_win_release(pkg=file.path(target, "$(RPKG)"))'

PHONY += check-winbuilder-oldrelease
check-winbuilder-oldrelease : r-source-release
	cd build; $(RCMD) -e 'target <- tempdir()' \
			  -e 'untar("$(RPKG)_$(PKG_VERSION).tar.gz", exdir=target)' \
			  -e 'devtools::check_win_oldrelease(pkg=file.path(target, "$(RPKG)"))'

PHONY += check-winbuilder
check-winbuilder : check-winbuilder-devel check-winbuilder-release check-winbuilder-oldrelease

#$(DIR_OBJ)/%.o: %.c $(INCS)
#    mkdir -p $(@D)
#    $(CC) -o $@ $(CFLAGS) -c $< $(INC_DIRS)

PHONY += clean
clean: clean-fixtures
	rm -rf _brms-cache/*
	rm -rf build/*
	rm -f man/*.Rd
	rm -f NAMESPACE
	rm -f inst/doc/$(RPKG).pdf
	rm -f src/$(RPKG).so
	rm -f src/*.o
	rm -f man/package-doc
	rm -f src/package-binary
	rm -f R/sysdata.rda
	rm -f demo/*.html
	rm -f vignettes/*.html
	rm -f vignettes/*.docx
	rm -rf .Rd2pdf*
	rm -f $(R_TEST_OBJS)
	rm -f $(R_TESTFAST_OBJS)
	rm -rf src
	rm -f R/stanmodels.R

PHONY += clean-all
# Deep clean: everything `clean` removes, plus the committed compact fixture
# *.report test-evidence files. Kept separate so day-to-day `clean` preserves
# the reports.
clean-all: clean
	rm -f $(COMPACT_FIXTURE_REPORT_OBJS)

clean-test:
	rm -f $(R_TEST_OBJS)
	rm -f $(R_TESTFAST_OBJS)

PHONY += doc
doc: $(DOC_OBJS)

PHONY += echoes
echoes:
	@echo "INC files: $(INCS)"
	@echo "SRC files: $(SRCS)"
	@echo "OBJ files: $(OBJS)"

PHONY += help
help:
	@echo "RBesT package development targets"
	@echo "=================================="
	@echo ""
	@echo "Build & install:"
	@echo "  r-source              Build source package (fast, no vignettes)"
	@echo "  r-source-release      Build release source package (with vignettes)"
	@echo "  binary                Compile Stan models and shared library"
	@echo "  derived               Generate NAMESPACE, binary, and docs"
	@echo "  dev-install           Install from source into build/installed/"
	@echo "  doc                   Generate Rd documentation"
	@echo ""
	@echo "Testing:"
	@echo "  testfast-all          Run all tests in fast/CRAN-like mode"
	@echo "  test-all              Run all tests with full fixtures"
	@echo "  retestfast-all        Clean test output and re-run fast tests"
	@echo "  retest-all            Clean test output and re-run full tests"
	@echo "  tests/testthat/test-FOO.Rtestfast   Run single test file (fast)"
	@echo "  tests/testthat/test-FOO.Rtest       Run single test file (full)"
	@echo ""
	@echo "Fixtures:"
	@echo "  test-fixtures         Build MCMC fixture .rds files"
	@echo "  compact-fixtures      Build compact fixture specs"
	@echo "  compact-fixture-report  Report compact fixture quality"
	@echo "  clean-fixtures        Remove regenerable fixture outputs (keeps *.report)"
	@echo ""
	@echo "Checks:"
	@echo "  r-source-check        R CMD check on source package"
	@echo "  r-source-release-check  R CMD check on release package"
	@echo "  revdepcheck           Reverse dependency check locally (REVDEP_WORKERS, REVDEP_JOBS)"
	@echo "  revdepcheck-ci        Reverse dependency check on GitHub runners (pushed tip of current branch)"
	@echo "  check-winbuilder      Submit to winbuilder (devel+release+old)"
	@echo ""
	@echo "Documentation:"
	@echo "  pkgdown               Build pkgdown site"
	@echo ""
	@echo "Housekeeping:"
	@echo "  clean                 Remove all generated artifacts (keeps *.report)"
	@echo "  clean-all             Remove all generated artifacts and *.report evidence"
	@echo "  clean-test            Remove test output files only"
	@echo ""
	@echo "Variables:"
	@echo "  FIXTURE_FORCE=true    Force rebuild of fixtures"
	@echo "  print-VARNAME         Print value of any Makefile variable"

##
# Debug target that allows you to print a variable
##
print-%  : ; @echo $* = $($*)


.PHONY : $(PHONY)
