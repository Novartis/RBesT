# makefile written using https://yuukidach.github.io/p/makefile-for-projects-with-subdirectories/ as template 

TARGET = r-source

OUTDIR = ./build

# includes all src dirs excluding R/
SRCDIR = ./demo ./inst/stan ./inst/stan/include ./man-roxygen
##DIR_OBJ = ./obj

OUTDIR_ABS=$(abspath $(OUTDIR))
PROJROOT_ABS=$(abspath .)

RPKG=$(patsubst ‘%’, %, $(word 2, $(shell grep ^Package: DESCRIPTION)))
INCS = 
R_PKG_SRCS = $(wildcard R/*.R inst/examples/*R)
R_SRCS = $(wildcard *.R $(foreach fd, $(SRCDIR), $(fd)/*.R))
R_TEST_SRCS = $(wildcard tests/testthat/test*.R)
R_TEST_OBJS = $(R_TEST_SRCS:.R=.Rtest)
R_TESTFAST_OBJS = $(R_TEST_SRCS:.R=.Rtestfast)
RMD_SRCS = $(wildcard *.Rmd $(foreach fd, $(SRCDIR), $(fd)/x*.Rmd))
STAN_SRCS = $(wildcard *.stan $(foreach fd, $(SRCDIR), $(fd)/*.stan))
SRCS = $(R_PKG_SRCS) $(R_SRCS) $(RMD_SRCS) $(STAN_SRCS)
NODIR_SRC = $(notdir $(SRCS))
BIN_OBJS = src/package-binary R/sysdata.rda
DOC_OBJS = man/package-doc inst/doc/$(RPKG).pdf
RCMD ?= R_PROFILE_USER="$(PROJROOT_ABS)/.Rprofile" "${R_HOME}/bin/R" -q

R_HOME ?= $(shell R RHOME)
PKG_VERSION ?= $(patsubst ‘%’, %, $(word 2, $(shell grep ^Version DESCRIPTION)))
GIT_TAG ?= v$(PKG_VERSION)

MD5 ?= md5sum
TMPDIR := $(realpath $(shell mktemp -d))

all : $(TARGET)

# tell makefile how to turn a Rmd into an md file
%.md : %.Rmd
	cd $(@D); echo running $(RCMD) -e "rmarkdown::render('$(<F)', output_format=rmarkdown::md_document(variant='markdown'))"
	cd $(@D); $(RCMD) -e "rmarkdown::render('$(<F)', output_format=rmarkdown::md_document(variant='markdown'))"

%.md : %.R
	cd $(@D); echo running $(RCMD) -e "rmarkdown::render('$(<F)', output_format=rmarkdown::md_document(variant='markdown'))"
	cd $(@D); $(RCMD) -q -e "rmarkdown::render('$(<F)', output_format=rmarkdown::md_document(variant='markdown'))"

# render an html via the respective md file
%.html : %.md
	cd $(@D); echo running $(RCMD) -e "rmarkdown::render('$(<F)', output_format=rmarkdown::html_document(self_contained=TRUE))"
	cd $(@D); $(RCMD) -e "rmarkdown::render('$(<F)', output_format=rmarkdown::html_document(self_contained=TRUE))"

tests/%.Rtest : tests/%.R
	NOT_CRAN=true $(RCMD) -e "devtools::load_all()" -e "test_file('$<')" > $@ 2>&1
	@printf "Test summary for $(<F): "
	@grep '^\[' $@ | tail -n 1

tests/%.Rtestfast : tests/%.R
	NOT_CRAN=false $(RCMD) -e "devtools::load_all()" -e "test_file('$<')" > $@ 2>&1
	@printf "Test summary for $(<F): "
	@grep '^\[' $@ | tail -n 1


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
	"${R_HOME}/bin/Rscript" -e 'roxygen2::roxygenize()'
	touch man/package-doc

inst/sbc/sbc_report.html : inst/sbc/calibration.rds

inst/sbc/calibration.rds :
	echo "Please run inst/sbc/make_reference_rankhist.R"
	exit 1

R/sysdata.rda: inst/sbc/calibration.rds
	"${R_HOME}/bin/R" --slave --file=tools/make-ds.R

inst/doc/$(RPKG).pdf : man/package-doc
	install -d inst/doc
	"${R_HOME}/bin/R" CMD Rd2pdf --batch --no-preview --force --output=inst/doc/$(RPKG).pdf .
	"${R_HOME}/bin/R" --vanilla --slave -e 'library(tools); tools::compactPDF("inst/doc/$(RPKG).pdf")'


NAMESPACE: man/package-doc


PHONY := $(TARGET)
$(TARGET): build/r-source-fast

build/r-source-fast : $(BIN_OBJS) $(DOC_OBJS) $(SRCS)
	install -d build
	git archive --format=tar.gz --prefix $(RPKG)-$(GIT_TAG)/ HEAD > build/$(RPKG)-$(GIT_TAG).tar.gz
	rm -rf build/$(RPKG)-$(GIT_TAG)
	cd build; tar x -C $(TMPDIR) -f $(RPKG)-$(GIT_TAG).tar.gz
	rm -f build/$(RPKG)-$(GIT_TAG).tar.gz
	cp -v NAMESPACE $(TMPDIR)/$(RPKG)-$(GIT_TAG)
	install -d $(TMPDIR)/$(RPKG)-$(GIT_TAG)/man
	cp -v man/*.Rd $(TMPDIR)/$(RPKG)-$(GIT_TAG)/man
	cd $(TMPDIR)/$(RPKG)-$(GIT_TAG); "${R_HOME}/bin/R" --slave --file=tools/make-ds.R
	cd $(TMPDIR); NOT_CRAN=false "${R_HOME}/bin/R" CMD build $(RPKG)-$(GIT_TAG) --no-build-vignettes --no-manual
	rm -rf $(TMPDIR)/$(RPKG)-$(GIT_TAG)
	mv $(TMPDIR)/$(RPKG)_$(PKG_VERSION).tar.gz build/$(RPKG)-source.tar.gz
	touch build/r-source-fast

build/r-source-release : $(BIN_OBJS) $(DOC_OBJS) $(SRCS) inst/sbc/sbc_report.html
	install -d build
	git archive --format=tar.gz --prefix $(RPKG)-$(GIT_TAG)/ HEAD > build/$(RPKG)-$(GIT_TAG).tar.gz
	rm -rf build/$(RPKG)-$(GIT_TAG)
	cd build; tar xf $(RPKG)-$(GIT_TAG).tar.gz
	cp -v NAMESPACE build/$(RPKG)-$(GIT_TAG)
	install -d build/$(RPKG)-$(GIT_TAG)/inst/doc
	cp -v inst/doc/$(RPKG).pdf build/$(RPKG)-$(GIT_TAG)/inst/doc
	cp -v inst/sbc/sbc_report.html build/$(RPKG)-$(GIT_TAG)/inst/sbc/sbc_report.html
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

build/installed/$(RPKG)/DESCRIPTION : build/r-source-fast
	rm -rf build/installed
	install -d build/installed
	cd build; $(RCMD) CMD INSTALL --library=./installed --no-docs --no-multiarch --no-test-load --no-clean-on-error $(RPKG)-source.tar.gz

PHONY += dev-install
dev-install: build/installed/$(RPKG)/DESCRIPTION

PHONY += test-all
test-all : $(R_TEST_OBJS)

PHONY += testfast-all
testfast-all : $(R_TESTFAST_OBJS)

PHONY += retestfast-all
retestfast-all : clean-test $(R_TESTFAST_OBJS)

PHONY += retest-all
retest-all : clean-test $(R_TEST_OBJS)

#$(DIR_OBJ)/%.o: %.c $(INCS)
#    mkdir -p $(@D)
#    $(CC) -o $@ $(CFLAGS) -c $< $(INC_DIRS)

PHONY += clean
clean:
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

##
# Debug target that allows you to print a variable
##
print-%  : ; @echo $* = $($*)


.PHONY = $(PHONY)
