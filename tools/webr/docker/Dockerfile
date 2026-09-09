# Native-architecture webR/Emscripten build environment for RBesT.
#
# This is "option C" of design/howto-build-rbest-webr.md (repo-root relative):
# a toolchain that runs natively on the host architecture (verified on
# arm64), avoiding both Rosetta and qemu. Unlike design/wasm-build/Dockerfile
# (repo-root relative) it does NOT use ghcr.io/r-wasm/webr as a runtime base
# -- that image is amd64-only. Instead the environment is assembled from
# three pieces:
#
#   1. host R, native, from the multi-arch rocker/r-ver image;
#   2. Emscripten, native, from the multi-arch emscripten/emsdk image, pinned to
#      the exact version the webR image ships;
#   3. the architecture-independent half of the webR tree, lifted out of the
#      amd64 webR image with COPY --from. Nothing is ever *run* in that stage,
#      so no emulation is involved.
#
# The same Dockerfile therefore builds on amd64 and arm64 with no branching:
# stages 1 and 2 resolve against the *build* platform via their manifest lists,
# and stage 3 is architecture-independent content.

# --- pins ------------------------------------------------------------------
#
# EMSDK_VERSION must match the Emscripten in WEBR_TAG's image. For v0.6.0 that
# is 5.0.7, read from /opt/emsdk/upstream/emscripten/emscripten-version.txt.
# Note that emscripten/emsdk only publishes arm64 for recent tags (4.0.8 is
# amd64-only), so an older webR pin would need `emsdk install` instead.
ARG WEBR_TAG=v0.6.0
ARG EMSDK_VERSION=5.0.7
ARG R_VERSION=4.6.0

# --- stage: webR tree (never executed, pinned to the published arch) --------
FROM --platform=linux/amd64 ghcr.io/r-wasm/webr:${WEBR_TAG} AS webr

# --- stage: Emscripten (resolves to the build platform's architecture) ------
FROM emscripten/emsdk:${EMSDK_VERSION} AS emsdk

# --- stage: the build environment ------------------------------------------
FROM rocker/r-ver:${R_VERSION}

ARG WEBR_TAG
ARG EMSDK_VERSION
ARG R_VERSION

# Corporate proxies: docker does *not* forward these into an image build, and
# the resulting error is misleading -- "package 'pak' is not available for this
# version of R" actually means "the CRAN index could not be downloaded".
ARG http_proxy=""
ARG https_proxy=""
ARG no_proxy=""
ARG HTTP_PROXY=""
ARG HTTPS_PROXY=""
ARG NO_PROXY=""
ENV http_proxy=${http_proxy} \
    https_proxy=${https_proxy} \
    no_proxy=${no_proxy} \
    HTTP_PROXY=${HTTP_PROXY} \
    HTTPS_PROXY=${HTTPS_PROXY} \
    NO_PROXY=${NO_PROXY}

ARG APT_PROXY=""

# --- Corporate TLS interception (optional) ---------------------------------
#
# Some networks run an intercepting proxy that terminates and re-signs every
# outbound TLS connection, so without corporate roots in the trust store every
# https fetch inside the build fails. Same mechanism as
# design/wasm-build/Dockerfile (repo-root relative): CA material is supplied
# via an optional, named build context (`corpcerts`), installed if present, a
# no-op if not. See compose.yaml
# (`build.additional_contexts`, `CORP_CERTS_DIR`) or invoke directly with
#   docker buildx build --build-context corpcerts=/path/to/certs .
COPY --from=corpcerts . /tmp/corporate-certs/

RUN set -eux; \
    if [ -n "${APT_PROXY}" ]; then \
      echo "Acquire::http::Proxy \"${APT_PROXY}\";" > /etc/apt/apt.conf.d/01proxy-build; \
      echo "Acquire::https::Proxy \"${APT_PROXY}\";" >> /etc/apt/apt.conf.d/01proxy-build; \
    fi; \
    for cert in /tmp/corporate-certs/*.pem; do \
      [ -e "$cert" ] || continue; \
      install -m 0644 "$cert" \
        "/usr/local/share/ca-certificates/$(basename "$cert" .pem).crt"; \
    done; \
    rm -rf /tmp/corporate-certs; \
    update-ca-certificates --fresh

ENV SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt \
    CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt \
    NODE_EXTRA_CA_CERTS=/etc/ssl/certs/ca-certificates.crt \
    REQUESTS_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt \
    GIT_SSL_CAINFO=/etc/ssl/certs/ca-certificates.crt

RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
      build-essential git ca-certificates curl xz-utils file \
      pkg-config python3 cmake make \
      libcurl4-openssl-dev libssl-dev libxml2-dev zlib1g-dev; \
    rm -rf /var/lib/apt/lists/*

# --- Emscripten, native for this architecture ------------------------------
#
# Relocated from /emsdk to /opt/emsdk because webR's wasm `Makeconf` bakes in
# absolute /opt/emsdk/... paths (CC=, PKG_CONFIG_LIBDIR=.../cache/sysroot/...).
# emsdk's own config is relocation-safe: .emscripten derives every path from
# dirname($EM_CONFIG).
COPY --from=emsdk /emsdk /opt/emsdk

ENV EMSDK=/opt/emsdk \
    EM_CONFIG=/opt/emsdk/.emscripten \
    EMSDK_QUIET=1 \
    PATH=/opt/emsdk:/opt/emsdk/upstream/emscripten:/opt/emsdk/upstream/bin:${PATH}

# --- The architecture-independent half of webR -----------------------------
#
# Required by rwasm, over and above the `wasm` tree that the design note
# mentions:
#   R/R-VERSION   read by rwasm:::.onLoad
#   R/build       the R source tree; webr-vars.mk adds its build/include and
#                 src/include to CPPFLAGS
#   host/bin      shell/perl/python *-config scripts, put on PATH by rwasm
#   tools         the pkg-config shim used via EM_PKG_CONFIG, and fortran.mk
COPY --from=webr /opt/webr/wasm /opt/webr/wasm
COPY --from=webr /opt/webr/R/R-VERSION /opt/webr/R/R-VERSION
COPY --from=webr /opt/webr/R/build /opt/webr/R/build
COPY --from=webr /opt/webr/host/bin /opt/webr/host/bin
COPY --from=webr /opt/webr/tools /opt/webr/tools
COPY --from=webr /opt/webr/libs /opt/webr/libs
COPY --from=webr /opt/webr/packages /opt/webr/packages

# The wasm sysroot built by webR (wasm artefacts, architecture-independent).
# Taking it wholesale guarantees the emscripten ports RBesT links against are
# byte-identical to the ones webR's own R was built with, instead of being
# rebuilt on demand from a possibly different port recipe.
COPY --from=webr /opt/emsdk/upstream/emscripten/cache /opt/emsdk/upstream/emscripten/cache

ENV WEBR_ROOT=/opt/webr \
    PKG_USE_BIOCONDUCTOR=false

# --- Native replacements for webR's host-R redirectors ---------------------
#
# $WEBR_ROOT/wasm/R-<ver>/lib/R/bin/Rscript is a *native amd64 ELF* that execv's
# /opt/webr/host/R-<ver>/lib/R, and bin/R is a bash wrapper with the same path
# hardcoded. Neither runs on arm64. They matter because webr-vars.mk sets
# R_HOME to the wasm R home, and RBesT's src/Makevars shells out to
# "$(R_HOME)/bin/Rscript" for StanHeaders/RcppParallel flags.
#
# Replace both with wrappers onto the native host R, and point the
# /opt/webr/host/R-<ver> path at it too for anything else that references it.
COPY webr-host-shims.sh /usr/local/bin/webr-host-shims.sh
RUN /usr/local/bin/webr-host-shims.sh

# --- rwasm and the native LinkingTo closure --------------------------------
ARG CRAN_MIRROR=https://cloud.r-project.org
ARG RWASM_REF=HEAD
ENV CRAN_MIRROR=${CRAN_MIRROR} \
    RWASM_REF=${RWASM_REF}

COPY proxy-env.R /usr/local/bin/proxy-env.R
COPY install-rwasm.R /usr/local/bin/install-rwasm.R
RUN R -q -f /usr/local/bin/install-rwasm.R

COPY --from=webrtools patch-stanheaders.R /usr/local/bin/patch-stanheaders.R

# The wasm cross-compile runs `R CMD INSTALL` under *host* R, so LinkingTo
# headers are resolved from the host library. rstan/StanHeaders come from CRAN,
# at versions derived from what repo.r-wasm.org publishes so that the host
# headers are the same release as the wasm StanHeaders binary the image ships.
COPY install-host-deps.R /usr/local/bin/install-host-deps.R
RUN R -q -f /usr/local/bin/install-host-deps.R

# --- Self-check ------------------------------------------------------------
COPY check-toolchain.sh /usr/local/bin/check-toolchain.sh
RUN /usr/local/bin/check-toolchain.sh

WORKDIR /out
COPY fetch-wasm-deps.R /usr/local/bin/fetch-wasm-deps.R
# The TBB patch helpers live one directory up, outside this build context; they
# come in through the `webrtools` additional context declared in compose.yaml.
# They cannot be read from the /src mount at run time either -- `tools/` is
# .Rbuildignore'd, so it is absent from a staged source package. The stub source
# itself lives under `inst/`, so it does travel in the staged package and is
# read from /src.
COPY --from=webrtools tbb-patch-common.R /usr/local/bin/tbb-patch-common.R
COPY --from=webrtools patch-rstan-tarball.R /usr/local/bin/patch-rstan-tarball.R
COPY build-rbest-wasm.R /usr/local/bin/build-rbest-wasm.R

CMD ["R", "-q", "-f", "/usr/local/bin/build-rbest-wasm.R"]
