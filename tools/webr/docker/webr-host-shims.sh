#!/bin/sh
#
# Replace webR's host-R redirectors with wrappers onto the native host R.
#
# In the published (amd64) webR image:
#
#   $WEBR_ROOT/wasm/R-<ver>/lib/R/bin/R        bash wrapper, R_HOME_DIR hardcoded
#                                              to /opt/webr/host/R-<ver>/lib/R
#   $WEBR_ROOT/wasm/R-<ver>/lib/R/bin/Rscript  *native amd64 ELF*, execv's the
#                                              same host R
#
# Both are architecture-specific and neither exists for arm64. They are on the
# critical path for RBesT because webr-vars.mk sets R_HOME to the wasm R home
# and RBesT's src/Makevars calls "$(R_HOME)/bin/Rscript" to obtain
# StanHeaders/RcppParallel flags.

set -eu

WEBR_ROOT="${WEBR_ROOT:-/opt/webr}"
R_VERSION="$(cat "${WEBR_ROOT}/R/R-VERSION")"
WASM_R_BIN="${WEBR_ROOT}/wasm/R-${R_VERSION}/lib/R/bin"
HOST_R_HOME="$(R RHOME)"

echo "[webr-host-shims] webR R version : ${R_VERSION}"
echo "[webr-host-shims] wasm R bin     : ${WASM_R_BIN}"
echo "[webr-host-shims] native R_HOME  : ${HOST_R_HOME}"

test -d "${WASM_R_BIN}"

rm -f "${WASM_R_BIN}/R" "${WASM_R_BIN}/Rscript" "${WASM_R_BIN}/Rscript.orig"

cat > "${WASM_R_BIN}/R" <<EOF
#!/bin/sh
# Native replacement for webR's amd64 host-R wrapper (see webr-host-shims.sh).
exec "${HOST_R_HOME}/bin/R" "\$@"
EOF

cat > "${WASM_R_BIN}/Rscript" <<EOF
#!/bin/sh
# Native replacement for webR's amd64 Rscript ELF (see webr-host-shims.sh).
exec "${HOST_R_HOME}/bin/Rscript" "\$@"
EOF

chmod 0755 "${WASM_R_BIN}/R" "${WASM_R_BIN}/Rscript"

# Anything still referring to the host tree by its absolute path gets the
# native R instead of a missing amd64 one.
mkdir -p "${WEBR_ROOT}/host/R-${R_VERSION}/lib"
ln -sfn "${HOST_R_HOME}" "${WEBR_ROOT}/host/R-${R_VERSION}/lib/R"

"${WASM_R_BIN}/Rscript" -e 'cat("[webr-host-shims] shim runs:", R.version.string, "\n")'
