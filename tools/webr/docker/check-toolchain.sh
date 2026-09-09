#!/bin/sh
#
# Fail the image build early and loudly if the assembled toolchain is not what
# option C assumes: native compilers, the pinned Emscripten version, a
# relocated-but-working emsdk, and an rwasm that finds both trees.

set -eu

fail() { echo "[check-toolchain] FAIL: $*" >&2; exit 1; }
ok()   { echo "[check-toolchain] ok: $*"; }

ARCH="$(uname -m)"
echo "[check-toolchain] build architecture: ${ARCH}"

# 1. The Emscripten binaries must be native, not the amd64 ones we deliberately
#    did not take from the webR image.
CLANG="${EMSDK}/upstream/bin/clang"
test -x "${CLANG}" || fail "no clang at ${CLANG}"
CLANG_KIND="$(file -bL "${CLANG}")"
echo "[check-toolchain] clang: ${CLANG_KIND}"
case "${ARCH}" in
  aarch64|arm64) echo "${CLANG_KIND}" | grep -q "aarch64\|ARM aarch64" || fail "clang is not aarch64: ${CLANG_KIND}" ;;
  x86_64)        echo "${CLANG_KIND}" | grep -q "x86-64"               || fail "clang is not x86-64: ${CLANG_KIND}" ;;
esac
ok "Emscripten LLVM is native for ${ARCH}"

# 2. Version must match the webR image's Emscripten, or the wasm side-module
#    ABI is not guaranteed to line up with the prebuilt wasm R.
EMCC_VERSION="$(emcc --version | head -1 | sed -E 's/.* ([0-9]+\.[0-9]+\.[0-9]+) .*/\1/')"
echo "[check-toolchain] emcc version: ${EMCC_VERSION}"
if [ -n "${EMSDK_VERSION:-}" ] && [ "${EMCC_VERSION}" != "${EMSDK_VERSION}" ]; then
  fail "emcc ${EMCC_VERSION} != pinned ${EMSDK_VERSION}"
fi
ok "Emscripten version pinned at ${EMCC_VERSION}"

# 3. The relocation from /emsdk to /opt/emsdk must not have broken node lookup.
NODE_JS="$(python3 - <<'PY'
import os, re
cfg = os.environ["EM_CONFIG"]
src = open(cfg).read()
ns = {"os": os}
exec(src, ns)
print(ns["NODE_JS"])
PY
)"
test -x "${NODE_JS}" || fail "emsdk node not found after relocation: ${NODE_JS}"
ok "emsdk node: ${NODE_JS} ($(${NODE_JS} --version))"

# 4. A real compile+link, to catch a broken sysroot before an hour of R builds.
TMP="$(mktemp -d)"
cat > "${TMP}/hello.cpp" <<'EOF'
#include <cstdio>
int main() { std::printf("emscripten ok\n"); return 0; }
EOF
emcc "${TMP}/hello.cpp" -o "${TMP}/hello.js" >/dev/null
"${NODE_JS}" "${TMP}/hello.js" | grep -q "emscripten ok" || fail "emcc hello world did not run"
rm -rf "${TMP}"
ok "emcc compiles and links a wasm program"

# 5. webR tree completeness, as required by rwasm and webr-vars.mk.
R_VERSION="$(cat "${WEBR_ROOT}/R/R-VERSION")"
for p in \
  "${WEBR_ROOT}/R/R-VERSION" \
  "${WEBR_ROOT}/R/build/R-${R_VERSION}/build/include" \
  "${WEBR_ROOT}/R/build/R-${R_VERSION}/src/include" \
  "${WEBR_ROOT}/wasm/R-${R_VERSION}/lib/R/etc/Makeconf" \
  "${WEBR_ROOT}/wasm/lib/pkgconfig" \
  "${WEBR_ROOT}/host/bin" \
  "${WEBR_ROOT}/tools/shims/pkg-config"
do
  test -e "${p}" || fail "missing from webR tree: ${p}"
done
ok "webR tree complete for R ${R_VERSION}"

# 6. rwasm must find both roots at load time.
R -q -e 'library(rwasm)' 2>&1 | sed 's/^/[check-toolchain] /'
R -q -s -e 'loadNamespace("rwasm"); stopifnot(!is.null(getOption("rwasm.webr_root")), !is.null(getOption("rwasm.emscripten_root")))' \
  || fail "rwasm did not pick up WEBR_ROOT / EMSDK"
ok "rwasm configured"

echo "[check-toolchain] all checks passed"
