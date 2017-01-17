#!/bin/bash

set -e

${SRC_DCTRL} ${BLD} --only=${MY_MODULE} configure
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} make

# this does nothing if all current tests are distributed already, but triggers full build if not
# -> builder will timeout -> manually run refresh_test_timings -> push results
# ${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ninja -v -j 1 refresh_test_timings

${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ninja -v -j 1 test_binaries_builder_${TESTS}
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ctest -V -j 2 -L "^builder_${TESTS}$"
${SUPERDIR}/.travis/init_sshkey.sh ${encrypted_95fb78800815_key} ${encrypted_95fb78800815_iv} keys/dune-community/dune-gdt-testlogs
travis_retry ${SUPERDIR}/scripts/bash/travis_upload_test_logs.bash ${DUNE_BUILD_DIR}/${MY_MODULE}/dune/gdt/test/ || echo Test upload failed

unset GH_TOKEN
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} make install | grep -v "Installing"
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} make package_source

