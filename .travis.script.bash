#!/bin/bash

set -e

${SRC_DCTRL} ${BLD} --only=${MY_MODULE} configure
if [[ "x${TESTS}" != "xheadercheck" ]]; then
    ${SRC_DCTRL} ${BLD} --only=${MY_MODULE} make
    ${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ninja -v -j 1 test_binaries_builder_${TESTS}
    ${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ctest -V -j 2 -L "^builder_${TESTS}$"
    #${SUPERDIR}/.travis/init_sshkey.sh ${encrypted_95fb78800815_key} ${encrypted_95fb78800815_iv} keys/dune-community/dune-gdt-testlogs
    #${SUPERDIR}/scripts/bash/travis_upload_test_logs.bash ${SUPERDIR}/${MY_MODULE}/${DUNE_BUILD_DIR}/dune/gdt/test/

else
    ${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ninja headercheck
fi

unset GH_TOKEN
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} make install | grep -v "Installing"
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} make package_source

