#!/bin/bash

${SRC_DCTRL} ${BLD} --only=${MY_MODULE} configure
if [[ "x${TESTS}" != "xheadercheck" ]]; then
    ${SRC_DCTRL} ${BLD} --only=${MY_MODULE} make
    travis_wait 50 ${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ninja -v test_binaries_builder_${TESTS}
    then travis_wait 50 ${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ctest -j 2 -L "^builder_${TESTS}$"
    ${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ninja headercheck
fi

${SUPERDIR}/.travis/init_sshkey.sh ${encrypted_95fb78800815_key} ${encrypted_95fb78800815_iv} keys/dune-community/dune-gdt-testlogs
unset GH_TOKEN
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} make install | grep -v "Installing"
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} make package_source
cp -ra ${DUNE_BUILD_DIR}/dune/gdt/test/ ${SUPERDIR}/${MY_MODULE}/test_dir
