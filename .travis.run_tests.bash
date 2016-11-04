#!/bin/bash

${SRC_DCTRL} ${BLD} --only=${MY_MODULE} configure
if [[ "x${TESTS}" != "xheadercheck" ]]; then
    ${SRC_DCTRL} ${BLD} --only=${MY_MODULE} make
    ${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ninja -v test_binaries_builder_${TESTS}
    ${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ctest -V -j 2 -L "^builder_${TESTS}$"
else
    ${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ninja headercheck
fi

unset GH_TOKEN
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} make install | grep -v "Installing"
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} make package_source
cp -ra ${DUNE_BUILD_DIR}/dune/gdt/test/ ${SUPERDIR}/${MY_MODULE}/test_dir
