#!/bin/bash

set -e
set -x

source ${SUPERDIR}/scripts/bash/retry_command.bash
cd ${SUPERDIR}
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} configure
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ${BUILD_CMD}
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ${BUILD_CMD} bindings
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ${BUILD_CMD} test_python

cd ${SUPERDIR}/${MY_MODULE}
${DUNE_BUILD_DIR}/${MY_MODULE}/dune-env pip install codecov
${DUNE_BUILD_DIR}/${MY_MODULE}/dune-env codecov -X gcov -F pytest -t ${CODECOV_TOKEN}

