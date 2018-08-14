#!/bin/bash

set -e
set -x

WAIT="${SUPERDIR}/scripts/bash/travis_wait_new.bash 45"
source ${SUPERDIR}/scripts/bash/retry_command.bash

${SRC_DCTRL} ${BLD} --only=${MY_MODULE} configure
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ${BUILD_CMD}
${WAIT} ${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ${BUILD_CMD} bindings
${WAIT} ${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ${BUILD_CMD} test_python

cd ${SUPERDIR}/${MY_MODULE}
${DUNE_BUILD_DIR}/${MY_MODULE}/dune-env pip install codecov
${DUNE_BUILD_DIR}/${MY_MODULE}/dune-env codecov -X gcov -F pytest -t ${CODECOV_TOKEN}

