#!/bin/bash
#
# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze    (2018)
#   Tobias Leibner (2018)
# ~~~

set -e
set -x

source ${SUPERDIR}/scripts/bash/retry_command.bash
cd ${SUPERDIR}
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} configure
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ${BUILD_CMD}
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ${BUILD_CMD} bindings
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ${BUILD_CMD} test_python

if [[ ${DRONE_BUILD_EVENT} != "push" ]] ; then
    exit 0
fi

cd ${SUPERDIR}/${MY_MODULE}
${DUNE_BUILD_DIR}/${MY_MODULE}/run-in-dune-env pip install codecov
${DUNE_BUILD_DIR}/${MY_MODULE}/run-in-dune-env codecov -X gcov -F pytest -t ${CODECOV_TOKEN}

