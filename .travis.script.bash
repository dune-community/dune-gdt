#!/bin/bash
#
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017)
#   Rene Milk       (2016 - 2017)

set -ex

WAIT="${SUPERDIR}/scripts/bash/travis_wait_new.bash 45"
source ${SUPERDIR}/scripts/bash/retry_command.bash

${SRC_DCTRL} ${BLD} --only=${MY_MODULE} configure
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} make

# this does nothing if all current tests are distributed already, but triggers full build if not
# -> builder will timeout -> manually run refresh_test_timings -> push results
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} bexec ninja -v -j 1 refresh_test_timings

free -h



#unset GH_TOKEN
#${SRC_DCTRL} ${BLD} --only=${MY_MODULE} make install | grep -v "Installing"
#${SRC_DCTRL} ${BLD} --only=${MY_MODULE} make package_source

