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
#   Felix Schindler (2013, 2016 - 2018)
#   René Fritze     (2016 - 2018)
#   Tobias Leibner  (2016 - 2018)
# ~~~

# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Rene Milk (2017 - 2018)
# ~~~

set -e
SWAP=$HOME/swap.img
MB_SIZE=${1:-4000}
sudo -E dd if=/dev/zero of=${SWAP} bs=1M count=${MB_SIZE}
sudo -E chown root:root ${SWAP}
sudo -E chmod 0600 ${SWAP}
sudo -E mkswap ${SWAP}
sudo -E swapon ${SWAP}
echo "echo 3 > /proc/sys/vm/drop_caches" | sudo -E sh

