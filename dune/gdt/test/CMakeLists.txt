# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2013, 2016 - 2017)
#   Rene Milk       (2016 - 2018)
#   Tobias Leibner  (2016 - 2017)
# ~~~

enable_testing()

set(DXT_TEST_TIMEOUT 1000)
if(DUNE_XT_WITH_PYTHON_BINDINGS)
  begin_testcases(dunegdt)
else()
  begin_testcases()
endif()

end_testcases()
