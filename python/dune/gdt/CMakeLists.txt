# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2018)
#   René Fritze     (2018)
#   Tobias Leibner  (2018)
# ~~~

file(GLOB_RECURSE header "*.hh")

dune_pybindxi_add_module(discretefunction EXCLUDE_FROM_ALL ${header} discretefunction.cc)
dune_pybindxi_add_module(usercode EXCLUDE_FROM_ALL ${header} usercode.cc)
dune_pybindxi_add_module(gamm_2019_talk_on_conservative_rb
                         EXCLUDE_FROM_ALL
                         ${header}
                         gamm-2019-talk-on-conservative-rb.cc)
