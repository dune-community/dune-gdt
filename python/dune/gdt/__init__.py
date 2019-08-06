# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017 - 2018)
#   René Fritze     (2016, 2018)
# ~~~

from dune.xt import _register_special_funcs

for mod_name in (
        'discretefunction',
        'gamm_2019_talk_on_conservative_rb',
        'tools',
        'usercode',
        ):
    _register_special_funcs(mod_name, 'dune.gdt')

