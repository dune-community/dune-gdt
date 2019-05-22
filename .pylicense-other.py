# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   René Fritze    (2018)
#   Tobias Leibner (2019)
# ~~~

name = 'This file is part of the dune-gdt project:'
url = 'https://github.com/dune-community/dune-gdt'
copyright_statement = 'Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.'
license = '''Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
          with "runtime exception" (http://www.dune-project.org/license.html)'''
prefix = '#'
lead_in = '# ~~~'
lead_out = '# ~~~'

include_patterns = ('*.txt', '*.cmake', '*.py', '*.sh', '*.bash', '*.dgf', '*.msh', '*.gdb', '*.cfg', '*.travis.*',
                    '*.gitignore', '*.mailmap', '*.gitattributes', '*gitignore-*', '*stamp-vc', '*dune.module',
                    '*Doxylocal', '*.clang-format', '*COPYING-CMAKE-SCRIPTS', '*README', '*LICENSE', '*mainpage',
                    '*switch-build_dir', '*dune-xt-common.pc.in', '*CMakeLists.txt')
exclude_patterns = ('*config.h.cmake', '*.vcsetup*', '*builder_definitions.cmake', '*.ci/shared/*')
