```
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2010, 2013 - 2014, 2016 - 2017)
#   Rene Milk       (2017)
```

[![Build Status](https://travis-ci.org/dune-community/dune-gdt.png?branch=master)](https://travis-ci.org/dune-community/dune-gdt)

dune-gdt is a DUNE (http://www.dune-project.org) module which provides a
generic discretization toolbox for grid-based numerical methods. It contains
building blocks - like local operators, local evaluations, local assemblers -
for discretization methods as well as generic interfaces for objects like
discrete function spaces and basefunction sets. Implementations are provided
using the main DUNE discretization modules, like dune-fem
(http://dune.mathematik.uni-freiburg.de/), dune-fem-localfunctions
(http://users.dune-project.org/projects/dune-fem-localfunctions) and
dune-pdelab (http://www.dune-project.org/pdelab/).

New users may best try out this module by using the git supermodule
dune-gdt-super (https://github.com/dune-community/dune-gdt-super).
Experienced DUNE users may go ahead. As usual, you will have to call
configure and make using dunecontrol
(see http://www.dune-project.org/doc/installation-notes.html), working
examples are located in 'dune/gdt/test/'...

If you want to start hacking go ahead and fork us on github.com:
https://github.com/dune-community/dune-gdt/
