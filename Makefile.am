# This file is part of the dune-gdt project:
#   http://users.dune-project.org/projects/dune-gdt
# Copyright holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

# $Id$

# we need the module file to be able to build via dunecontrol
EXTRA_DIST=dune.module

SUBDIRS = m4

if BUILD_DOCS
# TODO: set up documentation tree automatically
SUBDIRS += doc
endif

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-fem=$(DUNE_FEM_ROOT) --with-dune-stuff=$(DUNE_STUFF_ROOT) --with-dune-geometry=$(DUNE_GEOMETRY_ROOT) CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules

