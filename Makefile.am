# $Id$

# we need the module file to be able to build via dunecontrol
EXTRA_DIST=dune.module

SUBDIRS = m4 dune examples

if BUILD_DOCS
# TODO: set up documentation tree automatically
SUBDIRS += doc
endif

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-fem-localfunctions=$(DUNE_FEM_LOCALFUNCTIONS_ROOT) --with-dune-helper-tools=$(DUNE_HELPER_TOOLS_ROOT) --with-dune-geometry=$(DUNE_GEOMETRY_ROOT) CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules

