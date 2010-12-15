
SUBDIRS =

noinst_PROGRAMS = dune_fem_functionals

dune_fem_functionals_SOURCES = dune_fem_functionals.cc

dune_fem_functionals_CPPFLAGS = $(AM_CPPFLAGS) \
	$(DUNEMPICPPFLAGS) \
	$(UG_CPPFLAGS) \
	$(AMIRAMESH_CPPFLAGS) \
	$(ALBERTA_CPPFLAGS) \
	$(ALUGRID_CPPFLAGS)
# The libraries have to be given in reverse order (most basic libraries
# last).  Also, due to some misunderstanding, a lot of libraries include the
# -L option in LDFLAGS instead of LIBS -- so we have to include the LDFLAGS
# here as well.
dune_fem_functionals_LDADD = \
	$(DUNE_LDFLAGS) $(DUNE_LIBS) \
	$(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) \
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) \
	$(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) \
	$(UG_LDFLAGS) $(UG_LIBS) \
	$(DUNEMPILIBS)	\
	$(LDADD)
dune_fem_functionals_LDFLAGS = $(AM_LDFLAGS) \
	$(DUNEMPILDFLAGS) \
	$(UG_LDFLAGS) \
	$(AMIRAMESH_LDFLAGS) \
	$(ALBERTA_LDFLAGS) \
	$(ALUGRID_LDFLAGS) \
	$(DUNE_LDFLAGS)

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-fem=$(DUNE_FEM_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules

