
GRIDDIM=2
GRIDTYPE=YASPGRID
POLORDER=1

noinst_PROGRAMS = l2functional_test discretelinearfunctional_test

l2functional_test_SOURCES = l2functional_test.cc

l2functional_test_CPPFLAGS = $(AM_CPPFLAGS) \
	$(DUNEMPICPPFLAGS) \
	$(UG_CPPFLAGS) \
	$(AMIRAMESH_CPPFLAGS) \
	$(ALBERTA_CPPFLAGS) \
	$(ALUGRID_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) -D$(GRIDTYPE) -DPOLORDER=$(POLORDER) -Wall -O0 -DDEBUG -funroll-loops -g -ggdb -fno-strict-aliasing -std=c++0x
# The libraries have to be given in reverse order (most basic libraries
# last).  Also, due to some misunderstanding, a lot of libraries include the
# -L option in LDFLAGS instead of LIBS -- so we have to include the LDFLAGS
# here as well.
l2functional_test_LDADD = \
	$(DUNE_LDFLAGS) $(DUNE_LIBS) \
	$(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) \
	$(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) \
	$(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) \
	$(UG_LDFLAGS) $(UG_LIBS) \
	$(DUNEMPILIBS)	\
	$(LDADD)
l2functional_test_LDFLAGS = $(AM_LDFLAGS) \
	$(DUNEMPILDFLAGS) \
	$(UG_LDFLAGS) \
	$(AMIRAMESH_LDFLAGS) \
	$(ALBERTA_LDFLAGS) \
	$(ALUGRID_LDFLAGS) \
	$(DUNE_LDFLAGS)

discretelinearfunctional_test_SOURCES = discretelinearfunctional_test.cc

discretelinearfunctional_test_CPPFLAGS = $(AM_CPPFLAGS) \
        $(DUNEMPICPPFLAGS) \
        $(UG_CPPFLAGS) \
        $(AMIRAMESH_CPPFLAGS) \
        $(ALBERTA_CPPFLAGS) \
        $(ALUGRID_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) -D$(GRIDTYPE) -DPOLORDER=$(POLORDER) -Wall -O0 -DDEBUG -funroll-loops -g -ggdb -fno-strict-aliasing -std=c++0x
# The libraries have to be given in reverse order (most basic libraries
# last).  Also, due to some misunderstanding, a lot of libraries include the
# -L option in LDFLAGS instead of LIBS -- so we have to include the LDFLAGS
# here as well.
discretelinearfunctional_test_LDADD = \
        $(DUNE_LDFLAGS) $(DUNE_LIBS) \
        $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) \
        $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) \
        $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) \
        $(UG_LDFLAGS) $(UG_LIBS) \
        $(DUNEMPILIBS)  \
        $(LDADD)
discretelinearfunctional_test_LDFLAGS = $(AM_LDFLAGS) \
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

