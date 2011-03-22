# DO NOT CHANGE THESE DEFAULT VALUES
#
# use, for example, the following statement instead:
# make GRIDTYPE=ALUGRID_SIMPLEX GRIDDIM=3 POLORDER=2

GRIDTYPE = YASPGRID
GRIDDIM = 2
POLORDER=1

# programs to build 
check_PROGRAMS = test 

# compilation and linking flags
POISSON_CPPFLAGS = -DPOLORDER=$(POLORDER) -DFEMTIMER


# settings for 'poisson'
test_SOURCES = main.cc
test_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(POISSON_CPPFLAGS) $(DUNEMPICPPFLAGS)
test_LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILIBS)


# clean also temporary files 
CLEANFILES = manager.*.log *.dgf.* dgfparser.log

# install these headers 
testdir=$(includedir)/dune/fem/
test_HEADERS= constraints/constraints.hh constraints/localdefault.hh container/sparsitypattern.hh \
						 	container/factory.hh dofvector/dofvector.hh dofvector/localmatrix.hh operator/ellipticfiniteelement.hh \
						 	solver/femassembler.hh subspace/subspace.hh

# files included extra in dist 
EXTRA_DIST=*.dgf parameter 

include $(top_srcdir)/am/global-rules
