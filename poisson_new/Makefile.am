# DO NOT CHANGE THESE DEFAULT VALUES
#
# use, for example, the following statement instead:
# make GRIDTYPE=ALUGRID_SIMPLEX GRIDDIM=3 POLORDER=2

GRIDTYPE = YASPGRID
GRIDDIM = 2
POLORDER=1

# programs to build 
check_PROGRAMS = poisson adaptive constraint

# compilation and linking flags
POISSON_CPPFLAGS = -DPOLORDER=$(POLORDER) -DFEMTIMER


# settings for 'poisson'
poisson_SOURCES = main.cc
poisson_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(POISSON_CPPFLAGS) $(DUNEMPICPPFLAGS)
poisson_LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILIBS)

# settings for 'adaptive' (not parallel)
adaptive_SOURCES = mainadaptive.cc
adaptive_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(POISSON_CPPFLAGS)
adaptive_LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)

# settings for 'constraint'
constraint_SOURCES = mainconstraints.cc
constraint_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(POISSON_CPPFLAGS) $(DUNEMPICPPFLAGS)
constraint_LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILIBS)

# clean also temporary files 
CLEANFILES = manager.*.log *.dgf.* dgfparser.log

# install these headers 
poissondir=$(includedir)/tutorial/poisson
poisson_HEADERS= problemdata.hh algorithm.hh algorithmadaptive.hh estimator.hh \
laplaceoperator.hh functional.hh dirichletconstraints.hh neumannconstraints.hh constraints.hh

# files included extra in dist 
EXTRA_DIST=*.dgf parameter 

include $(top_srcdir)/am/global-rules
