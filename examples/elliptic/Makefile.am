
AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(DUNEMPILDFLAGS)
LDADD = $(ALL_PKG_LIBS) $(LOCAL_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS)

noinst_PROGRAMS = continuous_galerkin

CUSTOM_CPPFLAGS = -DALUGRID_SIMPLEX -DGRIDDIM=2 -DWORLDDIM=2 -Wall -O0 -DDEBUG -g -fdiagnostics-show-option -std=c++0x

continuous_galerkin_SOURCES = continuous_galerkin.cc 
continuous_galerkin_CPPFLAGS = $(AM_CPPFLAGS) $(EIGEN_CFLAGS) $(CUSTOM_CPPFLAGS) $(BOOST_CPPFLAGS)
continuous_galerkin_LDFLAGS = $(ALL_PKG_LDFLAGS) $(DUNEMPILDFLAGS) -lboost_system -lboost_filesystem
continuous_galerkin_LDADD = $(ALL_PKG_LIBS) $(LOCAL_LIBS)

include $(top_srcdir)/am/global-rules
