
AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(DUNEMPILDFLAGS)
LDADD = $(ALL_PKG_LIBS) $(LOCAL_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS)

noinst_PROGRAMS = continuousgalerkin

CUSTOM_CPPFLAGS = -DALUGRID_SIMPLEX -DGRIDDIM=2 -DWORLDDIM=2 -Wall -O0 -DDEBUG -g -fdiagnostics-show-option -std=c++0x

continuousgalerkin_SOURCES = continuousgalerkin.cc 
continuousgalerkin_CPPFLAGS = $(AM_CPPFLAGS) $(EIGEN_CFLAGS) $(CUSTOM_CPPFLAGS) $(BOOST_CPPFLAGS)
continuousgalerkin_LDFLAGS = $(ALL_LDFLAGS) $(BOOST_LDFLAGS) -lboost_system -lboost_filesystem
continuousgalerkin_LDADD = $(LDADD) $(BOOST_LDADD)

include $(top_srcdir)/am/global-rules
