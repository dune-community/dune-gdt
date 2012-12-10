
AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(DUNEMPILDFLAGS)
LDADD = $(ALL_PKG_LIBS) $(LOCAL_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS)

noinst_PROGRAMS = lagrange

CUSTOM_CPPFLAGS = -DSGRID -DGRIDDIM=1 -DWORLDDIM=1 -Wall -O0 -DDEBUG -g -fdiagnostics-show-option -std=c++0x

lagrange_SOURCES = lagrange.cc 
lagrange_CPPFLAGS = $(AM_CPPFLAGS) $(EIGEN_CFLAGS) $(CUSTOM_CPPFLAGS)
lagrange_LDFLAGS = $(ALL_LDFLAGS)
lagrange_LDADD = $(LDADD)

include $(top_srcdir)/am/global-rules
