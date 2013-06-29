# This file is part of the dune-gdt project:
#   http://users.dune-project.org/projects/dune-gdt
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(DUNEMPILDFLAGS)
LDADD = $(ALL_PKG_LIBS) $(LOCAL_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS)

noinst_PROGRAMS = continuousgalerkin discontinuousgalerkin_ip discontinuousgalerkin_swip

CUSTOM_CPPFLAGS = -DALUGRID_SIMPLEX -DGRIDDIM=2 -DWORLDDIM=2 -Wall -O0 -DDEBUG -g -fdiagnostics-show-option -std=c++0x

continuousgalerkin_SOURCES = continuousgalerkin.cc 
continuousgalerkin_CPPFLAGS = $(AM_CPPFLAGS) $(EIGEN_CFLAGS) $(CUSTOM_CPPFLAGS) $(BOOST_CPPFLAGS)
continuousgalerkin_LDFLAGS = $(ALL_LDFLAGS) $(BOOST_LDFLAGS) -lboost_system -lboost_filesystem
continuousgalerkin_LDADD = $(LDADD) $(BOOST_LDADD)

discontinuousgalerkin_ip_SOURCES = discontinuousgalerkin-ip.cc 
discontinuousgalerkin_ip_CPPFLAGS = $(AM_CPPFLAGS) $(EIGEN_CFLAGS) $(CUSTOM_CPPFLAGS) $(BOOST_CPPFLAGS)
discontinuousgalerkin_ip_LDFLAGS = $(ALL_LDFLAGS) $(BOOST_LDFLAGS) -lboost_system -lboost_filesystem
discontinuousgalerkin_ip_LDADD = $(LDADD) $(BOOST_LDADD)

discontinuousgalerkin_swip_SOURCES = discontinuousgalerkin-swip.cc 
discontinuousgalerkin_swip_CPPFLAGS = $(AM_CPPFLAGS) $(EIGEN_CFLAGS) $(CUSTOM_CPPFLAGS) $(BOOST_CPPFLAGS)
discontinuousgalerkin_swip_LDFLAGS = $(ALL_LDFLAGS) $(BOOST_LDFLAGS) -lboost_system -lboost_filesystem
discontinuousgalerkin_swip_LDADD = $(LDADD) $(BOOST_LDADD)

include $(top_srcdir)/am/global-rules
