# This file is part of the dune-gdt project:
#   http://users.dune-project.org/projects/dune-gdt
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

ADD_EXECUTABLE( examples_elliptic_continuousgalerkin "continuousgalerkin.cc" ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES( examples_elliptic_continuousgalerkin ${COMMON_LIBS} )

ADD_EXECUTABLE( examples_elliptic_discontinuousgalerkin_ip "discontinuousgalerkin-ip.cc" ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES( examples_elliptic_discontinuousgalerkin_ip ${COMMON_LIBS} )

ADD_EXECUTABLE( examples_elliptic_discontinuousgalerkin_swip "discontinuousgalerkin-swip.cc" ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES( examples_elliptic_discontinuousgalerkin_swip ${COMMON_LIBS} )

ADD_EXECUTABLE( examples_elliptic_convergence_study "convergence_study_main.cc" ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES( examples_elliptic_convergence_study ${COMMON_LIBS} )
