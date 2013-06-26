ADD_EXECUTABLE( examples_elliptic_continuousgalerkin "continuousgalerkin.cc" ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES( examples_elliptic_continuousgalerkin ${COMMON_LIBS} )

ADD_EXECUTABLE( examples_elliptic_discontinuousgalerkin_ip "discontinuousgalerkin-ip.cc" ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES( examples_elliptic_discontinuousgalerkin_ip ${COMMON_LIBS} )

