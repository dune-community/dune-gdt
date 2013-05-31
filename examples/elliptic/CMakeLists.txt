ADD_EXECUTABLE( examples_elliptic_continuousgalerkin "continuousgalerkin.cc" ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES( examples_elliptic_continuousgalerkin ${COMMON_LIBS} )

ADD_EXECUTABLE( examples_elliptic_discontinuousgalerkin "discontinuousgalerkin.cc" ${COMMON_HEADER} )
TARGET_LINK_LIBRARIES( examples_elliptic_discontinuousgalerkin ${COMMON_LIBS} )

