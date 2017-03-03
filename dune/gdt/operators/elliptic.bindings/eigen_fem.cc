#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include "../elliptic.bindings.hh"

namespace Dune {
namespace GDT {
namespace bindings {


// these lines have to match the corresponding ones in the .hh header file
#if HAVE_DUNE_FEM
#if HAVE_EIGEN
DUNE_GDT_OPERATORS_ELLIPTIC_BIND_FEM(template, YASP_2D_EQUIDISTANT_OFFSET, eigen_dense);
DUNE_GDT_OPERATORS_ELLIPTIC_BIND_FEM(template, YASP_2D_EQUIDISTANT_OFFSET, eigen_sparse);
#endif
#endif // HAVE_DUNE_FEM

#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID
#if HAVE_DUNE_FEM
#if HAVE_EIGEN
DUNE_GDT_OPERATORS_ELLIPTIC_BIND_FEM(template, ALU_2D_SIMPLEX_CONFORMING, eigen_dense);
DUNE_GDT_OPERATORS_ELLIPTIC_BIND_FEM(template, ALU_2D_SIMPLEX_CONFORMING, eigen_sparse);
#endif
#endif // HAVE_DUNE_FEM
#endif // HAVE_ALUGRID || HAVE_DUNE_ALUGRID


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
