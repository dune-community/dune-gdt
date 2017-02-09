#include "config.h"

#if HAVE_DUNE_PYBINDXI
#if HAVE_EIGEN

#include "../elliptic-ipdg.bindings.hh"

namespace Dune {
namespace GDT {
namespace bindings {


DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_GDT(template, YASP_2D_EQUIDISTANT_OFFSET, eigen_dense);
DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_GDT(template, YASP_2D_EQUIDISTANT_OFFSET, eigen_sparse);


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_EIGEN
#endif // HAVE_DUNE_PYBINDXI
