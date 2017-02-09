#include "config.h"

#if HAVE_DUNE_PYBINDXI
#if HAVE_DUNE_ISTL
#if HAVE_DUNE_PDELAB

#include "../elliptic-ipdg.bindings.hh"

namespace Dune {
namespace GDT {
namespace bindings {


DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_FEM(template, YASP_2D_EQUIDISTANT_OFFSET, istl_sparse);


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PDELAB
#endif // HAVE_DUNE_ISTL
#endif // HAVE_DUNE_PYBINDXI
