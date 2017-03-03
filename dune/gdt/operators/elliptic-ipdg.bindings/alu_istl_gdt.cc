#include "config.h"

#if HAVE_DUNE_PYBINDXI
#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID
#if HAVE_DUNE_ISTL

#include "../elliptic-ipdg.bindings.hh"

namespace Dune {
namespace GDT {
namespace bindings {


DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_GDT(template, ALU_2D_SIMPLEX_CONFORMING, istl_sparse);


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_ISTL
#endif // HAVE_ALUGRID || HAVE_DUNE_ALUGRID
#endif // HAVE_DUNE_PYBINDXI
