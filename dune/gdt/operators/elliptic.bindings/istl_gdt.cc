#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include "../elliptic.bindings.hh"

namespace Dune {
namespace GDT {
namespace bindings {


// these lines have to match the corresponding ones in the .hh header file
#if HAVE_DUNE_ISTL
DUNE_GDT_OPERATORS_ELLIPTIC_BIND_GDT(template, YASP_2D_EQUIDISTANT_OFFSET, istl_sparse);
#endif

#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID
#if HAVE_DUNE_ISTL
DUNE_GDT_OPERATORS_ELLIPTIC_BIND_GDT(template, ALU_2D_SIMPLEX_CONFORMING, istl_sparse);
#endif
#endif // HAVE_ALUGRID || HAVE_DUNE_ALUGRID


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
