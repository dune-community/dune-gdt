#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include "../elliptic.bindings.hh"

namespace Dune {
namespace GDT {
namespace bindings {


// these lines have to match the corresponding ones in the .hh header file
DUNE_GDT_OPERATORS_ELLIPTIC_BIND_GDT(template, YASP_2D_EQUIDISTANT_OFFSET, common_dense);

#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID
DUNE_GDT_OPERATORS_ELLIPTIC_BIND_GDT(template, ALU_2D_SIMPLEX_CONFORMING, common_dense);
#endif // HAVE_ALUGRID || HAVE_DUNE_ALUGRID


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
