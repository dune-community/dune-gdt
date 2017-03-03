#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include "fv.bindings.hh"

namespace Dune {
namespace GDT {
namespace bindings {


// these lines have to match the corresponding ones in the .hh header file
DUNE_GDT_SPACES_FV_BIND_GDT(template, YASP_2D_EQUIDISTANT_OFFSET);

#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID
DUNE_GDT_SPACES_FV_BIND_GDT(template, ALU_2D_SIMPLEX_CONFORMING);
#endif // HAVE_ALUGRID || HAVE_DUNE_ALUGRID

} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
