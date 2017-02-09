#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include "dg.bindings.hh"

namespace Dune {
namespace GDT {
namespace bindings {


// these lines have to match the corresponding ones in the .hh header file
#if HAVE_DUNE_FEM
DUNE_GDT_SPACES_DG_BIND_FEM(template, YASP_2D_EQUIDISTANT_OFFSET);
#endif

#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID
#if HAVE_DUNE_FEM
DUNE_GDT_SPACES_DG_BIND_FEM(template, ALU_2D_SIMPLEX_CONFORMING);
#endif
#endif // HAVE_ALUGRID || HAVE_DUNE_ALUGRID


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
