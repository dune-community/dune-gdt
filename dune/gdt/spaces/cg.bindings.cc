#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include "cg.bindings.hh"

namespace Dune {
namespace GDT {
namespace bindings {


// these lines have to match the corresponding ones in the .hh header file
#if HAVE_DUNE_FEM
DUNE_GDT_SPACES_CG_BIND_FEM(template, YASP_2D_EQUIDISTANT_OFFSET);
#endif
#if HAVE_DUNE_PDELAB
DUNE_GDT_SPACES_CG_BIND_PDELAB(template, YASP_2D_EQUIDISTANT_OFFSET);
#endif

#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID
#if HAVE_DUNE_FEM
DUNE_GDT_SPACES_CG_BIND_FEM(template, ALU_2D_SIMPLEX_CONFORMING);
#endif
#if HAVE_DUNE_PDELAB
DUNE_GDT_SPACES_CG_BIND_PDELAB(template, ALU_2D_SIMPLEX_CONFORMING);
#endif
#endif // HAVE_ALUGRID || HAVE_DUNE_ALUGRID


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
