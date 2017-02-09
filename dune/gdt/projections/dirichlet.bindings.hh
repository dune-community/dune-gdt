// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_PROJECTIONS_DIRICHLET_BINDINGS_HH
#define DUNE_GDT_PROJECTIONS_DIRICHLET_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/grid/grids.bindings.hh>
#include <dune/xt/la/container.bindings.hh>

#include <dune/gdt/spaces.bindings.hh>

#include "dirichlet.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <class GV, class S, class R, class F = double>
class DirichletProjectionLocalizableOperator
{
public:
  typedef GDT::DirichletProjectionLocalizableOperator<GV, S, R, F> type;

private:
  typedef XT::Grid::Walker<GV> BaseType;

public:
  typedef pybind11::class_<type, BaseType> bound_type;

  static bound_type bind(pybind11::module& m, const std::string& space_id, const std::string& la_id)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    bound_type c(m, std::string("DirichletProjectionLocalizableOperator__" + space_id + "__" + la_id).c_str());
    c.def("apply", [](type& self) { self.apply(); });

    m.def(std::string("make_localizable_dirichlet_projection_operator").c_str(),
          [](const XT::Grid::BoundaryInfo<typename XT::Grid::Intersection<GV>::Type>& boundary_info,
             const S& source,
             R& range) {
            return make_localizable_dirichlet_projection_operator(
                       range.space().grid_view(), boundary_info, source, range)
                .release();
          },
          "boundary_info"_a,
          "source"_a,
          "range"_a,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>());

    return c;
  } // ... bind(...)
}; // class DirichletProjectionLocalizableOperator


#define DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_FEM(_prefix, _GRID, _LA)                                                   \
  _prefix DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_(CG_SPACE(_GRID, leaf, fem, 1, 1, 1), _LA);                              \
  _prefix DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_(CG_SPACE(_GRID, level, fem, 1, 1, 1), _LA)

#define DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_PDELAB(_prefix, _GRID, _LA)                                                \
  _prefix DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_(CG_SPACE(_GRID, leaf, pdelab, 1, 1, 1), _LA);                           \
  _prefix DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_(CG_SPACE(_GRID, level, pdelab, 1, 1, 1), _LA)

#define DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_(_SPACE, _LA)                                                              \
  class DirichletProjectionLocalizableOperator<                                                                        \
      typename _SPACE::GridViewType,                                                                                   \
      XT::Functions::LocalizableFunctionInterface<typename _SPACE::EntityType,                                         \
                                                  typename _SPACE::DomainFieldType,                                    \
                                                  _SPACE::dimDomain,                                                   \
                                                  typename _SPACE::RangeFieldType,                                     \
                                                  _SPACE::dimRange,                                                    \
                                                  _SPACE::dimRangeCols>,                                               \
      DiscreteFunction<_SPACE, _LA>>


// these lines have to match the corresponding ones in the .cc source file
#if HAVE_DUNE_FEM
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_FEM(extern template, YASP_2D_EQUIDISTANT_OFFSET, COMMON_DENSE_VECTOR);
#if HAVE_EIGEN
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_FEM(extern template, YASP_2D_EQUIDISTANT_OFFSET, EIGEN_DENSE_VECTOR);
#endif
#if HAVE_DUNE_ISTL
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_FEM(extern template, YASP_2D_EQUIDISTANT_OFFSET, ISTL_DENSE_VECTOR);
#endif
#endif // HAVE_DUNE_FEM
#if HAVE_DUNE_PDELAB
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_PDELAB(extern template, YASP_2D_EQUIDISTANT_OFFSET, COMMON_DENSE_VECTOR);
#if HAVE_EIGEN
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_PDELAB(extern template, YASP_2D_EQUIDISTANT_OFFSET, EIGEN_DENSE_VECTOR);
#endif
#if HAVE_DUNE_ISTL
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_PDELAB(extern template, YASP_2D_EQUIDISTANT_OFFSET, ISTL_DENSE_VECTOR);
#endif
#endif // HAVE_DUNE_PDELAB

#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID
#if HAVE_DUNE_FEM
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_FEM(extern template, ALU_2D_SIMPLEX_CONFORMING, COMMON_DENSE_VECTOR);
#if HAVE_EIGEN
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_FEM(extern template, ALU_2D_SIMPLEX_CONFORMING, EIGEN_DENSE_VECTOR);
#endif
#if HAVE_DUNE_ISTL
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_FEM(extern template, ALU_2D_SIMPLEX_CONFORMING, ISTL_DENSE_VECTOR);
#endif
#endif // HAVE_DUNE_FEM
#if HAVE_DUNE_PDELAB
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_PDELAB(extern template, ALU_2D_SIMPLEX_CONFORMING, COMMON_DENSE_VECTOR);
#if HAVE_EIGEN
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_PDELAB(extern template, ALU_2D_SIMPLEX_CONFORMING, EIGEN_DENSE_VECTOR);
#endif
#if HAVE_DUNE_ISTL
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_PDELAB(extern template, ALU_2D_SIMPLEX_CONFORMING, ISTL_DENSE_VECTOR);
#endif
#endif // HAVE_DUNE_PDELAB
#endif // HAVE_ALUGRID || HAVE_DUNE_ALUGRID


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_PROJECTIONS_DIRICHLET_BINDINGS_HH
