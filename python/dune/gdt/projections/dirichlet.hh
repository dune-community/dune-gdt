// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef PYTHON_DUNE_GDT_PROJECTIONS_DIRICHLET_BINDINGS_HH
#define PYTHON_DUNE_GDT_PROJECTIONS_DIRICHLET_BINDINGS_HH

#include <dune/pybindxi/pybind11.h>

#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/la/container.bindings.hh>

#include <python/dune/gdt/spaces/interface.hh>
#include <dune/gdt/type_traits.hh>

#include <dune/gdt/projections/dirichlet.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class SP, class V>
class DirichletProjectionLocalizableOperator
{
  typedef typename SP::type SpaceType;
  static_assert(is_space<SpaceType>::value, "");
  static_assert(XT::LA::is_vector<V>::value, "");
  typedef typename SpaceType::GridLayerType GridLayerType;
  typedef typename XT::Functions::LocalizableFunctionInterface<typename SpaceType::EntityType,
                                                               typename SpaceType::DomainFieldType,
                                                               SpaceType::dimDomain,
                                                               typename SpaceType::RangeFieldType,
                                                               SpaceType::dimRange,
                                                               SpaceType::dimRangeCols>
      SourceType;
  typedef DiscreteFunction<SpaceType, V> RangeType;
  typedef XT::Grid::Walker<GridLayerType> BaseType;

public:
  typedef GDT::DirichletProjectionLocalizableOperator<GridLayerType, SourceType, RangeType, double> type;
  typedef pybind11::class_<type, BaseType> bound_type;

  static bound_type bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName =
        XT::Common::to_camel_case("dirichlet_projection_localizable_operator_" + space_name<SP>::value() + "_"
                                  + XT::LA::bindings::container_name<V>::value());

    bound_type c(m, ClassName.c_str());
    c.def("apply", [](type& self) { self.apply(); });

    m.def(std::string("make_localizable_dirichlet_projection_operator").c_str(),
          [](const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info,
             const SourceType& source,
             RangeType& range) {
            return make_localizable_dirichlet_projection_operator(
                       range.space().grid_layer(), boundary_info, source, range)
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


// begin: this is what we need for the .so

#define _DUNE_GDT_PROJECTIONS_DIRICHLET_BIND(_m, _GRID, _layer, _backend, _r, _rC, _la)                                \
  Dune::GDT::bindings::                                                                                                \
      DirichletProjectionLocalizableOperator<Dune::GDT::CgSpaceProvider<_GRID,                                         \
                                                                        Dune::XT::Grid::Layers::_layer,                \
                                                                        Dune::GDT::Backends::_backend,                 \
                                                                        1,                                             \
                                                                        double,                                        \
                                                                        _r,                                            \
                                                                        _rC>,                                          \
                                             typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::   \
                                                 VectorType>::bind(_m)

#define _DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_GDT(_m, _la)                                                              \
  _DUNE_GDT_PROJECTIONS_DIRICHLET_BIND(_m, GDT_BINDINGS_GRID, leaf, gdt, _la);


#define _DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_BACKENDS(_m, _la) _DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_GDT(_m, _la);

#define _DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_COMMON(_m) _DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_BACKENDS(_m, common_dense)

#if HAVE_EIGEN
#define _DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_EIGEN(_m) _DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_BACKENDS(_m, eigen_dense)
#else
#define _DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_EIGEN(_m)
#endif

#define _DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_ISTL(_m) _DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_BACKENDS(_m, istl_dense)

#define DUNE_GDT_PROJECTIONS_DIRICHLET_BIND(_m) _DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_ISTL(_m)
//  _DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_COMMON(_m);                                                                     \
//  _DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_EIGEN(_m);                                                                      \

// end: this is what we need for the .so


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // PYTHON_DUNE_GDT_PROJECTIONS_DIRICHLET_BINDINGS_HH
