// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BINDINGS_HH
#define DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/bindings.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/la/type_traits.hh>
#include <dune/xt/la/container.bindings.hh>

#include <dune/gdt/spaces.bindings.hh>
#include <dune/gdt/playground/spaces/block.hh>
#include <dune/gdt/playground/spaces/restricted.hh>
#include <dune/gdt/type_traits.hh>

#include "default.hh"

namespace Dune {
namespace GDT {
namespace bindings {
namespace internal {


template <class S, class V>
class ConstDiscreteFunction
{
  static_assert(is_space<S>::value, "");
  static_assert(XT::LA::is_vector<V>::value, "");
  typedef XT::Grid::extract_grid_t<typename S::GridLayerType> G;

public:
  typedef GDT::ConstDiscreteFunction<S, V> type;

private:
  typedef XT::Functions::LocalizableFunctionInterface<typename S::EntityType,
                                                      typename S::DomainFieldType,
                                                      S::dimDomain,
                                                      typename S::RangeFieldType,
                                                      S::dimRange,
                                                      S::dimRangeCols>
      BaseType;

public:
  typedef pybind11::class_<type, BaseType> bound_type;

  static bound_type bind(pybind11::module& m, const std::string& space_name)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = XT::Common::to_camel_case("const_discrete_function_" + space_name + "_"
                                                     + XT::LA::bindings::container_name<V>::value());

    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def(py::init<const S&, V&, const std::string>(),
          "space"_a,
          "vector"_a,
          "name"_a = "gdt.constdiscretefunction",
          py::keep_alive<1, 2>(),
          py::keep_alive<1, 3>());
    c.def("space", [](type& self) { return self.space(); });
    c.def("vector_copy", [](type& self) { return self.vector(); });
    c.def("visualize",
          [](type& self, const std::string filename, const bool subsampling) {
            return self.visualize(filename, subsampling);
          },
          "filename"_a,
          "subsampling"_a = (S::polOrder > 1));
    // these two are copied from <dune/xt/functions/interfaces.pbh>, would be nicer to inherit them
    c.def("visualize",
          [](const type& self,
             const XT::Grid::GridProvider<G>& grid_provider,
             const std::string& layer,
             const ssize_t lvl,
             const std::string& path,
             const bool subsampling) {
            const auto level = XT::Common::numeric_cast<int>(lvl);
            if (layer == "leaf")
              self.visualize(grid_provider.leaf_view(), path, subsampling);
            else if (layer == "level")
              self.visualize(grid_provider.template layer<XT::Grid::Layers::level, XT::Grid::Backends::view>(level),
                             path,
                             subsampling);
            else
              DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
                         "Given layer has to be one of ('leaf', 'level'), is '" << layer << "'!");
          },
          "grid_provider"_a,
          "layer"_a = "leaf",
          "level"_a = -1,
          "path"_a,
          "subsampling"_a = true);
#if HAVE_DUNE_FEM
    c.def("visualize",
          [](const type& self,
             const XT::Grid::GridProvider<G, XT::Grid::DD::SubdomainGrid<G>>& dd_grid_provider,
             const std::string& layer,
             const ssize_t lvl_or_sbdmn,
             const std::string& path,
             const bool subsampling) {
            const auto level_or_subdomain = XT::Common::numeric_cast<int>(lvl_or_sbdmn);
            if (layer == "leaf")
              self.visualize(dd_grid_provider.leaf_view(), path, subsampling);
            else if (layer == "level")
              self.visualize(dd_grid_provider.template layer<XT::Grid::Layers::level, XT::Grid::Backends::view>(
                                 level_or_subdomain),
                             path,
                             subsampling);
            else if (layer == "dd_subdomain")
              self.visualize(dd_grid_provider.template layer<XT::Grid::Layers::dd_subdomain, XT::Grid::Backends::part>(
                                 level_or_subdomain),
                             path,
                             subsampling);
            else if (layer == "dd_subdomain_oversampled")
              self.visualize(
                  dd_grid_provider.template layer<XT::Grid::Layers::dd_subdomain_oversampled, XT::Grid::Backends::part>(
                      level_or_subdomain),
                  path,
                  subsampling);
            else
              DUNE_THROW(
                  XT::Common::Exceptions::wrong_input_given,
                  "Given layer has to be one of ('leaf', 'level', 'dd_subdomain', 'dd_subdomain_oversampled'), is '"
                      << layer
                      << "'!");
          },
          "dd_grid_provider"_a,
          "layer"_a = "leaf",
          "level_or_subdomain"_a = -1,
          "path"_a,
          "subsampling"_a = true);
#endif // HAVE_DUNE_FEM
    return c;
  } // ... bind(...)
}; // class ConstDiscreteFunction


template <class S, class V>
class DiscreteFunction
{
  static_assert(is_space<S>::value, "");
  static_assert(XT::LA::is_vector<V>::value, "");

  typedef GDT::ConstDiscreteFunction<S, V> BaseType;

public:
  typedef GDT::DiscreteFunction<S, V> type;
  typedef pybind11::class_<type, BaseType> bound_type;

  static bound_type bind(pybind11::module& m, const std::string& space_name)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    bindings::internal::ConstDiscreteFunction<S, V>::bind(m, space_name);

    const auto ClassName = XT::Common::to_camel_case("discrete_function_" + space_name + "_"
                                                     + XT::LA::bindings::container_name<V>::value());

    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def(py::init<const S&, V&, const std::string>(),
          "space"_a,
          "vector"_a,
          "name"_a = "gdt.discretefunction",
          py::keep_alive<1, 2>(),
          py::keep_alive<1, 3>());
    c.def(
        py::init<const S&, const std::string>(), "space"_a, "name"_a = "gdt.discretefunction", py::keep_alive<1, 2>());

    m.def(
        std::string("make_discrete_function").c_str(),
        [](const S& space, V& vector, const std::string& name) { return make_discrete_function(space, vector, name); },
        "space"_a,
        "vector"_a,
        "name"_a = "gdt.discretefunction",
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>());
    m.def(std::string("make_discrete_function").c_str(),
          [](const S& space, const std::string& name) { return make_discrete_function<V>(space, name); },
          "space"_a,
          "name"_a = "gdt.discretefunction",
          py::keep_alive<0, 1>());

    return c;
  } // ... bind(...)

}; // class DiscreteFunction


} // namespace internal


template <class SP, class V>
class DiscreteFunction
{
  typedef typename SP::type S;
  static_assert(is_space<S>::value, "");
  static_assert(XT::LA::is_vector<V>::value, "");
  typedef GDT::ConstDiscreteFunction<S, V> BaseType;

public:
  typedef GDT::DiscreteFunction<S, V> type;
  typedef pybind11::class_<type, BaseType> bound_type;

  template <XT::Grid::Backends backend, XT::Grid::Layers layer>
  static void addbind_restricted(pybind11::module& m, const std::string sp_name)
  {
    try { // we might not be the first to add this
      internal::DiscreteFunction<GDT::RestrictedSpace<S,
                                                      typename XT::Grid::
                                                          Layer<XT::Grid::extract_grid_t<typename S::GridLayerType>,
                                                                layer,
                                                                backend>::type>,
                                 V>::bind(m,
                                          sp_name + "_restricted_to_" + XT::Grid::bindings::layer_name<layer>::value()
                                              + "_"
                                              + XT::Grid::bindings::backend_name<backend>::value());
    } catch (std::runtime_error&) {
    }
  } // ... addbind_restricted(...)

  static bound_type bind(pybind11::module& m)
  {
    const auto sp_name = space_name<SP>::value();
    auto c = internal::DiscreteFunction<S, V>::bind(m, sp_name);

    addbind_restricted<XT::Grid::Backends::part, XT::Grid::Layers::adaptive_leaf>(m, sp_name);
    addbind_restricted<XT::Grid::Backends::part, XT::Grid::Layers::dd_subdomain>(m, sp_name);
    addbind_restricted<XT::Grid::Backends::part, XT::Grid::Layers::dd_subdomain_boundary>(m, sp_name);
    addbind_restricted<XT::Grid::Backends::part, XT::Grid::Layers::dd_subdomain_coupling>(m, sp_name);
    addbind_restricted<XT::Grid::Backends::part, XT::Grid::Layers::dd_subdomain_oversampled>(m, sp_name);
    addbind_restricted<XT::Grid::Backends::part, XT::Grid::Layers::leaf>(m, sp_name);
    addbind_restricted<XT::Grid::Backends::part, XT::Grid::Layers::level>(m, sp_name);
    addbind_restricted<XT::Grid::Backends::view, XT::Grid::Layers::dd_subdomain>(m, sp_name);
    addbind_restricted<XT::Grid::Backends::view, XT::Grid::Layers::dd_subdomain_oversampled>(m, sp_name);
    addbind_restricted<XT::Grid::Backends::view, XT::Grid::Layers::leaf>(m, sp_name);
    addbind_restricted<XT::Grid::Backends::view, XT::Grid::Layers::level>(m, sp_name);

    return c;
  }
}; // class DiscreteFunction


} // namespace bindings
} // namespace GDT
} // namespace Dune


// begin: this is what we need for the .so


#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND(_m, _G, _g_layer, _s_type, _s_backend, _p, _r, _rC, _la)               \
  Dune::GDT::bindings::                                                                                                \
      DiscreteFunction<Dune::GDT::SpaceProvider<_G,                                                                    \
                                                Dune::XT::Grid::Layers::_g_layer,                                      \
                                                Dune::GDT::SpaceType::_s_type,                                         \
                                                Dune::GDT::Backends::_s_backend,                                       \
                                                _p,                                                                    \
                                                double,                                                                \
                                                _r,                                                                    \
                                                _rC>,                                                                  \
                       typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::VectorType>::bind(_m)

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_COMMON(_m, _G, _g_layer, _s_type, _s_backend, _p, _r, _rC)             \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND(_m, _G, _g_layer, _s_type, _s_backend, _p, _r, _rC, common_dense)

#if HAVE_EIGEN
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_EIGEN(_m, _G, _g_layer, _s_type, _s_backend, _p, _r, _rC)              \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND(_m, _G, _g_layer, _s_type, _s_backend, _p, _r, _rC, eigen_dense)
#else
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_EIGEN(_m, _G, _g_layer, _s_type, _s_backend, _p, _r, _rC)
#endif

#if HAVE_DUNE_ISTL
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ISTL(_m, _G, _g_layer, _s_type, _s_backend, _p, _r, _rC)               \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND(_m, _G, _g_layer, _s_type, _s_backend, _p, _r, _rC, istl_dense)
#else
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ISTL(_m, _G, _g_layer, _s_type, _s_backend, _p, _r, _rC)
#endif

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ALL_LA(_m, _G, _g_layer, _s_type, _s_backend, _p, _r, _rC)             \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ISTL(_m, _G, _g_layer, _s_type, _s_backend, _p, _r, _rC)
//  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_COMMON(_m, _G, _g_layer, _s_type, _s_backend, _p, _r, _rC);                  \
//  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_EIGEN(_m, _G, _g_layer, _s_type, _s_backend, _p, _r, _rC);                   \

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_YASP(_m, _g_layer, _s_type, _s_backend, _p, _r, _rC)
//  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ALL_LA(                                                                      \
//      _m, YASP_1D_EQUIDISTANT_OFFSET, _g_layer, _s_type, _s_backend, _p, _r, _rC);                                     \
//  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ALL_LA(                                                                      \
//      _m, YASP_2D_EQUIDISTANT_OFFSET, _g_layer, _s_type, _s_backend, _p, _r, _rC)

#if HAVE_DUNE_ALUGRID
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ALU(_m, _g_layer, _s_type, _s_backend, _p, _r, _rC)                    \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ALL_LA(                                                                      \
      _m, ALU_2D_SIMPLEX_CONFORMING, _g_layer, _s_type, _s_backend, _p, _r, _rC)
#else
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ALU(_m, _g_layer, _s_type, _s_backend, _p, _r, _rC)
#endif

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ALL_GRIDS(_m, _g_layer, _s_type, _s_backend, _p, _r, _rC)              \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_YASP(_m, _g_layer, _s_type, _s_backend, _p, _r, _rC);                        \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ALU(_m, _g_layer, _s_type, _s_backend, _p, _r, _rC)

#if HAVE_DUNE_FEM
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM(_m)                                                                \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ALL_GRIDS(_m, leaf, cg, fem, 1, 1, 1);                                       \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ALL_GRIDS(_m, level, cg, fem, 1, 1, 1);                                      \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ALL_GRIDS(_m, dd_subdomain, cg, fem, 1, 1, 1);                               \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ALL_GRIDS(_m, dd_subdomain, block_cg, fem, 1, 1, 1);                         \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ALL_GRIDS(_m, leaf, dg, fem, 1, 1, 1);                                       \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ALL_GRIDS(_m, level, dg, fem, 1, 1, 1);                                      \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ALL_GRIDS(_m, dd_subdomain, dg, fem, 1, 1, 1);                               \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ALL_GRIDS(_m, dd_subdomain, block_dg, fem, 1, 1, 1)
#else
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM(_m)
#endif

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT(_m)                                                                \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ALL_GRIDS(_m, leaf, fv, gdt, 0, 1, 1);                                       \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ALL_GRIDS(_m, level, fv, gdt, 0, 1, 1)

#if HAVE_DUNE_PDELAB
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB(_m)                                                             \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ALU(_m, leaf, rt, pdelab, 0, 2, 1);                                          \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ALU(_m, level, rt, pdelab, 0, 2, 1)
#else
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB(_m)
#endif

#define DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND(_m)                                                                     \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM(_m);                                                                     \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT(_m);                                                                     \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB(_m)

// end: this is what we need for the .so


#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BINDINGS_HH
