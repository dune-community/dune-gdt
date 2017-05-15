// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BINDINGS_HH
#define DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/la/container.hh>
#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/layers.bindings.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>

#include <dune/gdt/spaces.hh>
#include <dune/gdt/spaces.bindings.hh>

#include "integrands/elliptic-ipdg.hh"
#include "operators/integrals.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <class DF, typename DT, class SP, XT::Grid::Layers layer, LocalEllipticIpdgIntegrands::Method method>
class LocalEllipticIpdgInnerIntegralOperator
{
  static_assert(XT::Functions::is_localizable_function<DF>::value, "");
  static_assert(XT::Functions::is_localizable_function<DT>::value || std::is_same<DT, void>::value, "");
  typedef typename SP::type S;
  static_assert(is_space<S>::value, "");
  typedef XT::Grid::extract_grid_t<typename S::GridLayerType> G;
  typedef typename S::BaseFunctionSetType B;
  typedef XT::Grid::extract_intersection_t<typename XT::Grid::Layer<G,
                                                                    layer,
                                                                    S::layer_backend
#if HAVE_DUNE_FEM
                                                                    ,
                                                                    XT::Grid::DD::SubdomainGrid<G>
#endif
                                                                    >::type>
      I;

public:
  typedef LocalCouplingIntegralOperator<LocalEllipticIpdgIntegrands::Inner<DF, DT, method>, B, I> type;
  typedef pybind11::class_<type> bound_type;

private:
  template <bool single_diffusion = std::is_same<DT, void>::value,
            bool scalar = (DF::dimRange == 1 && DF::dimRangeCols == 1),
            bool anything = false>
  struct diffusion_switch
  {
    static std::string suffix()
    {
      return "diffusion_factor_and_tensor";
    }

    static void addbind_factory_methods(pybind11::module& m)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      const std::string method_name = "make_local_elliptic_" + LocalEllipticIpdgIntegrands::method_name<method>::value()
                                      + "_inner_integral_operator_" + space_name<SP>::value_wo_grid() + "_on_"
                                      + XT::Grid::bindings::layer_name<layer>::value() + "_"
                                      + XT::Grid::bindings::backend_name<SP::layer_backend>::value() + "_intersection";

      m.def(method_name.c_str(),
            [](const DF& diffusion_factor, const DT& diffusion_tensor, const size_t over_integrate) {
              return type(over_integrate, diffusion_factor, diffusion_tensor);
            },
            "diffusion_factor"_a,
            "diffusion_tensor"_a,
            "over_integrate"_a = 0,
            py::keep_alive<0, 1>(),
            py::keep_alive<0, 2>());
    } // ... addbind_factory_methods(...)
  }; // struct diffusion_switch

  struct diffusion_switch_scalar_base
  {
    static void addbind_factory_methods(pybind11::module& m)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      const std::string method_name = "make_local_elliptic_" + LocalEllipticIpdgIntegrands::method_name<method>::value()
                                      + "_inner_integral_operator_" + space_name<SP>::value_wo_grid() + "_on_"
                                      + XT::Grid::bindings::layer_name<layer>::value() + "_"
                                      + XT::Grid::bindings::backend_name<SP::layer_backend>::value() + "_intersection";

      m.def(method_name.c_str(),
            [](const DF& diffusion, const size_t over_integrate) { return type(over_integrate, diffusion); },
            "diffusion"_a,
            "over_integrate"_a = 0,
            py::keep_alive<0, 1>());
    } // ... addbind_factory_methods(...)
  }; // struct diffusion_switch_scalar_base

  template <bool anything>
  struct diffusion_switch<true, true, anything> : public diffusion_switch_scalar_base
  {

    static std::string suffix()
    {
      return "single_diffusion_factor";
    }
  };

  template <bool anything>
  struct diffusion_switch<true, false, anything> : public diffusion_switch_scalar_base
  {

    static std::string suffix()
    {
      return "single_diffusion_tensor";
    }
  };

public:
  static bound_type bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName =
        XT::Common::to_camel_case("local_elliptic_" + LocalEllipticIpdgIntegrands::method_name<method>::value() + "_"
                                  + diffusion_switch<>::suffix()
                                  + "_inner_integral_operator_"
                                  + space_name<SP>::value()
                                  + "_on_"
                                  + XT::Grid::bindings::layer_name<layer>::value()
                                  + "_intersection");

    bound_type c(m, ClassName.c_str());

    diffusion_switch<>::addbind_factory_methods(m);

    return c;
  } // ... bind(...)
}; // class LocalEllipticIpdgInnerIntegralOperator


} // namespace bindings
} // namespace GDT
} // namespace Dune


// begin: this is what we need for the .so

#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_1D(                                                               \
    _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, _method)                    \
  Dune::GDT::bindings::                                                                                                \
      LocalEllipticIpdgInnerIntegralOperator<Dune::XT::Functions::                                                     \
                                                 LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<        \
                                                                                  typename Dune::XT::Grid::            \
                                                                                      Layer<_G,                        \
                                                                                            Dune::XT::Grid::Layers::   \
                                                                                                _g_layer,              \
                                                                                            Dune::XT::Grid::Backends:: \
                                                                                                _g_backend,            \
                                                                                            Dune::XT::Grid::DD::       \
                                                                                                SubdomainGrid<_G>>::   \
                                                                                          type>,                       \
                                                                              double,                                  \
                                                                              _d,                                      \
                                                                              double,                                  \
                                                                              1,                                       \
                                                                              1>,                                      \
                                             Dune::XT::Functions::                                                     \
                                                 LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<        \
                                                                                  typename Dune::XT::Grid::            \
                                                                                      Layer<_G,                        \
                                                                                            Dune::XT::Grid::Layers::   \
                                                                                                _g_layer,              \
                                                                                            Dune::XT::Grid::Backends:: \
                                                                                                _g_backend,            \
                                                                                            Dune::XT::Grid::DD::       \
                                                                                                SubdomainGrid<_G>>::   \
                                                                                          type>,                       \
                                                                              double,                                  \
                                                                              _d,                                      \
                                                                              double,                                  \
                                                                              _d,                                      \
                                                                              _d>,                                     \
                                             Dune::GDT::SpaceProvider<_G,                                              \
                                                                      Dune::XT::Grid::Layers::_g_layer,                \
                                                                      Dune::GDT::SpaceType::_s_type,                   \
                                                                      Dune::GDT::Backends::_s_backend,                 \
                                                                      _p,                                              \
                                                                      _R,                                              \
                                                                      _r,                                              \
                                                                      _rC>,                                            \
                                             Dune::XT::Grid::Layers::_i_layer_type,                                    \
                                             Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method>::bind(_m);       \
  Dune::GDT::bindings::                                                                                                \
      LocalEllipticIpdgInnerIntegralOperator<Dune::XT::Functions::                                                     \
                                                 LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<        \
                                                                                  typename Dune::XT::Grid::            \
                                                                                      Layer<_G,                        \
                                                                                            Dune::XT::Grid::Layers::   \
                                                                                                _g_layer,              \
                                                                                            Dune::XT::Grid::Backends:: \
                                                                                                _g_backend,            \
                                                                                            Dune::XT::Grid::DD::       \
                                                                                                SubdomainGrid<_G>>::   \
                                                                                          type>,                       \
                                                                              double,                                  \
                                                                              _d,                                      \
                                                                              double,                                  \
                                                                              1,                                       \
                                                                              1>,                                      \
                                             void,                                                                     \
                                             Dune::GDT::SpaceProvider<_G,                                              \
                                                                      Dune::XT::Grid::Layers::_g_layer,                \
                                                                      Dune::GDT::SpaceType::_s_type,                   \
                                                                      Dune::GDT::Backends::_s_backend,                 \
                                                                      _p,                                              \
                                                                      _R,                                              \
                                                                      _r,                                              \
                                                                      _rC>,                                            \
                                             Dune::XT::Grid::Layers::_i_layer_type,                                    \
                                             Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method>::bind(_m)

#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_DD(                                                               \
    _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, _method)                    \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_1D(                                                                     \
      _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, _method);                 \
  Dune::GDT::bindings::                                                                                                \
      LocalEllipticIpdgInnerIntegralOperator<Dune::XT::Functions::                                                     \
                                                 LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<        \
                                                                                  typename Dune::XT::Grid::            \
                                                                                      Layer<_G,                        \
                                                                                            Dune::XT::Grid::Layers::   \
                                                                                                _g_layer,              \
                                                                                            Dune::XT::Grid::Backends:: \
                                                                                                _g_backend,            \
                                                                                            Dune::XT::Grid::DD::       \
                                                                                                SubdomainGrid<_G>>::   \
                                                                                          type>,                       \
                                                                              double,                                  \
                                                                              _d,                                      \
                                                                              double,                                  \
                                                                              _d,                                      \
                                                                              _d>,                                     \
                                             void,                                                                     \
                                             Dune::GDT::SpaceProvider<_G,                                              \
                                                                      Dune::XT::Grid::Layers::_g_layer,                \
                                                                      Dune::GDT::SpaceType::_s_type,                   \
                                                                      Dune::GDT::Backends::_s_backend,                 \
                                                                      _p,                                              \
                                                                      _R,                                              \
                                                                      _r,                                              \
                                                                      _rC>,                                            \
                                             Dune::XT::Grid::Layers::_i_layer_type,                                    \
                                             Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method>::bind(_m);

#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_METHODS_1D(                                                       \
    _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC)                             \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_1D(                                                                     \
      _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, sipdg);                   \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_1D(                                                                     \
      _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, swipdg);                  \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_1D(                                                                     \
      _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, swipdg_affine_factor);    \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_1D(                                                                     \
      _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, swipdg_affine_tensor)

#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_METHODS_DD(                                                       \
    _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC)                             \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_DD(                                                                     \
      _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, sipdg);                   \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_DD(                                                                     \
      _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, swipdg);                  \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_DD(                                                                     \
      _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, swipdg_affine_factor);    \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_DD(                                                                     \
      _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, swipdg_affine_tensor)

#if HAVE_DUNE_ALUGRID
#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_ALU(                                                              \
    _m, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC)                                     \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_METHODS_DD(                                                             \
      _m, ALU_2D_SIMPLEX_CONFORMING, 2, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC)
#else
#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_ALU(                                                              \
    _m, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC)
#endif

#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_YASP(                                                             \
    _m, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC)                                     \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_METHODS_1D(                                                             \
      _m, YASP_1D_EQUIDISTANT_OFFSET, 1, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC);   \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_METHODS_DD(                                                             \
      _m, YASP_2D_EQUIDISTANT_OFFSET, 2, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC)

#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_ALL_GRIDS(                                                        \
    _m, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC)                                     \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_ALU(                                                                    \
      _m, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC);                                  \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_YASP(                                                                   \
      _m, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC)

#if HAVE_DUNE_FEM
#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_PARTS(_m, _s_type, _s_backend, _p, _R, _r, _rC)                   \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_ALL_GRIDS(                                                              \
      _m, adaptive_leaf, part, _s_type, _s_backend, adaptive_leaf, _p, _R, _r, _rC);                                   \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_ALL_GRIDS(_m, leaf, part, _s_type, _s_backend, leaf, _p, _R, _r, _rC);  \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_ALL_GRIDS(                                                              \
      _m, level, part, _s_type, _s_backend, level, _p, _R, _r, _rC);                                                   \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_ALL_GRIDS(                                                              \
      _m, dd_subdomain, part, _s_type, _s_backend, dd_subdomain, _p, _R, _r, _rC);                                     \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_ALL_GRIDS(                                                              \
      _m, dd_subdomain, part, _s_type, _s_backend, dd_subdomain_coupling, _p, _R, _r, _rC);                            \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_ALL_GRIDS(                                                              \
      _m, dd_subdomain, part, _s_type, _s_backend, dd_subdomain_boundary, _p, _R, _r, _rC)
#else
#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_PARTS(_m, _s_type, _s_backend, _p, _R, _r, _rC)
#endif

//#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_VIEWS(_m, _s_type, _s_backend, _p, _R, _r, _rC)                   \
//  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_ALL_GRIDS(_m, leaf, view, _s_type, _s_backend, leaf, _p, _R, _r, _rC);  \
//  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_ALL_GRIDS(_m, level, view, _s_type, _s_backend, level, _p, _R, _r, _rC)

#if HAVE_DUNE_FEM
#define DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_FEM(_m)                                                            \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_PARTS(_m, cg, fem, 1, double, 1, 1);                                    \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_PARTS(_m, dg, fem, 1, double, 1, 1)
#else
#define DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_FEM(_m)
#endif

#define DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND(_m) DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_FEM(_m)

// end: this is what we need for the .so

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BINDINGS_HH
