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

#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>

#include "integrands/elliptic-ipdg.hh"
#include "operators/integrals.hh"
#include "assembler.bindings.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <class G, class DF, typename DT, LocalEllipticIpdgIntegrands::Method method>
class LocalEllipticIpdgInnerIntegralOperator
{
  static_assert(XT::Functions::is_localizable_function<DF>::value, "");
  static_assert(XT::Functions::is_localizable_function<DT>::value || std::is_same<DT, void>::value, "");

public:
  typedef LocalCouplingIntegralOperator<LocalEllipticIpdgIntegrands::Inner<DF, DT, method>> type;
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
                                      + "_inner_integral_operator";

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
                                      + "_inner_integral_operator";

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
                                  + "_inner_integral_operator"
                                  + XT::Grid::bindings::grid_name<G>::value());

    bound_type c(m, ClassName.c_str());

    diffusion_switch<>::addbind_factory_methods(m);

    // the matching assembler
    bindings::LocalCouplingTwoFormAssembler<type>::bind(m, ClassName);

    return c;
  } // ... bind(...)
}; // class LocalEllipticIpdgInnerIntegralOperator


} // namespace bindings
} // namespace GDT
} // namespace Dune


// begin: this is what we need for the .so

#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_1D(_m, _G, _d, _g_layer, _g_backend, _method)                     \
  Dune::GDT::bindings::                                                                                                \
      LocalEllipticIpdgInnerIntegralOperator<_G,                                                                       \
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
                                             Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method>::bind(_m);       \
  Dune::GDT::bindings::                                                                                                \
      LocalEllipticIpdgInnerIntegralOperator<_G,                                                                       \
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
                                                                              1,                                       \
                                                                              1>,                                      \
                                             void,                                                                     \
                                             Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method>::bind(_m)

#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_DD(_m, _G, _d, _g_layer, _g_backend, _method)                     \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_1D(_m, _G, _d, _g_layer, _g_backend, _method);                          \
  Dune::GDT::bindings::                                                                                                \
      LocalEllipticIpdgInnerIntegralOperator<_G,                                                                       \
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
                                             void,                                                                     \
                                             Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method>::bind(_m);

#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_METHODS_1D(_m, _G, _d, _g_layer, _g_backend)                      \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_1D(_m, _G, _d, _g_layer, _g_backend, sipdg);                            \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_1D(_m, _G, _d, _g_layer, _g_backend, swipdg);                           \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_1D(_m, _G, _d, _g_layer, _g_backend, swipdg_affine_factor);             \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_1D(_m, _G, _d, _g_layer, _g_backend, swipdg_affine_tensor)

#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_METHODS_DD(_m, _G, _d, _g_layer, _g_backend)                      \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_DD(_m, _G, _d, _g_layer, _g_backend, sipdg);                            \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_DD(_m, _G, _d, _g_layer, _g_backend, swipdg);                           \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_DD(_m, _G, _d, _g_layer, _g_backend, swipdg_affine_factor);             \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_DD(_m, _G, _d, _g_layer, _g_backend, swipdg_affine_tensor)

#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_YASP(_m)                                                          \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_METHODS_1D(_m, YASP_1D_EQUIDISTANT_OFFSET, 1, leaf, view);              \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_METHODS_DD(_m, YASP_2D_EQUIDISTANT_OFFSET, 2, leaf, view)

#if HAVE_DUNE_ALUGRID
#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_ALU(_m)                                                           \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_METHODS_DD(_m, ALU_2D_SIMPLEX_CONFORMING, 2, level, view)
#else
#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_ALU(_m)
#endif

#define DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND(_m)                                                                \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_YASP(_m);                                                               \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_ALU(_m)


// end: this is what we need for the .so

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BINDINGS_HH
