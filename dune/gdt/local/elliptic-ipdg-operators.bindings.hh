// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
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
#include <dune/gdt/type_traits.hh>

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
  typedef XT::Grid::extract_intersection_t<
      typename XT::Grid::Layer<G, layer, S::layer_backend, XT::Grid::DD::SubdomainGrid<G>>::type>
      I;

public:
  typedef GDT::LocalCouplingTwoFormInterface<B, I> InterfaceType;
  typedef LocalCouplingIntegralOperator<LocalEllipticIpdgIntegrands::Inner<DF, DT, method>, B, I> type;
  typedef pybind11::class_<type, InterfaceType> bound_type;

private:
  template <bool intersection_matches_layer = (layer == SP::grid_layer), bool anything = false>
  struct intersection_postfix
  {
    static std::string value()
    {
      return "";
    }
  };

  template <bool anything>
  struct intersection_postfix<false, anything>
  {
    static std::string value()
    {
      return "_" + XT::Grid::bindings::layer_name<layer>::value() + "_intersection";
    }
  };

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
      using XT::Common::to_string;

      const std::string method_name =
          "make_local_elliptic_" + LocalEllipticIpdgIntegrands::method_name<method>::value()
          + "_inner_integral_operator_" + to_string(size_t(S::dimRange)) + "x" + to_string(size_t(S::dimRangeCols))
          + "_p" + to_string(int(S::polOrder)) + "_" + space_type_name<SP::space_type>::value() + "_"
          + backend_name<SP::space_backend>::value() + "_space" + intersection_postfix<>::value();

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
      using XT::Common::to_string;

      const std::string method_name =
          "make_local_elliptic_" + LocalEllipticIpdgIntegrands::method_name<method>::value()
          + "_inner_integral_operator_" + to_string(size_t(S::dimRange)) + "x" + to_string(size_t(S::dimRangeCols))
          + "_p" + to_string(int(S::polOrder)) + "_" + space_type_name<SP::space_type>::value() + "_"
          + backend_name<SP::space_backend>::value() + "_space" + intersection_postfix<>::value();

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
  static bound_type bind(pybind11::module& m, const std::string& layer_name)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;
    using XT::Common::to_string;

    // bind interface, guard since we might not be the first to do so for this intersection
    try {
      const auto InterfaceName = XT::Common::to_camel_case(
          "local_coupling_two_form_interface_" + XT::Grid::bindings::grid_name<G>::value() + layer_name + "_to_"
          + to_string(size_t(S::dimRange))
          + "x"
          + to_string(size_t(S::dimRangeCols))
          + "_p"
          + to_string(int(S::polOrder))
          + "_"
          + space_type_name<SP::space_type>::value()
          + "_"
          + backend_name<SP::space_backend>::value()
          + "_space"
          + intersection_postfix<>::value());
      py::class_<InterfaceType>(m, InterfaceName.c_str(), InterfaceName.c_str(), py::metaclass());
    } catch (std::runtime_error&) {
    }

    const auto ClassName =
        XT::Common::to_camel_case("local_elliptic_" + LocalEllipticIpdgIntegrands::method_name<method>::value() + "_"
                                  + diffusion_switch<>::suffix()
                                  + "_inner_integral_operator_"
                                  + XT::Grid::bindings::grid_name<G>::value()
                                  + layer_name
                                  + "_to_"
                                  + to_string(size_t(S::dimRange))
                                  + "x"
                                  + to_string(size_t(S::dimRangeCols))
                                  + "_p"
                                  + to_string(int(S::polOrder))
                                  + "_"
                                  + space_type_name<SP::space_type>::value()
                                  + "_"
                                  + backend_name<SP::space_backend>::value()
                                  + "_space"
                                  + intersection_postfix<>::value());

    bound_type c(m, ClassName.c_str());

    diffusion_switch<>::addbind_factory_methods(m);

    return c;
  } // ... bind(...)
}; // class LocalEllipticIpdgInnerIntegralOperator


template <class DF, typename DT, class SP, XT::Grid::Layers layer, LocalEllipticIpdgIntegrands::Method method>
class LocalEllipticIpdgBoundaryIntegralOperator
{
  static_assert(XT::Functions::is_localizable_function<DF>::value, "");
  static_assert(XT::Functions::is_localizable_function<DT>::value || std::is_same<DT, void>::value, "");
  typedef typename SP::type S;
  static_assert(is_space<S>::value, "");
  typedef XT::Grid::extract_grid_t<typename S::GridLayerType> G;
  typedef typename S::BaseFunctionSetType B;
  typedef XT::Grid::extract_intersection_t<
      typename XT::Grid::Layer<G, layer, S::layer_backend, XT::Grid::DD::SubdomainGrid<G>>::type>
      I;

public:
  typedef GDT::LocalBoundaryTwoFormInterface<B, I> InterfaceType;
  typedef LocalBoundaryIntegralOperator<LocalEllipticIpdgIntegrands::BoundaryLHS<DF, DT, method>, B, I> type;
  typedef pybind11::class_<type, InterfaceType> bound_type;

private:
  template <bool intersection_matches_layer = (layer == SP::grid_layer), bool anything = false>
  struct intersection_postfix
  {
    static std::string value()
    {
      return "";
    }
  };

  template <bool anything>
  struct intersection_postfix<false, anything>
  {
    static std::string value()
    {
      return "_" + XT::Grid::bindings::layer_name<layer>::value() + "_intersection";
    }
  };

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
      using XT::Common::to_string;

      const std::string method_name =
          "make_local_elliptic_" + LocalEllipticIpdgIntegrands::method_name<method>::value()
          + "_boundary_integral_operator_" + to_string(size_t(S::dimRange)) + "x" + to_string(size_t(S::dimRangeCols))
          + "_p" + to_string(int(S::polOrder)) + "_" + space_type_name<SP::space_type>::value() + "_"
          + backend_name<SP::space_backend>::value() + "_space" + intersection_postfix<>::value();

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
      using XT::Common::to_string;

      const std::string method_name =
          "make_local_elliptic_" + LocalEllipticIpdgIntegrands::method_name<method>::value()
          + "_boundary_integral_operator_" + to_string(size_t(S::dimRange)) + "x" + to_string(size_t(S::dimRangeCols))
          + "_p" + to_string(int(S::polOrder)) + "_" + space_type_name<SP::space_type>::value() + "_"
          + backend_name<SP::space_backend>::value() + "_space" + intersection_postfix<>::value();

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
  static bound_type bind(pybind11::module& m, const std::string& layer_name)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;
    using XT::Common::to_string;

    // bind interface, guard since we might not be the first to do so for this intersection
    try {
      const auto InterfaceName = XT::Common::to_camel_case(
          "local_boundary_two_form_interface_" + XT::Grid::bindings::grid_name<G>::value() + layer_name + "_to_"
          + to_string(size_t(S::dimRange))
          + "x"
          + to_string(size_t(S::dimRangeCols))
          + "_p"
          + to_string(int(S::polOrder))
          + "_"
          + space_type_name<SP::space_type>::value()
          + "_"
          + backend_name<SP::space_backend>::value()
          + "_space"
          + intersection_postfix<>::value());
      py::class_<InterfaceType>(m, InterfaceName.c_str(), InterfaceName.c_str(), py::metaclass());
    } catch (std::runtime_error&) {
    }

    const auto ClassName =
        XT::Common::to_camel_case("local_elliptic_" + LocalEllipticIpdgIntegrands::method_name<method>::value() + "_"
                                  + diffusion_switch<>::suffix()
                                  + "_boundary_integral_operator_"
                                  + XT::Grid::bindings::grid_name<G>::value()
                                  + layer_name
                                  + "_to_"
                                  + to_string(size_t(S::dimRange))
                                  + "x"
                                  + to_string(size_t(S::dimRangeCols))
                                  + "_p"
                                  + to_string(int(S::polOrder))
                                  + "_"
                                  + space_type_name<SP::space_type>::value()
                                  + "_"
                                  + backend_name<SP::space_backend>::value()
                                  + "_space"
                                  + intersection_postfix<>::value());

    bound_type c(m, ClassName.c_str());

    diffusion_switch<>::addbind_factory_methods(m);

    return c;
  } // ... bind(...)
}; // class LocalEllipticIpdgBoundaryIntegralOperator


} // namespace bindings
} // namespace GDT
} // namespace Dune


// begin: this is what we need for the .so

#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_1D(                                                                \
    _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, _layer_name, _method)        \
  Dune::GDT::bindings::                                                                                                 \
      LocalEllipticIpdgInnerIntegralOperator<Dune::XT::Functions::                                                      \
                                                 LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                  typename Dune::XT::Grid::             \
                                                                                      Layer<_G,                         \
                                                                                            Dune::XT::Grid::Layers::    \
                                                                                                _g_layer,               \
                                                                                            Dune::XT::Grid::Backends::  \
                                                                                                _g_backend,             \
                                                                                            Dune::XT::Grid::DD::        \
                                                                                                SubdomainGrid<_G>>::    \
                                                                                          type>,                        \
                                                                              double,                                   \
                                                                              _d,                                       \
                                                                              double,                                   \
                                                                              1,                                        \
                                                                              1>,                                       \
                                             Dune::XT::Functions::                                                      \
                                                 LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                  typename Dune::XT::Grid::             \
                                                                                      Layer<_G,                         \
                                                                                            Dune::XT::Grid::Layers::    \
                                                                                                _g_layer,               \
                                                                                            Dune::XT::Grid::Backends::  \
                                                                                                _g_backend,             \
                                                                                            Dune::XT::Grid::DD::        \
                                                                                                SubdomainGrid<_G>>::    \
                                                                                          type>,                        \
                                                                              double,                                   \
                                                                              _d,                                       \
                                                                              double,                                   \
                                                                              _d,                                       \
                                                                              _d>,                                      \
                                             Dune::GDT::SpaceProvider<_G,                                               \
                                                                      Dune::XT::Grid::Layers::_g_layer,                 \
                                                                      Dune::GDT::SpaceType::_s_type,                    \
                                                                      Dune::GDT::Backends::_s_backend,                  \
                                                                      _p,                                               \
                                                                      _R,                                               \
                                                                      _r,                                               \
                                                                      _rC>,                                             \
                                             Dune::XT::Grid::Layers::_i_layer_type,                                     \
                                             Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method>::                 \
          bind(_m, _layer_name);                                                                                        \
  Dune::GDT::bindings::                                                                                                 \
      LocalEllipticIpdgInnerIntegralOperator<Dune::XT::Functions::                                                      \
                                                 LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                  typename Dune::XT::Grid::             \
                                                                                      Layer<_G,                         \
                                                                                            Dune::XT::Grid::Layers::    \
                                                                                                _g_layer,               \
                                                                                            Dune::XT::Grid::Backends::  \
                                                                                                _g_backend,             \
                                                                                            Dune::XT::Grid::DD::        \
                                                                                                SubdomainGrid<_G>>::    \
                                                                                          type>,                        \
                                                                              double,                                   \
                                                                              _d,                                       \
                                                                              double,                                   \
                                                                              1,                                        \
                                                                              1>,                                       \
                                             void,                                                                      \
                                             Dune::GDT::SpaceProvider<_G,                                               \
                                                                      Dune::XT::Grid::Layers::_g_layer,                 \
                                                                      Dune::GDT::SpaceType::_s_type,                    \
                                                                      Dune::GDT::Backends::_s_backend,                  \
                                                                      _p,                                               \
                                                                      _R,                                               \
                                                                      _r,                                               \
                                                                      _rC>,                                             \
                                             Dune::XT::Grid::Layers::_i_layer_type,                                     \
                                             Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method>::                 \
          bind(_m, _layer_name);                                                                                        \
  Dune::GDT::bindings::                                                                                                 \
      LocalEllipticIpdgBoundaryIntegralOperator<Dune::XT::Functions::                                                   \
                                                    LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<      \
                                                                                     typename Dune::XT::Grid::          \
                                                                                         Layer<_G,                      \
                                                                                               Dune::XT::Grid::         \
                                                                                                   Layers::_g_layer,    \
                                                                                               Dune::XT::Grid::         \
                                                                                                   Backends::           \
                                                                                                       _g_backend,      \
                                                                                               Dune::XT::Grid::DD::     \
                                                                                                   SubdomainGrid<_G>>:: \
                                                                                             type>,                     \
                                                                                 double,                                \
                                                                                 _d,                                    \
                                                                                 double,                                \
                                                                                 1,                                     \
                                                                                 1>,                                    \
                                                Dune::XT::Functions::                                                   \
                                                    LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<      \
                                                                                     typename Dune::XT::Grid::          \
                                                                                         Layer<_G,                      \
                                                                                               Dune::XT::Grid::         \
                                                                                                   Layers::_g_layer,    \
                                                                                               Dune::XT::Grid::         \
                                                                                                   Backends::           \
                                                                                                       _g_backend,      \
                                                                                               Dune::XT::Grid::DD::     \
                                                                                                   SubdomainGrid<_G>>:: \
                                                                                             type>,                     \
                                                                                 double,                                \
                                                                                 _d,                                    \
                                                                                 double,                                \
                                                                                 _d,                                    \
                                                                                 _d>,                                   \
                                                Dune::GDT::SpaceProvider<_G,                                            \
                                                                         Dune::XT::Grid::Layers::_g_layer,              \
                                                                         Dune::GDT::SpaceType::_s_type,                 \
                                                                         Dune::GDT::Backends::_s_backend,               \
                                                                         _p,                                            \
                                                                         _R,                                            \
                                                                         _r,                                            \
                                                                         _rC>,                                          \
                                                Dune::XT::Grid::Layers::_i_layer_type,                                  \
                                                Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method>::              \
          bind(_m, _layer_name);                                                                                        \
  Dune::GDT::bindings::                                                                                                 \
      LocalEllipticIpdgBoundaryIntegralOperator<Dune::XT::Functions::                                                   \
                                                    LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<      \
                                                                                     typename Dune::XT::Grid::          \
                                                                                         Layer<_G,                      \
                                                                                               Dune::XT::Grid::         \
                                                                                                   Layers::_g_layer,    \
                                                                                               Dune::XT::Grid::         \
                                                                                                   Backends::           \
                                                                                                       _g_backend,      \
                                                                                               Dune::XT::Grid::DD::     \
                                                                                                   SubdomainGrid<_G>>:: \
                                                                                             type>,                     \
                                                                                 double,                                \
                                                                                 _d,                                    \
                                                                                 double,                                \
                                                                                 1,                                     \
                                                                                 1>,                                    \
                                                void,                                                                   \
                                                Dune::GDT::SpaceProvider<_G,                                            \
                                                                         Dune::XT::Grid::Layers::_g_layer,              \
                                                                         Dune::GDT::SpaceType::_s_type,                 \
                                                                         Dune::GDT::Backends::_s_backend,               \
                                                                         _p,                                            \
                                                                         _R,                                            \
                                                                         _r,                                            \
                                                                         _rC>,                                          \
                                                Dune::XT::Grid::Layers::_i_layer_type,                                  \
                                                Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method>::              \
          bind(_m, _layer_name)

#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_DD(                                                                \
    _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, _layer_name, _method)        \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_1D(                                                                      \
      _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, _layer_name, _method);     \
  Dune::GDT::bindings::                                                                                                 \
      LocalEllipticIpdgInnerIntegralOperator<Dune::XT::Functions::                                                      \
                                                 LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                  typename Dune::XT::Grid::             \
                                                                                      Layer<_G,                         \
                                                                                            Dune::XT::Grid::Layers::    \
                                                                                                _g_layer,               \
                                                                                            Dune::XT::Grid::Backends::  \
                                                                                                _g_backend,             \
                                                                                            Dune::XT::Grid::DD::        \
                                                                                                SubdomainGrid<_G>>::    \
                                                                                          type>,                        \
                                                                              double,                                   \
                                                                              _d,                                       \
                                                                              double,                                   \
                                                                              _d,                                       \
                                                                              _d>,                                      \
                                             void,                                                                      \
                                             Dune::GDT::SpaceProvider<_G,                                               \
                                                                      Dune::XT::Grid::Layers::_g_layer,                 \
                                                                      Dune::GDT::SpaceType::_s_type,                    \
                                                                      Dune::GDT::Backends::_s_backend,                  \
                                                                      _p,                                               \
                                                                      _R,                                               \
                                                                      _r,                                               \
                                                                      _rC>,                                             \
                                             Dune::XT::Grid::Layers::_i_layer_type,                                     \
                                             Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method>::                 \
          bind(_m, _layer_name);                                                                                        \
  Dune::GDT::bindings::                                                                                                 \
      LocalEllipticIpdgBoundaryIntegralOperator<Dune::XT::Functions::                                                   \
                                                    LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<      \
                                                                                     typename Dune::XT::Grid::          \
                                                                                         Layer<_G,                      \
                                                                                               Dune::XT::Grid::         \
                                                                                                   Layers::_g_layer,    \
                                                                                               Dune::XT::Grid::         \
                                                                                                   Backends::           \
                                                                                                       _g_backend,      \
                                                                                               Dune::XT::Grid::DD::     \
                                                                                                   SubdomainGrid<_G>>:: \
                                                                                             type>,                     \
                                                                                 double,                                \
                                                                                 _d,                                    \
                                                                                 double,                                \
                                                                                 _d,                                    \
                                                                                 _d>,                                   \
                                                void,                                                                   \
                                                Dune::GDT::SpaceProvider<_G,                                            \
                                                                         Dune::XT::Grid::Layers::_g_layer,              \
                                                                         Dune::GDT::SpaceType::_s_type,                 \
                                                                         Dune::GDT::Backends::_s_backend,               \
                                                                         _p,                                            \
                                                                         _R,                                            \
                                                                         _r,                                            \
                                                                         _rC>,                                          \
                                                Dune::XT::Grid::Layers::_i_layer_type,                                  \
                                                Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method>::              \
          bind(_m, _layer_name)

#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_METHODS_1D(                                                       \
    _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, _layer_name)                \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_1D(                                                                     \
      _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, _layer_name, sipdg);      \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_1D(                                                                     \
      _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, _layer_name, swipdg);     \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_1D(_m,                                                                  \
                                                  _G,                                                                  \
                                                  _d,                                                                  \
                                                  _g_layer,                                                            \
                                                  _g_backend,                                                          \
                                                  _s_type,                                                             \
                                                  _s_backend,                                                          \
                                                  _i_layer_type,                                                       \
                                                  _p,                                                                  \
                                                  _R,                                                                  \
                                                  _r,                                                                  \
                                                  _rC,                                                                 \
                                                  _layer_name,                                                         \
                                                  swipdg_affine_factor);                                               \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_1D(_m,                                                                  \
                                                  _G,                                                                  \
                                                  _d,                                                                  \
                                                  _g_layer,                                                            \
                                                  _g_backend,                                                          \
                                                  _s_type,                                                             \
                                                  _s_backend,                                                          \
                                                  _i_layer_type,                                                       \
                                                  _p,                                                                  \
                                                  _R,                                                                  \
                                                  _r,                                                                  \
                                                  _rC,                                                                 \
                                                  _layer_name,                                                         \
                                                  swipdg_affine_tensor)

#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_METHODS_DD(                                                       \
    _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, _layer_name)                \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_DD(                                                                     \
      _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, _layer_name, sipdg);      \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_DD(                                                                     \
      _m, _G, _d, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, _layer_name, swipdg);     \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_DD(_m,                                                                  \
                                                  _G,                                                                  \
                                                  _d,                                                                  \
                                                  _g_layer,                                                            \
                                                  _g_backend,                                                          \
                                                  _s_type,                                                             \
                                                  _s_backend,                                                          \
                                                  _i_layer_type,                                                       \
                                                  _p,                                                                  \
                                                  _R,                                                                  \
                                                  _r,                                                                  \
                                                  _rC,                                                                 \
                                                  _layer_name,                                                         \
                                                  swipdg_affine_factor);                                               \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_DD(_m,                                                                  \
                                                  _G,                                                                  \
                                                  _d,                                                                  \
                                                  _g_layer,                                                            \
                                                  _g_backend,                                                          \
                                                  _s_type,                                                             \
                                                  _s_backend,                                                          \
                                                  _i_layer_type,                                                       \
                                                  _p,                                                                  \
                                                  _R,                                                                  \
                                                  _r,                                                                  \
                                                  _rC,                                                                 \
                                                  _layer_name,                                                         \
                                                  swipdg_affine_tensor)

#if HAVE_DUNE_ALUGRID
#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_ALU(                                                              \
    _m, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, _layer_name)                        \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_METHODS_DD(_m,                                                          \
                                                          ALU_2D_SIMPLEX_CONFORMING,                                   \
                                                          2,                                                           \
                                                          _g_layer,                                                    \
                                                          _g_backend,                                                  \
                                                          _s_type,                                                     \
                                                          _s_backend,                                                  \
                                                          _i_layer_type,                                               \
                                                          _p,                                                          \
                                                          _R,                                                          \
                                                          _r,                                                          \
                                                          _rC,                                                         \
                                                          _layer_name)
#else
#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_ALU(                                                              \
    _m, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, _layer_name)
#endif

#define _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_YASP(                                                             \
    _m, _g_layer, _g_backend, _s_type, _s_backend, _i_layer_type, _p, _R, _r, _rC, _layer_name)                        \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_METHODS_1D(_m,                                                          \
                                                          YASP_1D_EQUIDISTANT_OFFSET,                                  \
                                                          1,                                                           \
                                                          _g_layer,                                                    \
                                                          _g_backend,                                                  \
                                                          _s_type,                                                     \
                                                          _s_backend,                                                  \
                                                          _i_layer_type,                                               \
                                                          _p,                                                          \
                                                          _R,                                                          \
                                                          _r,                                                          \
                                                          _rC,                                                         \
                                                          _layer_name);                                                \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_METHODS_DD(_m,                                                          \
                                                          YASP_2D_EQUIDISTANT_OFFSET,                                  \
                                                          2,                                                           \
                                                          _g_layer,                                                    \
                                                          _g_backend,                                                  \
                                                          _s_type,                                                     \
                                                          _s_backend,                                                  \
                                                          _i_layer_type,                                               \
                                                          _p,                                                          \
                                                          _R,                                                          \
                                                          _r,                                                          \
                                                          _rC,                                                         \
                                                          _layer_name);                                                \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_METHODS_DD(_m,                                                          \
                                                          YASP_3D_EQUIDISTANT_OFFSET,                                  \
                                                          3,                                                           \
                                                          _g_layer,                                                    \
                                                          _g_backend,                                                  \
                                                          _s_type,                                                     \
                                                          _s_backend,                                                  \
                                                          _i_layer_type,                                               \
                                                          _p,                                                          \
                                                          _R,                                                          \
                                                          _r,                                                          \
                                                          _rC,                                                         \
                                                          _layer_name)

//_DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_ALU(                                                                    \
//    _m, leaf, part, _s_type, _s_backend, leaf, _p, _R, _r, _rC, "_leaf_");                                           \
//_DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_ALU(                                                                    \
//    _m, level, part, _s_type, _s_backend, level, _p, _R, _r, _rC, "_level_");                                        \
//_DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_YASP(_m, leaf, part, _s_type, _s_backend, leaf, _p, _R, _r, _rC, "");   \
//_DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_YASP(                                                                   \
//    _m, dd_subdomain, part, _s_type, _s_backend, dd_subdomain_coupling, _p, _R, _r, _rC, "_dd_subdomain_")

#define DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND(_m)                                                                \
  _DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BIND_ALU(                                                                    \
      _m, dd_subdomain, view, dg, gdt, dd_subdomain_coupling, 1, double, 1, 1, "_dd_subdomain_")

// end: this is what we need for the .so

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_LOCAL_ELLIPTIC_IPDG_OPERATORS_BINDINGS_HH
