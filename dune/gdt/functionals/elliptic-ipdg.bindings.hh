// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BINDINGS_HH
#define DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/la/container.bindings.hh>

#include <dune/gdt/spaces/interface.bindings.hh>
#include <dune/gdt/type_traits.hh>

#include "elliptic-ipdg.hh"
#include "base.bindings.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <class DI,
          class DF,
          typename DT, // may be void
          class SP,
          LocalEllipticIpdgIntegrands::Method method,
          class V /* = typename XT::LA::Container<typename R::RangeFieldType>::VectorType,
          class GL = typename RP::type::GridLayerType,
          class F = typename RP::type::RangeFieldType*/>
class EllipticIpdgDirichletVectorFunctional
{
  typedef typename SP::type S;
  static_assert(is_space<S>::value, "");

public:
  typedef GDT::EllipticIpdgDirichletVectorFunctional<DI, DF, DT, S, method, V /*, GL, F*/> type;
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

    template <class C>
    static void addbind_factory_methods(pybind11::module& m)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      const auto method_name =
          "make_elliptic_" + LocalEllipticIpdgIntegrands::method_name<method>::value() + "_dirichlet_vector_functional";

      m.def(
          std::string(method_name + "_" + XT::LA::bindings::container_name<V>::value()).c_str(),
          [](const DI& dirichlet,
             const DF& diffusion_factor,
             const DT& diffusion_tensor,
             const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<typename S::GridLayerType>>& boundary_info,
             const S& space,
             const size_t over_integrate) {
            return make_elliptic_ipdg_dirichlet_vector_functional<V, method>(
                       dirichlet, diffusion_factor, diffusion_tensor, boundary_info, space, over_integrate)
                .release(); //   <- b.c. EllipticIpdgDirichletVectorFunctional is not movable, returning the raw pointer
          }, //                                                                lets pybind11 correctly manage the memory
          "dirichlet"_a,
          "diffusion_factor"_a,
          "diffusion_tensor"_a,
          "boundary_info"_a,
          "space"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>(),
          py::keep_alive<0, 4>(),
          py::keep_alive<0, 5>());

      m.def(
          method_name.c_str(),
          [](const DI& dirichlet,
             const DF& diffusion_factor,
             const DT& diffusion_tensor,
             const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<typename S::GridLayerType>>& boundary_info,
             V& vector,
             const S& space,
             const size_t over_integrate) {
            return make_elliptic_ipdg_dirichlet_vector_functional<method>(
                       dirichlet, diffusion_factor, diffusion_tensor, boundary_info, vector, space, over_integrate)
                .release(); //                                                                     <- s.a. for release()
          },
          "dirichlet"_a,
          "diffusion_factor"_a,
          "diffusion_tensor"_a,
          "boundary_info"_a,
          "vector"_a,
          "space"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>(),
          py::keep_alive<0, 4>(),
          py::keep_alive<0, 5>(),
          py::keep_alive<0, 6>());
    } // ... addbind_factory_methods(...)
  }; // struct diffusion_switch

  struct diffusion_switch_scalar_base
  {
    template <class C>
    static void addbind_factory_methods(pybind11::module& m)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      const auto method_name =
          "make_elliptic_" + LocalEllipticIpdgIntegrands::method_name<method>::value() + "_dirichlet_vector_functional";

      m.def(
          std::string(method_name + "_" + XT::LA::bindings::container_name<V>::value()).c_str(),
          [](const DI& dirichlet,
             const DF& diffusion,
             const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<typename S::GridLayerType>>& boundary_info,
             const S& space,
             const size_t over_integrate) {
            return make_elliptic_ipdg_dirichlet_vector_functional<V, method>(
                       dirichlet, diffusion, boundary_info, space, over_integrate)
                .release(); //   <- b.c. EllipticIpdgDirichletVectorFunctional is not movable, returning the raw pointer
          }, //                                                                lets pybind11 correctly manage the memory
          "dirichlet"_a,
          "diffusion"_a,
          "boundary_info"_a,
          "space"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>(),
          py::keep_alive<0, 4>());

      m.def(
          method_name.c_str(),
          [](const DI& dirichlet,
             const DF& diffusion,
             const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<typename S::GridLayerType>>& boundary_info,
             V& vector,
             const S& space,
             const size_t over_integrate) {
            return make_elliptic_ipdg_dirichlet_vector_functional<method>(
                       dirichlet, diffusion, boundary_info, vector, space, over_integrate)
                .release(); //                                                                     <- s.a. for release()
          },
          "dirichlet"_a,
          "diffusion"_a,
          "boundary_info"_a,
          "vector"_a,
          "space"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>(),
          py::keep_alive<0, 4>(),
          py::keep_alive<0, 5>());
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

    const auto ClassName = XT::Common::to_camel_case(
        "elliptic_" + LocalEllipticIpdgIntegrands::method_name<method>::value() + "_dirichlet_vector_functional_"
        + space_name<SP>::value()
        + "_"
        + XT::LA::bindings::container_name<V>::value()
        + "_"
        + diffusion_switch<>::suffix());

    auto c = VectorFunctionalBase<type>::bind(m, ClassName.c_str());

    diffusion_switch<>::template addbind_factory_methods<type>(m);

    return c;
  } // ... bind(...)
}; // class EllipticIpdgDirichletVectorFunctional


} // namespace bindings
} // namespace GDT
} // namespace Dune


// begin: this is what we need for the lib

#define _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_1D(                                                               \
    _prefix, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la, _method)                                         \
  _prefix class Dune::GDT::bindings::                                                                                  \
      EllipticIpdgDirichletVectorFunctional<Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             1>,                                       \
                                            Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             1>,                                       \
                                            Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             1>,                                       \
                                            Dune::GDT::SpaceProvider<_GRID,                                            \
                                                                     Dune::XT::Grid::Layers::_layer,                   \
                                                                     Dune::GDT::SpaceType::_s_type,                    \
                                                                     Dune::GDT::Backends::_s_backend,                  \
                                                                     _p,                                               \
                                                                     double,                                           \
                                                                     1,                                                \
                                                                     1>,                                               \
                                            Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method,                   \
                                            typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::    \
                                                VectorType>;                                                           \
  _prefix class Dune::GDT::bindings::                                                                                  \
      EllipticIpdgDirichletVectorFunctional<Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             1>,                                       \
                                            Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             1>,                                       \
                                            void,                                                                      \
                                            Dune::GDT::SpaceProvider<_GRID,                                            \
                                                                     Dune::XT::Grid::Layers::_layer,                   \
                                                                     Dune::GDT::SpaceType::_s_type,                    \
                                                                     Dune::GDT::Backends::_s_backend,                  \
                                                                     _p,                                               \
                                                                     double,                                           \
                                                                     1,                                                \
                                                                     1>,                                               \
                                            Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method,                   \
                                            typename Dune::XT::LA::Container<double,                                   \
                                                                             Dune::XT::LA::Backends::_la>::VectorType>

#define _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_D(                                                                \
    _prefix, _d, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la, _method)                                     \
  _prefix class Dune::GDT::bindings::                                                                                  \
      EllipticIpdgDirichletVectorFunctional<Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             _d,                                       \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             1>,                                       \
                                            Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             _d,                                       \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             1>,                                       \
                                            Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             _d,                                       \
                                                                             double,                                   \
                                                                             _d,                                       \
                                                                             _d>,                                      \
                                            Dune::GDT::SpaceProvider<_GRID,                                            \
                                                                     Dune::XT::Grid::Layers::_layer,                   \
                                                                     Dune::GDT::SpaceType::_s_type,                    \
                                                                     Dune::GDT::Backends::_s_backend,                  \
                                                                     _p,                                               \
                                                                     double,                                           \
                                                                     1,                                                \
                                                                     1>,                                               \
                                            Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method,                   \
                                            typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::    \
                                                VectorType>;                                                           \
  _prefix class Dune::GDT::bindings::                                                                                  \
      EllipticIpdgDirichletVectorFunctional<Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             _d,                                       \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             1>,                                       \
                                            Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             _d,                                       \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             1>,                                       \
                                            void,                                                                      \
                                            Dune::GDT::SpaceProvider<_GRID,                                            \
                                                                     Dune::XT::Grid::Layers::_layer,                   \
                                                                     Dune::GDT::SpaceType::_s_type,                    \
                                                                     Dune::GDT::Backends::_s_backend,                  \
                                                                     _p,                                               \
                                                                     double,                                           \
                                                                     1,                                                \
                                                                     1>,                                               \
                                            Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method,                   \
                                            typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::    \
                                                VectorType>;                                                           \
  _prefix class Dune::GDT::bindings::                                                                                  \
      EllipticIpdgDirichletVectorFunctional<Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             _d,                                       \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             1>,                                       \
                                            Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             _d,                                       \
                                                                             double,                                   \
                                                                             _d,                                       \
                                                                             _d>,                                      \
                                            void,                                                                      \
                                            Dune::GDT::SpaceProvider<_GRID,                                            \
                                                                     Dune::XT::Grid::Layers::_layer,                   \
                                                                     Dune::GDT::SpaceType::_s_type,                    \
                                                                     Dune::GDT::Backends::_s_backend,                  \
                                                                     _p,                                               \
                                                                     double,                                           \
                                                                     1,                                                \
                                                                     1>,                                               \
                                            Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method,                   \
                                            typename Dune::XT::LA::Container<double,                                   \
                                                                             Dune::XT::LA::Backends::_la>::VectorType>

#define _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_METHODS_1D(                                                       \
    _prefix, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la)                                                  \
  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_1D(                                                                     \
      _prefix, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la, sipdg);                                        \
  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_1D(                                                                     \
      _prefix, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la, swipdg);                                       \
  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_1D(                                                                     \
      _prefix, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la, swipdg_affine_factor);                         \
  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_1D(                                                                     \
      _prefix, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la, swipdg_affine_tensor)

#define _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_METHODS_D(                                                        \
    _prefix, _d, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la)                                              \
  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_D(                                                                      \
      _prefix, _d, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la, sipdg);                                    \
  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_D(                                                                      \
      _prefix, _d, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la, swipdg);                                   \
  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_D(                                                                      \
      _prefix, _d, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la, swipdg_affine_factor);                     \
  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_D(                                                                      \
      _prefix, _d, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la, swipdg_affine_tensor)

//#if HAVE_ALBERTA
//#define DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_ALBERTA(_prefix, _layer, _g_backend, _s_type, _s_backend, _p, _la) \
//  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_METHODS_D(                                                                \
//      _prefix, 2, ALBERTA_2D, _layer, _g_backend, _s_type, _s_backend, _p, _la)
//#else
#define DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_ALBERTA(_prefix, _layer, _g_backend, _s_type, _s_backend, _p, _la)
//#endif

#if HAVE_DUNE_ALUGRID
#define DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_ALU(_prefix, _layer, _g_backend, _s_type, _s_backend, _p, _la)     \
  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_METHODS_D(                                                              \
      _prefix, 2, ALU_2D_SIMPLEX_CONFORMING, _layer, _g_backend, _s_type, _s_backend, _p, _la)
#else
#define DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_ALU(_prefix, _layer, _g_backend, _s_type, _s_backend, _p, _la)
#endif

//#if HAVE_DUNE_UGGRID || HAVE_UG
//#define DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_UG(_prefix, _layer, _g_backend, _s_type, _s_backend, _p, _la)      \
//  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_METHODS_D(                                                                \
//      _prefix, 2, UG_2D, _layer, _g_backend, _s_type, _s_backend, _p, _la)
//#else
#define DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_UG(_prefix, _layer, _g_backend, _s_type, _s_backend, _p, _la)
//#endif

#define DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_YASP(_prefix, _layer, _g_backend, _s_type, _s_backend, _p, _la)    \
  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_METHODS_1D(                                                             \
      _prefix, YASP_1D_EQUIDISTANT_OFFSET, _layer, _g_backend, _s_type, _s_backend, _p, _la);                          \
  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_METHODS_D(                                                              \
      _prefix, 2, YASP_2D_EQUIDISTANT_OFFSET, _layer, _g_backend, _s_type, _s_backend, _p, _la)

// alu_istl.cc
#if HAVE_DUNE_ALUGRID && HAVE_DUNE_ISTL
DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_ALU(extern template, leaf, view, cg, gdt, 1, istl_dense);
DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_ALU(extern template, level, view, cg, gdt, 1, istl_dense);
// DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_ALU(extern template, dd_subdomain, part, cg, gdt, 1, istl_dense);
#endif

// yasp_istl.cc
#if HAVE_DUNE_ISTL
// DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_YASP(extern template, leaf, view, cg, gdt, 1, istl_dense);
// DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_YASP(extern template, level, view, cg, gdt, 1, istl_dense);
// DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_YASP(extern template, dd_subdomain, part, cg, gdt, 1, istl_dense);
//#endif

// end: this is what we need for the lib


// begin: this is what we need for the .so

#define _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_1D(                                                                   \
    _m, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la, _method)                                              \
  Dune::GDT::bindings::                                                                                                \
      EllipticIpdgDirichletVectorFunctional<Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             1>,                                       \
                                            Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             1>,                                       \
                                            Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             1>,                                       \
                                            Dune::GDT::SpaceProvider<_GRID,                                            \
                                                                     Dune::XT::Grid::Layers::_layer,                   \
                                                                     Dune::GDT::SpaceType::_s_type,                    \
                                                                     Dune::GDT::Backends::_s_backend,                  \
                                                                     _p,                                               \
                                                                     double,                                           \
                                                                     1,                                                \
                                                                     1>,                                               \
                                            Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method,                   \
                                            typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::    \
                                                VectorType>::bind(_m);                                                 \
  Dune::GDT::bindings::                                                                                                \
      EllipticIpdgDirichletVectorFunctional<Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             1>,                                       \
                                            Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             1>,                                       \
                                            void,                                                                      \
                                            Dune::GDT::SpaceProvider<_GRID,                                            \
                                                                     Dune::XT::Grid::Layers::_layer,                   \
                                                                     Dune::GDT::SpaceType::_s_type,                    \
                                                                     Dune::GDT::Backends::_s_backend,                  \
                                                                     _p,                                               \
                                                                     double,                                           \
                                                                     1,                                                \
                                                                     1>,                                               \
                                            Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method,                   \
                                            typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::    \
                                                VectorType>::bind(_m)

#define _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_D(                                                                    \
    _m, _d, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la, _method)                                          \
  Dune::GDT::bindings::                                                                                                \
      EllipticIpdgDirichletVectorFunctional<Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             _d,                                       \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             1>,                                       \
                                            Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             _d,                                       \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             1>,                                       \
                                            Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             _d,                                       \
                                                                             double,                                   \
                                                                             _d,                                       \
                                                                             _d>,                                      \
                                            Dune::GDT::SpaceProvider<_GRID,                                            \
                                                                     Dune::XT::Grid::Layers::_layer,                   \
                                                                     Dune::GDT::SpaceType::_s_type,                    \
                                                                     Dune::GDT::Backends::_s_backend,                  \
                                                                     _p,                                               \
                                                                     double,                                           \
                                                                     1,                                                \
                                                                     1>,                                               \
                                            Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method,                   \
                                            typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::    \
                                                VectorType>::bind(_m);                                                 \
  Dune::GDT::bindings::                                                                                                \
      EllipticIpdgDirichletVectorFunctional<Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             _d,                                       \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             1>,                                       \
                                            Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             _d,                                       \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             1>,                                       \
                                            void,                                                                      \
                                            Dune::GDT::SpaceProvider<_GRID,                                            \
                                                                     Dune::XT::Grid::Layers::_layer,                   \
                                                                     Dune::GDT::SpaceType::_s_type,                    \
                                                                     Dune::GDT::Backends::_s_backend,                  \
                                                                     _p,                                               \
                                                                     double,                                           \
                                                                     1,                                                \
                                                                     1>,                                               \
                                            Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method,                   \
                                            typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::    \
                                                VectorType>::bind(_m);                                                 \
  Dune::GDT::bindings::                                                                                                \
      EllipticIpdgDirichletVectorFunctional<Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             _d,                                       \
                                                                             double,                                   \
                                                                             1,                                        \
                                                                             1>,                                       \
                                            Dune::XT::Functions::                                                      \
                                                LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<         \
                                                                                 typename Dune::XT::Grid::             \
                                                                                     Layer<_GRID,                      \
                                                                                           Dune::XT::Grid::Layers::    \
                                                                                               _layer,                 \
                                                                                           Dune::XT::Grid::Backends::  \
                                                                                               _g_backend,             \
                                                                                           Dune::XT::Grid::DD::        \
                                                                                               SubdomainGrid<_GRID>>:: \
                                                                                         type>,                        \
                                                                             double,                                   \
                                                                             _d,                                       \
                                                                             double,                                   \
                                                                             _d,                                       \
                                                                             _d>,                                      \
                                            void,                                                                      \
                                            Dune::GDT::SpaceProvider<_GRID,                                            \
                                                                     Dune::XT::Grid::Layers::_layer,                   \
                                                                     Dune::GDT::SpaceType::_s_type,                    \
                                                                     Dune::GDT::Backends::_s_backend,                  \
                                                                     _p,                                               \
                                                                     double,                                           \
                                                                     1,                                                \
                                                                     1>,                                               \
                                            Dune::GDT::LocalEllipticIpdgIntegrands::Method::_method,                   \
                                            typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::    \
                                                VectorType>::bind(_m)

#define _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_METHODS_1D(                                                           \
    _m, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la)                                                       \
  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_1D(_m, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la, sipdg);     \
  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_1D(_m, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la, swipdg);    \
  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_1D(                                                                         \
      _m, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la, swipdg_affine_factor);                              \
  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_1D(                                                                         \
      _m, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la, swipdg_affine_tensor)

#define _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_METHODS_D(                                                            \
    _m, _d, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la)                                                   \
  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_D(                                                                          \
      _m, _d, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la, swipdg_affine_factor)
//_DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_D(_m, _d, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la, sipdg);  \
//_DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_D(_m, _d, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la, swipdg); \
//_DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_D(                                                                          \
//    _m, _d, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la, swipdg_affine_tensor)

//#if HAVE_ALBERTA
//#define DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_ALBERTA(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la)          \
//  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_METHODS_D(_m, 2, ALBERTA_2D, _layer, _g_backend, _s_type, _s_backend, _p, _la)
//#else
#define DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_ALBERTA(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la)
//#endif

#if HAVE_DUNE_ALUGRID
#define DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_ALU(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la)              \
  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_METHODS_D(                                                                  \
      _m, 2, ALU_2D_SIMPLEX_CONFORMING, _layer, _g_backend, _s_type, _s_backend, _p, _la)
#else
#define DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_ALU(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la)
#endif

//#if HAVE_DUNE_UGGRID || HAVE_UG
//#define DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_UG(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la)               \
//  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_METHODS_D(_m, 2, UG_2D, _layer, _g_backend, _s_type, _s_backend, _p, _la)
//#else
#define DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_UG(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la)
//#endif

#define DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_YASP(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la)             \
  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_METHODS_1D(                                                                 \
      _m, YASP_1D_EQUIDISTANT_OFFSET, _layer, _g_backend, _s_type, _s_backend, _p, _la);                               \
  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_METHODS_D(                                                                  \
      _m, 2, YASP_2D_EQUIDISTANT_OFFSET, _layer, _g_backend, _s_type, _s_backend, _p, _la)


// end: this is what we need for the .so

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BINDINGS_HH
