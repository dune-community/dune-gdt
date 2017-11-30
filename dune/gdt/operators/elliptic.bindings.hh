// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_OPERATORS_ELLIPTIC_BINDINGS_HH
#define DUNE_GDT_OPERATORS_ELLIPTIC_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/la/container.bindings.hh>

#include <dune/gdt/spaces.bindings.hh>
#include <dune/gdt/type_traits.hh>

#include "base.bindings.hh"
#include "elliptic.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <class DF,
          typename DT, // may be void
          class RP,
          class M /* = typename XT::LA::Container<typename R::RangeFieldType>::MatrixType,
          class GL = typename R::GridLayerType,
          class S = R,
          class F = typename R::RangeFieldType*/>
class EllipticMatrixOperator
{
  typedef typename RP::type R;
  static_assert(is_space<R>::value, "");

public:
  typedef GDT::EllipticMatrixOperator<DF, DT, R, M /*, GL, S, F*/> type;
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

      const std::string method_name = "make_elliptic_matrix_operator_" + XT::LA::bindings::container_name<M>::value();

      m.def(
          method_name.c_str(),
          [](const DF& diffusion_factor, const DT& diffusion_tensor, const R& space, const size_t over_integrate) {
            return make_elliptic_matrix_operator<M>(diffusion_factor, diffusion_tensor, space, over_integrate)
                .release(); //    <- b.c. EllipticMatrixOperator is not movable, returning the raw pointer lets pybind11
          }, //                                                                              correctly manage the memory
          "diffusion_factor"_a,
          "diffusion_tensor"_a,
          "space"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>());

      m.def(
          std::string(method_name).c_str(),
          [](const DF& diffusion_factor,
             const DT& diffusion_tensor,
             M& matrix,
             const R& space,
             const size_t over_integrate) {
            return make_elliptic_matrix_operator(diffusion_factor, diffusion_tensor, matrix, space, over_integrate)
                .release(); //                                                                     <- s.a. for release()
          },
          "diffusion_factor"_a,
          "diffusion_tensor"_a,
          "matrix"_a,
          "space"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>());
    } // ... addbind_factory_methods(...)
  }; // struct diffusion_switch

  struct diffusion_switch_scalar_base
  {
    template <class C>
    static void addbind_factory_methods(pybind11::module& m)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      const std::string method_name = "make_elliptic_matrix_operator_" + XT::LA::bindings::container_name<M>::value();

      m.def(method_name.c_str(),
            [](const DF& diffusion, const R& space, const size_t over_integrate) {
              return make_elliptic_matrix_operator<M>(diffusion, space, over_integrate).release(); // <- s.a.
            },
            "diffusion"_a,
            "space"_a,
            "over_integrate"_a = 0,
            py::keep_alive<0, 1>(),
            py::keep_alive<0, 2>());

      m.def(std::string(method_name).c_str(),
            [](const DF& diffusion, M& matrix, const R& space, const size_t over_integrate) {
              return make_elliptic_matrix_operator(diffusion, matrix, space, over_integrate).release(); // <- s.a.
            },
            "diffusion"_a,
            "matrix"_a,
            "space"_a,
            "over_integrate"_a = 0,
            py::keep_alive<0, 1>(),
            py::keep_alive<0, 2>());
    } // ... addbind_factory_methods(...)
  }; // struct diffusion_switch<..., void>

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
        "elliptic_matrix_operator_" + space_name<RP>::value() + "_" + XT::LA::bindings::container_name<M>::value() + "_"
        + diffusion_switch<>::suffix());

    auto c = MatrixOperatorBase<type>::bind(m, ClassName);

    diffusion_switch<>::template addbind_factory_methods<type>(m);

    return c;
  } // ... bind(...)
}; // class EllipticMatrixOperator


} // namespace bindings
} // namespace GDT
} // namespace Dune


// begin: this is what we need for the .so

#define _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_1D(_m, _d, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la)          \
  Dune::GDT::bindings::                                                                                                \
      EllipticMatrixOperator<Dune::XT::Functions::                                                                     \
                                 LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<                        \
                                                                  typename Dune::XT::Grid::                            \
                                                                      Layer<_GRID,                                     \
                                                                            Dune::XT::Grid::Layers::_layer,            \
                                                                            Dune::XT::Grid::Backends::_g_backend,      \
                                                                            Dune::XT::Grid::DD::                       \
                                                                                SubdomainGrid<_GRID>>::type>,          \
                                                              double,                                                  \
                                                              _d,                                                      \
                                                              double,                                                  \
                                                              1,                                                       \
                                                              1>,                                                      \
                             Dune::XT::Functions::                                                                     \
                                 LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<                        \
                                                                  typename Dune::XT::Grid::                            \
                                                                      Layer<_GRID,                                     \
                                                                            Dune::XT::Grid::Layers::_layer,            \
                                                                            Dune::XT::Grid::Backends::_g_backend,      \
                                                                            Dune::XT::Grid::DD::                       \
                                                                                SubdomainGrid<_GRID>>::type>,          \
                                                              double,                                                  \
                                                              _d,                                                      \
                                                              double,                                                  \
                                                              _d,                                                      \
                                                              _d>,                                                     \
                             Dune::GDT::SpaceProvider<_GRID,                                                           \
                                                      Dune::XT::Grid::Layers::_layer,                                  \
                                                      Dune::GDT::SpaceType::_s_type,                                   \
                                                      Dune::GDT::Backends::_s_backend,                                 \
                                                      _p,                                                              \
                                                      double,                                                          \
                                                      1,                                                               \
                                                      1>,                                                              \
                             typename Dune::XT::LA::Container<double,                                                  \
                                                              Dune::XT::LA::Backends::_la>::MatrixType>::bind(_m);     \
  Dune::GDT::bindings::                                                                                                \
      EllipticMatrixOperator<Dune::XT::Functions::                                                                     \
                                 LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<                        \
                                                                  typename Dune::XT::Grid::                            \
                                                                      Layer<_GRID,                                     \
                                                                            Dune::XT::Grid::Layers::_layer,            \
                                                                            Dune::XT::Grid::Backends::_g_backend,      \
                                                                            Dune::XT::Grid::DD::                       \
                                                                                SubdomainGrid<_GRID>>::type>,          \
                                                              double,                                                  \
                                                              _d,                                                      \
                                                              double,                                                  \
                                                              1,                                                       \
                                                              1>,                                                      \
                             void,                                                                                     \
                             Dune::GDT::SpaceProvider<_GRID,                                                           \
                                                      Dune::XT::Grid::Layers::_layer,                                  \
                                                      Dune::GDT::SpaceType::_s_type,                                   \
                                                      Dune::GDT::Backends::_s_backend,                                 \
                                                      _p,                                                              \
                                                      double,                                                          \
                                                      1,                                                               \
                                                      1>,                                                              \
                             typename Dune::XT::LA::Container<double,                                                  \
                                                              Dune::XT::LA::Backends::_la>::MatrixType>::bind(_m)

#define _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_D(_m, _d, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la)           \
  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_1D(_m, _d, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la);               \
  Dune::GDT::bindings::                                                                                                \
      EllipticMatrixOperator<Dune::XT::Functions::                                                                     \
                                 LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<                        \
                                                                  typename Dune::XT::Grid::                            \
                                                                      Layer<_GRID,                                     \
                                                                            Dune::XT::Grid::Layers::_layer,            \
                                                                            Dune::XT::Grid::Backends::_g_backend,      \
                                                                            Dune::XT::Grid::DD::                       \
                                                                                SubdomainGrid<_GRID>>::type>,          \
                                                              double,                                                  \
                                                              _d,                                                      \
                                                              double,                                                  \
                                                              _d,                                                      \
                                                              _d>,                                                     \
                             void,                                                                                     \
                             Dune::GDT::SpaceProvider<_GRID,                                                           \
                                                      Dune::XT::Grid::Layers::_layer,                                  \
                                                      Dune::GDT::SpaceType::_s_type,                                   \
                                                      Dune::GDT::Backends::_s_backend,                                 \
                                                      _p,                                                              \
                                                      double,                                                          \
                                                      1,                                                               \
                                                      1>,                                                              \
                             typename Dune::XT::LA::Container<double,                                                  \
                                                              Dune::XT::LA::Backends::_la>::MatrixType>::bind(_m)

//#if HAVE_ALBERTA
//#define _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_ALBERTA(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la)              \
//  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_D(_m, 2, ALBERTA_2D, _layer, _g_backend, _s_type, _s_backend, _p, _la)
//#else
#define _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_ALBERTA(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la)
//#endif

#if HAVE_DUNE_ALUGRID
#define _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_ALU(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la)                    \
  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_D(                                                                                 \
      _m, 2, ALU_2D_SIMPLEX_CONFORMING, _layer, _g_backend, _s_type, _s_backend, _p, _la)
#else
#define _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_ALU(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la)
#endif

//#if HAVE_DUNE_UGGRID || HAVE_UG
//#define _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_UG(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la)                   \
//  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_D(_m, 2, UG_2D, _layer, _g_backend, _s_type, _s_backend, _p, _la)
//#else
#define _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_UG(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la)
//#endif

#define _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_YASP(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la)
//  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_1D(                                                                                \
//      _m, 1, YASP_1D_EQUIDISTANT_OFFSET, _layer, _g_backend, _s_type, _s_backend, _p, _la);                            \
//  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_D(                                                                                 \
//      _m, 2, YASP_2D_EQUIDISTANT_OFFSET, _layer, _g_backend, _s_type, _s_backend, _p, _la)

#define _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_GRIDS(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la)                  \
  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_ALBERTA(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la);                     \
  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_ALU(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la);                         \
  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_UG(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la);                          \
  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_YASP(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la)

#if HAVE_DUNE_FEM
#define _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_FEM(_m, _la)                                                                 \
  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_GRIDS(_m, leaf, part, cg, fem, 1, _la);                                            \
  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_GRIDS(_m, level, part, cg, fem, 1, _la);                                           \
  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_GRIDS(_m, dd_subdomain, part, cg, fem, 1, _la);                                    \
  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_GRIDS(_m, leaf, part, dg, fem, 1, _la);                                            \
  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_GRIDS(_m, level, part, dg, fem, 1, _la);                                           \
  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_GRIDS(_m, dd_subdomain, part, dg, fem, 1, _la)
#define DUNE_GDT_OPERATORS_ELLIPTIC_BIND_FEM_COMMON(_m)
//_DUNE_GDT_OPERATORS_ELLIPTIC_BIND_FEM(_m, common_dense)
//#if HAVE_EIGEN
//#define DUNE_GDT_OPERATORS_ELLIPTIC_BIND_FEM_EIGEN(_m)                                                               \
//  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_FEM(_m, eigen_dense);                                                            \
//  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_FEM(_m, eigen_sparse)
//#else
#define DUNE_GDT_OPERATORS_ELLIPTIC_BIND_FEM_EIGEN(_m)
//#endif
#if HAVE_DUNE_ISTL
#define DUNE_GDT_OPERATORS_ELLIPTIC_BIND_FEM_ISTL(_m) _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_FEM(_m, istl_sparse)
#else
#define DUNE_GDT_OPERATORS_ELLIPTIC_BIND_FEM_ISTL(_m)
#endif
#else
#define DUNE_GDT_OPERATORS_ELLIPTIC_BIND_FEM_COMMON(_m)
#define DUNE_GDT_OPERATORS_ELLIPTIC_BIND_FEM_EIGEN(_m)
#define DUNE_GDT_OPERATORS_ELLIPTIC_BIND_FEM_ISTL(_m)
#endif

//#if HAVE_DUNE_PDELAB
//#define _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_PDELAB(_m, _la)                                                            \
//  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_ALBERTA(_m, leaf, view, cg, pdelab, 1, _la);                                     \
//  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_ALU(_m, leaf, view, cg, pdelab, 1, _la);                                         \
//  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_YASP(_m, leaf, view, cg, pdelab, 1, _la);                                        \
//  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_ALBERTA(_m, level, view, cg, pdelab, 1, _la);                                    \
//  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_ALU(_m, level, view, cg, pdelab, 1, _la);                                        \
//  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_YASP(_m, level, view, cg, pdelab, 1, _la)
//#define DUNE_GDT_OPERATORS_ELLIPTIC_BIND_PDELAB_COMMON(_m) _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_PDELAB(_m, common_dense)
//#if HAVE_EIGEN
//#define DUNE_GDT_OPERATORS_ELLIPTIC_BIND_PDELAB_EIGEN(_m)                                                            \
//  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_PDELAB(_m, eigen_dense);                                                         \
//  _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_PDELAB(_m, eigen_sparse)
//#else
//#define DUNE_GDT_OPERATORS_ELLIPTIC_BIND_PDELAB_EIGEN(_m)
//#endif
//#if HAVE_DUNE_ISTL
//#define DUNE_GDT_OPERATORS_ELLIPTIC_BIND_PDELAB_ISTL(_m) _DUNE_GDT_OPERATORS_ELLIPTIC_BIND_PDELAB(_m, istl_sparse)
//#else
//#define DUNE_GDT_OPERATORS_ELLIPTIC_BIND_PDELAB_ISTL(_m)
//#endif
//#else
#define DUNE_GDT_OPERATORS_ELLIPTIC_BIND_PDELAB_COMMON(_m)
#define DUNE_GDT_OPERATORS_ELLIPTIC_BIND_PDELAB_EIGEN(_m)
#define DUNE_GDT_OPERATORS_ELLIPTIC_BIND_PDELAB_ISTL(_m)
//#endif

// end: this is what we need for the .so


#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_OPERATORS_ELLIPTIC_BINDINGS_HH
