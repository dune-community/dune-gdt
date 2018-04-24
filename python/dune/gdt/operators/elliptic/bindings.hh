// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   Rene Milk       (2018)

#ifndef PYTHON_DUNE_GDT_OPERATORS_ELLIPTIC_BINDINGS_HH
#define PYTHON_DUNE_GDT_OPERATORS_ELLIPTIC_BINDINGS_HH

#include <dune/pybindxi/pybind11.h>

#include <python/dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/type_traits.hh>
#include <python/dune/xt/la/container.bindings.hh>

#include <python/dune/gdt/spaces/interface.hh>
#include <dune/gdt/type_traits.hh>

#include <python/dune/gdt/operators/base.hh>
#include <dune/gdt/operators/elliptic.hh>

namespace Dune {
namespace GDT {
namespace bindings {
namespace internal {


template <class DF,
          typename DT, // may be void
          class R,
          class M,
          class GL /*,
          class S,
          class F*/>
class EllipticMatrixOperator
{
  static_assert(is_space<R>::value, "");
  static_assert(XT::Grid::is_layer<GL>::value, "");

public:
  typedef GDT::EllipticMatrixOperator<DF, DT, R, M, GL /*, S, F*/> type;
  using bound_type = typename bindings::MatrixOperatorBase<type>::bound_type;

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
  static bound_type bind(pybind11::module& m, const std::string& space_name, const std::string& grid_layer_name)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = XT::Common::to_camel_case(
        "elliptic_matrix_operator_" + space_name + "_" + XT::LA::bindings::container_name<M>::value() + "_"
        + diffusion_switch<>::suffix());

    bound_type c = bindings::MatrixOperatorBase<type>::bind(m, ClassName, space_name, space_name, grid_layer_name);

    diffusion_switch<>::template addbind_factory_methods<type>(m);

    return c;
  } // ... bind(...)
}; // class EllipticMatrixOperator


} // namespace internal


template <class G,
          XT::Grid::Layers gl,
          XT::Grid::Backends gl_backend,
          bool with_tensor,
          XT::LA::Backends la,
          Backends space_backend,
          SpaceType space_type,
          XT::Grid::Layers space_layer,
          int space_polorder>
class EllipticMatrixOperator
{
  using GL = typename XT::Grid::Layer<G, gl, gl_backend, XT::Grid::DD::SubdomainGrid<G>>::type;
  using E = XT::Grid::extract_entity_t<GL>;
  using D = typename G::ctype;
  static const constexpr size_t d = G::dimension;
  using R = double;
  using DF = XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1>;
  using DT = std::conditional_t<with_tensor, XT::Functions::LocalizableFunctionInterface<E, D, d, R, d, d>, void>;
  using SP = SpaceProvider<G, space_layer, space_type, space_backend, space_polorder, R, 1, 1>;
  using S = typename SP::type;
  using M = typename Dune::XT::LA::Container<R, la>::MatrixType;

  using binder = internal::EllipticMatrixOperator<DF, DT, S, M, GL>;

public:
  using bound_type = typename binder::bound_type;

  static bound_type bind(pybind11::module& m)
  {
    return binder::bind(m,
                        space_name<SP>::value(),
                        XT::Grid::bindings::layer_name<gl>::value()
                            + XT::Grid::bindings::backend_name<gl_backend>::value());
  }
}; // class EllipticMatrixOperator

} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // PYTHON_DUNE_GDT_OPERATORS_ELLIPTIC_BINDINGS_HH
