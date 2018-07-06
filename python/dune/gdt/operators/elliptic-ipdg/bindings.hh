// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   Rene Milk       (2018)

#ifndef PYTHON_DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BINDINGS_HH
#define PYTHON_DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BINDINGS_HH

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/string.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/type_traits.hh>
#include <python/dune/xt/la/container.bindings.hh>

#include <python/dune/gdt/spaces/interface.hh>
#include <dune/gdt/type_traits.hh>

#include <dune/gdt/operators/elliptic-ipdg.hh>
#include <python/dune/gdt/operators/base.hh>

namespace Dune {
namespace GDT {
namespace bindings {
namespace internal {


template <class DF,
          typename DT, // may be void
          class R,
          LocalEllipticIpdgIntegrands::Method method,
          class M,
          class GL /*,
          class S = R,
          class F = typename RP::type::RangeFieldType*/>
class EllipticIpdgMatrixOperator
{
  static_assert(is_space<R>::value, "");
  static_assert(XT::Grid::is_layer<GL>::value, "");

public:
  typedef GDT::EllipticIpdgMatrixOperator<DF, DT, R, method, M, GL /*, S, F*/> type;
  using bound_type = typename MatrixOperatorBase<type>::bound_type;

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
          "make_elliptic_" + LocalEllipticIpdgIntegrands::method_name<method>::value() + "_matrix_operator";

      m.def(std::string(method_name + "_" + XT::LA::bindings::container_name<M>::value()).c_str(),
            [](const DF& diffusion_factor,
               const DT& diffusion_tensor,
               const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GL>>& boundary_info,
               const R& space,
               const size_t over_integrate) {
              return new type(over_integrate, boundary_info, diffusion_factor, diffusion_tensor, space);
            },
            "diffusion_factor"_a,
            "diffusion_tensor"_a,
            "boundary_info"_a,
            "space"_a,
            "over_integrate"_a = 0,
            py::keep_alive<0, 1>(),
            py::keep_alive<0, 2>(),
            py::keep_alive<0, 3>(),
            py::keep_alive<0, 4>(),
            py::return_value_policy::take_ownership);

      m.def(method_name.c_str(),
            [](const DF& diffusion_factor,
               const DT& diffusion_tensor,
               const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GL>>& boundary_info,
               M& matrix,
               const R& space,
               const size_t over_integrate) {
              return new type(over_integrate, boundary_info, diffusion_factor, diffusion_tensor, matrix, space);
            },
            "diffusion_factor"_a,
            "diffusion_tensor"_a,
            "boundary_info"_a,
            "matrix"_a,
            "space"_a,
            "over_integrate"_a = 0,
            py::keep_alive<0, 1>(),
            py::keep_alive<0, 2>(),
            py::keep_alive<0, 3>(),
            py::keep_alive<0, 4>(),
            py::keep_alive<0, 5>(),
            py::return_value_policy::take_ownership);
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
          "make_elliptic_" + LocalEllipticIpdgIntegrands::method_name<method>::value() + "_matrix_operator";

      m.def(
          std::string(method_name + "_" + XT::LA::bindings::container_name<M>::value()).c_str(),
          [](const DF& diffusion,
             const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GL>>& boundary_info,
             const R& space,
             const size_t over_integrate) {
            return make_elliptic_ipdg_matrix_operator<M, method>(diffusion, boundary_info, space, over_integrate)
                .release(); //                                                                     <- s.a. for release()
          },
          "diffusion"_a,
          "boundary_info"_a,
          "space"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>());

      m.def(
          method_name.c_str(),
          [](const DF& diffusion,
             const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GL>>& boundary_info,
             M& matrix,
             const R& space,
             const size_t over_integrate) {
            return make_elliptic_ipdg_matrix_operator<method>(diffusion, boundary_info, matrix, space, over_integrate)
                .release(); //                                                                     <- s.a. for release()
          },
          "diffusion"_a,
          "boundary_info"_a,
          "matrix"_a,
          "space"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>(),
          py::keep_alive<0, 4>());
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
  static bound_type bind(pybind11::module& m, std::string class_name)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;
    MatrixOperatorBase<type>::bind_bases(m);

    bound_type c(m, class_name.c_str(), class_name.c_str());
    MatrixOperatorBase<type>::bind(c);
    diffusion_switch<>::template addbind_factory_methods<type>(m);

    c.def("assemble", [](type* self) { self->assemble(false); });
    return c;
  } // ... bind(...)

  //! this is only needed for RS2017 op which cannot call some ctor variant in the fac bindings
  static bound_type bind_no_factories(pybind11::module& m, std::string class_name)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;
    MatrixOperatorBase<type>::bind_bases(m);

    bound_type c(m, class_name.c_str(), class_name.c_str());
    MatrixOperatorBase<type>::bind(c);

    c.def("assemble", [](type* self) { self->assemble(false); });
    return c;
  } // ... bind(...)
}; // EllipticIpdgMatrixOperator


} // namespace internal


template <class G,
          XT::Grid::Layers gl,
          XT::Grid::Backends gl_backend,
          bool with_tensor,
          LocalEllipticIpdgIntegrands::Method method,
          XT::LA::Backends la,
          Backends space_backend,
          SpaceType space_type,
          XT::Grid::Layers space_layer,
          int space_polorder>
class EllipticIpdgMatrixOperator
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

  using binder = internal::EllipticIpdgMatrixOperator<DF, DT, S, method, M, GL>;

public:
  using bound_type = typename binder::bound_type;

  static bound_type bind(pybind11::module& m)
  {
    using df = typename binder::
        template diffusion_switch<std::is_same<DT, void>::value, (DF::dimRange == 1 && DF::dimRangeCols == 1), false>;
    const std::string class_name = XT::Common::to_camel_case(
        "elliptic_" + LocalEllipticIpdgIntegrands::method_name<method>::value() + "_matrix_operator_"
        + space_name<SP>::value()
        + "_"
        + XT::LA::bindings::container_name<M>::value()
        + "_"
        + df::suffix());

    return bind(m, class_name);
  }

  static bound_type bind(pybind11::module& m, std::string name)
  {
    return binder::bind(m, name);
  }
}; // class EllipticIpdgMatrixOperator


} // naemspace bindings
} // namespace GDT
} // namespace Dune


#endif // PYTHON_DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BINDINGS_HH
