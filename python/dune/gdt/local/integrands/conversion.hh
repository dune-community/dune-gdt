// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_GDT_LOCAL_INTEGRANDS_CONVERSION_HH
#define PYTHON_DUNE_GDT_LOCAL_INTEGRANDS_CONVERSION_HH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/gridprovider/provider.hh>

#include <dune/gdt/local/integrands/conversion.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>


namespace Dune {
namespace GDT {
namespace bindings {


template <class E,
          size_t t_r = 1,
          size_t t_rC = 1,
          class TF = double,
          class F = double,
          size_t a_r = t_r,
          size_t a_rC = t_rC,
          class AF = TF>
class LocalBinaryToUnaryElementIntegrand
{
  using G = std::decay_t<XT::Grid::extract_grid_t<E>>;
  using GP = XT::Grid::GridProvider<G>;
  static const size_t d = G::dimension;

public:
  using type = GDT::LocalBinaryToUnaryElementIntegrand<E, t_r, t_rC, TF, F, a_r, a_rC, AF>;
  using base_type = GDT::LocalUnaryElementIntegrandInterface<E, t_r, t_rC, TF, F>;
  using bound_type = pybind11::class_<type, base_type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "local_binary_to_unary_element_integrand",
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value(),
                         const std::string& layer_id = "")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id;
    class_name += "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    std::string test_string = "";
    test_string += "_" + XT::Common::to_string(t_r) + "d";
    if (t_rC > 1)
      test_string += "x" + XT::Common::to_string(t_rC) + "d";
    if (!std::is_same<TF, double>::value)
      test_string += "_" + XT::Common::Typename<TF>::value(/*fail_wo_typeid=*/true);
    test_string += "_test_basis";
    std::string ansatz_string = "";
    ansatz_string += "_" + XT::Common::to_string(a_r) + "d";
    if (a_rC > 1)
      ansatz_string += "x" + XT::Common::to_string(a_rC) + "d";
    if (!std::is_same<AF, double>::value)
      ansatz_string += "_" + XT::Common::Typename<AF>::value(/*fail_wo_typeid=*/true);
    ansatz_string += "_ansatz_basis";
    class_name += test_string;
    if (!test_string.empty() && !ansatz_string.empty())
      class_name += "_x";
    class_name += ansatz_string;
    class_name += "_to_scalar";
    if (!std::is_same<F, double>::value)
      class_name += "_" + XT::Common::Typename<F>::value(/*fail_wo_typeid=*/true);
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def(py::init<const typename type::LocalBinaryElementIntegrandType&,
                   const XT::Functions::GridFunctionInterface<E, a_r, a_rC, AF>&,
                   const std::string&>(),
          "local_binary_integrand"_a,
          "inducing_function_as_ansatz_basis"_a,
          "logging_prefix"_a = "",
          py::keep_alive<1, 2>(),
          py::keep_alive<1, 3>());

    m.def(XT::Common::to_camel_case(class_id).c_str(),
          [](const typename type::LocalBinaryElementIntegrandType& local_binary_integrand,
             const XT::Functions::GridFunctionInterface<E, a_r, a_rC, AF>& inducing_function_as_ansatz_basis,
             const std::string& logging_prefix) {
            return type(local_binary_integrand, inducing_function_as_ansatz_basis, logging_prefix);
          },
          "local_binary_integrand"_a,
          "inducing_function_as_ansatz_basis"_a,
          "logging_prefix"_a = "",
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());

    return c;
  } // ... bind(...)
}; // class LocalBinaryToUnaryElementIntegrand


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // PYTHON_DUNE_GDT_LOCAL_INTEGRANDS_CONVERSION_HH
