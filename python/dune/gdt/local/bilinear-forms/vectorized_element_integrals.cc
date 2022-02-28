// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#include "config.h"

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/functional.h>
#include <dune/pybindxi/numpy.h>
#include <dune/pybindxi/stl.h>

#include <dune/common/unused.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/grids.hh>

#include <dune/gdt/local/bilinear-forms/integrals.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/numpy.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/parameter.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/traits.hh>


namespace Dune {
namespace GDT {
namespace bindings {


template <class E,
          size_t t_r = 1,
          size_t t_rC = 1,
          class TR = double,
          class F = double,
          size_t a_r = t_r,
          size_t a_rC = t_rC,
          class AR = TR>
class VectorizedLocalElementIntegralBilinearForm
  : public LocalElementBilinearFormInterface<E, t_r, t_rC, TR, F, a_r, a_rC, AR>
{
  using ThisType = VectorizedLocalElementIntegralBilinearForm;
  using BaseType = LocalElementBilinearFormInterface<E, t_r, t_rC, TR, F, a_r, a_rC, AR>;

public:
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  // would have liked to use & for the bases (instead of *), but: https://github.com/pybind/pybind11/issues/1123
  using IntegrandType = std::function<pybind11::array_t<double>(const E&,
                                                                pybind11::array_t<double>,
                                                                const LocalTestBasisType*,
                                                                const LocalAnsatzBasisType*,
                                                                const XT::Common::Parameter&)>;

  using bound_type = pybind11::class_<ThisType, BaseType>;

  static bound_type
  bind(pybind11::module& m,
       const std::string& layer_id = "",
       const std::string& grid_id = XT::Grid::bindings::grid_name<XT::Grid::extract_grid_t<E>>::value(),
       const std::string& class_id = "vectorized_local_element_integral_bilinear_form")
  {
    static_assert(t_rC * a_rC == 1, "Not implemented for matrix-valued bases yet!");

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
    if (!std::is_same<TR, double>::value)
      test_string += "_" + XT::Common::Typename<TR>::value(/*fail_wo_typeid=*/true);
    test_string += "_test_basis";
    std::string ansatz_string = "";
    ansatz_string += "_" + XT::Common::to_string(a_r) + "d";
    if (a_rC > 1)
      ansatz_string += "x" + XT::Common::to_string(a_rC) + "d";
    if (!std::is_same<AR, double>::value)
      ansatz_string += "_" + XT::Common::Typename<AR>::value(/*fail_wo_typeid=*/true);
    ansatz_string += "_ansatz_basis";
    class_name += test_string;
    if (!test_string.empty() && !ansatz_string.empty())
      class_name += "_x";
    class_name += ansatz_string;
    class_name += "_to_scalar";
    if (!std::is_same<F, double>::value)
      class_name += "_" + XT::Common::Typename<F>::value(/*fail_wo_typeid=*/true);
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), class_id.c_str());
    c.def(py::init<IntegrandType, const int, const XT::Common::ParameterType&, const std::string&>(),
          "vectorized_binary_element_integrand"_a,
          "integrand_order"_a,
          "integrand_parameters"_a = XT::Common::ParameterType(),
          "logging_prefix"_a = "");

    // factory
    using G = XT::Grid::extract_grid_t<E>;
    using GP = XT::Grid::GridProvider<G>;
    if constexpr (t_r == 1 && a_r == 1)
      m.def(
          XT::Common::to_camel_case(class_id).c_str(),
          [](const GP& /*grid*/,
             IntegrandType integrand,
             const int integrand_order,
             XT::Grid::bindings::Dimension<1> /*test_dim_range*/,
             XT::Grid::bindings::Dimension<1> /*ansatz_dim_range*/,
             const XT::Common::ParameterType& integrand_parameters,
             const std::string& logging_prefix = "") {
            return new ThisType(integrand, integrand_order, integrand_parameters, logging_prefix);
          },
          "grid"_a,
          "vectorized_binary_element_integrand"_a,
          "integrand_order"_a,
          "test_dim_range"_a = XT::Grid::bindings::Dimension<1>(),
          "ansatz_dim_range"_a = XT::Grid::bindings::Dimension<1>(),
          "integrand_parameters"_a = XT::Common::ParameterType(),
          "logging_prefix"_a = "");
    else
      m.def(
          XT::Common::to_camel_case(class_id).c_str(),
          [](const GP& /*grid*/,
             IntegrandType integrand,
             const int integrand_order,
             XT::Grid::bindings::Dimension<t_r> /*test_dim_range*/,
             XT::Grid::bindings::Dimension<a_r> /*ansatz_dim_range*/,
             const XT::Common::ParameterType& integrand_parameters,
             const std::string& logging_prefix = "") {
            return new ThisType(integrand, integrand_order, integrand_parameters, logging_prefix);
          },
          "grid"_a,
          "vectorized_element_integrand"_a,
          "integrand_order"_a,
          "test_dim_range"_a,
          "ansatz_dim_range"_a,
          "integrand_parameters"_a = XT::Common::ParameterType(),
          "logging_prefix"_a = "");
    return c;
  } // ... bind(...)

  VectorizedLocalElementIntegralBilinearForm(
      IntegrandType integrand,
      const int integrand_order = 0,
      const XT::Common::ParameterType& integrand_parameters = {},
      const std::string& logging_prefix = "",
      const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
    : BaseType(integrand_parameters,
               logging_prefix.empty() ? "VectorizedLocalElementIntegralBilinearForm" : logging_prefix,
               logging_state)
    , integrand_(integrand)
    , integrand_order_(integrand_order)
  {
    LOG_(info) << "VectorizedLocalElementIntegralBilinearForm(this=" << this << ", integrand=" << &integrand
               << ", integrand_order=" << integrand_order << ", integrand_parameters=" << integrand_parameters << ")"
               << std::endl;
  }

  VectorizedLocalElementIntegralBilinearForm(const ThisType& other) = default;

  VectorizedLocalElementIntegralBilinearForm(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy() const override final
  {
    LOG_(debug) << "copy()" << std::endl;
    return std::make_unique<ThisType>(*this);
  }

  using BaseType::apply2;

  void apply2(const LocalTestBasisType& test_basis,
              const LocalAnsatzBasisType& ansatz_basis,
              DynamicMatrix<F>& result,
              const XT::Common::Parameter& param = {}) const override final
  {
    namespace py = pybind11;
    LOG_(debug) << "apply2(test_basis.size()=" << test_basis.size(param)
                << ", ansatz_basis.size()=" << ansatz_basis.size(param) << ", param=" << param << ")" << std::endl;
    // prepare integand
    const auto& element = ansatz_basis.element();
    assert(test_basis.element() == element && "This must not happen!");
    //    integrand_->bind(element);
    // prepare storage
    const size_t rows = test_basis.size(param);
    const size_t cols = ansatz_basis.size(param);
    if (result.rows() < rows || result.cols() < cols)
      result.resize(rows, cols);
    result *= 0;
    // collect quadrature data
    // - therefore, count how many quadrature points there are
    size_t num_quadrature_points = 0;
    for (const auto& DUNE_UNUSED(quadrature_point) : QuadratureRules<D, d>::rule(element.type(), integrand_order_))
      ++num_quadrature_points;
    std::vector<double> integration_factor(num_quadrature_points);
    std::vector<double> quadrature_weight(num_quadrature_points);
    // - and store them in a numpy.ndarray
    py::array_t<double> quadrature_points(/*shape=*/{num_quadrature_points, size_t(d)});
    auto access_to_quadrature_points = quadrature_points.mutable_unchecked<2>();
    size_t pp = 0;
    for (const auto& quadrature_point : QuadratureRules<D, d>::rule(element.type(), integrand_order_)) {
      integration_factor[pp] = element.geometry().integrationElement(quadrature_point.position());
      quadrature_weight[pp] = quadrature_point.weight();
      for (size_t ii = 0; ii < d; ++ii)
        access_to_quadrature_points(pp, ii) = quadrature_point.position()[ii];
      pp += 1;
    }
    // evaluate integrand
    const auto values = integrand_(element, quadrature_points, &test_basis, &ansatz_basis, param);
    const auto& access_to_values = XT::Common::bindings::access_array</*ndim=*/3>(
        /*array=*/values,
        /*required_shape=*/{num_quadrature_points, rows, cols},
        /*array_name=*/"Result of vectorized_binary_element_integrand",
        /*required_shape_docs=*/"(num_quadrature_points, test_basis.size(mu), ansatz_basis.size(mu))");
    // compute integral
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        for (size_t qq = 0; qq < num_quadrature_points; ++qq)
          result[ii][jj] += access_to_values(qq, ii, jj) * integration_factor[qq] * quadrature_weight[qq];

    LOG_(debug) << "  result = " << result << std::endl;
  } // ... apply2(...)

private:
  const IntegrandType integrand_;
  const int integrand_order_;
}; // class VectorizedLocalElementIntegralBilinearForm


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct VectorizedLocalElementIntegralBilinearForm_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GV = typename G::LeafGridView;
  using E = Dune::XT::Grid::extract_entity_t<GV>;
  static const constexpr size_t d = G::dimension;
  using F = double;

  static void bind(pybind11::module& m)
  {
    using Dune::GDT::bindings::VectorizedLocalElementIntegralBilinearForm;

    VectorizedLocalElementIntegralBilinearForm<E>::bind(m);
    if (d > 1) {
      VectorizedLocalElementIntegralBilinearForm<E, 1, 1, F, F, d, 1, F>::bind(m);
      //      VectorizedLocalElementIntegralBilinearForm<E, 1, 1, F, F, d, d, F>::bind(m);
      VectorizedLocalElementIntegralBilinearForm<E, d, 1, F, F, 1, 1, F>::bind(m);
      VectorizedLocalElementIntegralBilinearForm<E, d, 1, F, F, d, 1, F>::bind(m);
      //      VectorizedLocalElementIntegralBilinearForm<E, d, 1, F, F, d, d, F>::bind(m);
      //      VectorizedLocalElementIntegralBilinearForm<E, d, d, F, F, 1, 1, F>::bind(m);
      //      VectorizedLocalElementIntegralBilinearForm<E, d, d, F, F, d, 1, F>::bind(m);
      //      VectorizedLocalElementIntegralBilinearForm<E, d, d, F, F, d, d, F>::bind(m);
    }
    // add your extra dimensions here
    // ...
    VectorizedLocalElementIntegralBilinearForm_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct VectorizedLocalElementIntegralBilinearForm_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_local_bilinear_forms_vectorized_element_integrals, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");
  py::module::import("dune.gdt._local_bilinear_forms_element_interface");

  VectorizedLocalElementIntegralBilinearForm_for_all_grids<XT::Grid::bindings::AvailableGridTypes>::bind(m);
}
