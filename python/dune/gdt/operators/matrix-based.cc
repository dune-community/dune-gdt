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
#include <dune/pybindxi/stl.h>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/la/type_traits.hh>

#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/tools/sparsity-pattern.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/parameter.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/walker.hh>
#include <python/dune/xt/la/traits.hh>
#include <python/dune/gdt/operators/interfaces.hh>


namespace Dune {
namespace GDT {
namespace bindings {


template <class M, class MatrixTag, class SparsityTag, class GV, size_t s_r = 1, size_t r_r = s_r>
class MatrixOperator
{
  using G = std::decay_t<XT::Grid::extract_grid_t<GV>>;
  static const size_t d = G::dimension;
  using GP = XT::Grid::GridProvider<G>;

public:
  using type = GDT::MatrixOperator<M, GV, s_r, 1, r_r>;
  using base_operator_type = GDT::OperatorInterface<M, GV, s_r, 1, r_r>;
  using base_walker_type = XT::Grid::Walker<GV>;
  using bound_type = pybind11::class_<type, base_operator_type, base_walker_type>;

private:
  using SS = typename type::SourceSpaceType;
  using RS = typename type::RangeSpaceType;
  using E = typename type::E;
  using F = typename type::F;

  template <bool needs_sparsity_tag = !std::is_same<SparsityTag, void>::value, bool anything = true>
  struct addbind /*<true, ...>*/
  {
    static void factory(pybind11::module& m, const std::string& FactoryName)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      m.def(FactoryName.c_str(),
            [](GP& grid,
               const SS& source_space,
               const RS& range_space,
               const MatrixTag&,
               const SparsityTag&,
               const XT::LA::SparsityPatternDefault& pattern) {
              return type(grid.leaf_view(), source_space, range_space, pattern);
            },
            "grid"_a,
            "source_space"_a,
            "range_space"_a,
            "unused_matrixbackend_tag"_a,
            "unused_sparsity_tag"_a,
            "sparsity_pattern"_a,
            py::keep_alive<0, 1>(),
            py::keep_alive<0, 2>(),
            py::keep_alive<0, 3>());
      m.def(FactoryName.c_str(),
            [](GP& grid, const SS& source_space, const RS& range_space, const MatrixTag&, const SparsityTag&) {
              return type(grid.leaf_view(),
                          source_space,
                          range_space,
                          make_element_and_intersection_sparsity_pattern(range_space, source_space, grid.leaf_view()));
            },
            "grid"_a,
            "source_space"_a,
            "range_space"_a,
            "unused_matrixbackend_tag"_a,
            "unused_sparsity_tag"_a,
            py::keep_alive<0, 1>(),
            py::keep_alive<0, 2>(),
            py::keep_alive<0, 3>());
    }
  }; // struct addbind<true, ...>

  template <bool anything>
  struct addbind<false, anything>
  {
    static void factory(pybind11::module& m, const std::string& FactoryName)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      m.def(FactoryName.c_str(),
            [](GP& grid,
               const SS& source_space,
               const RS& range_space,
               const MatrixTag&,
               const XT::LA::SparsityPatternDefault& pattern) {
              return type(grid.leaf_view(), source_space, range_space, pattern);
            },
            "grid"_a,
            "source_space"_a,
            "range_space"_a,
            "unused_matrixbackend_tag"_a,
            "sparsity_pattern"_a,
            py::keep_alive<0, 1>(),
            py::keep_alive<0, 2>(),
            py::keep_alive<0, 3>());
      m.def(FactoryName.c_str(),
            [](GP& grid, const SS& source_space, const RS& range_space, const MatrixTag&) {
              return type(grid.leaf_view(),
                          source_space,
                          range_space,
                          make_element_and_intersection_sparsity_pattern(range_space, source_space, grid.leaf_view()));
            },
            "grid"_a,
            "source_space"_a,
            "range_space"_a,
            "unused_matrixbackend_tag"_a,
            py::keep_alive<0, 1>(),
            py::keep_alive<0, 2>(),
            py::keep_alive<0, 3>());
    }
  }; // struct addbind<false, ...>

public:
  static bound_type bind(pybind11::module& m,
                         const std::string& matrix_id,
                         const std::string& class_id = "matrix_operator",
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value(),
                         const std::string& layer_id = "")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id;
    class_name += "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    class_name += "_" + XT::Common::to_string(r_r) + "d_range_space";
    class_name += "_" + XT::Common::to_string(s_r) + "d_source_space";
    class_name += "_" + matrix_id + "_matrix";
    const auto ClassName = XT::Common::to_camel_case(class_name);

    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def(py::init([](GP& grid, const SS& source_space, const RS& range_space, M& matrix) {
            return new type(grid.leaf_view(), source_space, range_space, matrix);
          }),
          "grid"_a,
          "source_space"_a,
          "range_space"_a,
          "matrix"_a,
          py::keep_alive<1, 2>(),
          py::keep_alive<1, 3>(),
          py::keep_alive<1, 4>(),
          py::keep_alive<1, 5>());
    c.def(
        py::init(
            [](GP& grid, const SS& source_space, const RS& range_space, const XT::LA::SparsityPatternDefault& pattern) {
              return new type(grid.leaf_view(), source_space, range_space, pattern);
            }),
        "grid"_a,
        "source_space"_a,
        "range_space"_a,
        "sparsity_pattern"_a,
        py::keep_alive<1, 2>(),
        py::keep_alive<1, 3>(),
        py::keep_alive<1, 4>());
    c.def(py::init([](GP& grid, const SS& source_space, const RS& range_space) {
            return new type(
                grid.leaf_view(),
                source_space,
                range_space,
                make_element_and_intersection_sparsity_pattern(range_space, source_space, grid.leaf_view()));
          }),
          "grid"_a,
          "source_space"_a,
          "range_space"_a,
          py::keep_alive<1, 2>(),
          py::keep_alive<1, 3>(),
          py::keep_alive<1, 4>());

    // doing this so complicated to get an actual reference instead of a copy
    c.def_property("matrix", (const M& (type::*)() const) & type::matrix, (M & (type::*)()) & type::matrix);

    // methods from walker base, to allow for overloads
    XT::Grid::bindings::Walker<G>::addbind_methods(c);

    // methods from operator base, to allow for overloads
    bindings::OperatorInterface<M, G, s_r, 1, r_r, 1>::addbind_methods(c);

    // additional methods
    c.def("clear", [](type& self) { self.clear(); });
    c.def("append",
          [](type& self,
             const LocalElementBilinearFormInterface<E, r_r, 1, F, F, s_r, 1, F>& local_bilinear_form,
             const XT::Common::Parameter& param,
             const XT::Grid::ElementFilter<GV>& filter) { self.append(local_bilinear_form, param, filter); },
          "local_element_bilinear_form"_a,
          "param"_a = XT::Common::Parameter(),
          "element_filter"_a = XT::Grid::ApplyOn::AllElements<GV>());
    c.def("assemble",
          [](type& self, const bool use_tbb) { self.assemble(use_tbb); },
          "parallel"_a = false,
          py::call_guard<py::gil_scoped_release>());

    // factories
    const auto FactoryName = XT::Common::to_camel_case(class_id);
    m.def(FactoryName.c_str(),
          [](GP& grid, const SS& source_space, const RS& range_space, M& matrix) {
            return type(grid.leaf_view(), source_space, range_space, matrix);
          },
          "grid"_a,
          "source_space"_a,
          "range_space"_a,
          "matrix"_a,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>(),
          py::keep_alive<0, 4>());
    addbind<>::factory(m, FactoryName);
    m.def(XT::Common::to_camel_case(matrix_id + "_" + class_id).c_str(),
          [](GP& grid, const SS& source_space, const RS& range_space, const XT::LA::SparsityPatternDefault& pattern) {
            return type(grid.leaf_view(), source_space, range_space, pattern);
          },
          "grid"_a,
          "source_space"_a,
          "range_space"_a,
          "sparsity_pattern"_a,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>());

    return c;
  } // ... bind(...)
}; // class MatrixOperator


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class M, class MT, class ST, class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct MatrixOperator_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m, const std::string& matrix_id)
  {
    Dune::GDT::bindings::MatrixOperator<M, MT, ST, GV>::bind(m, matrix_id);
    if (d > 1) {
      Dune::GDT::bindings::MatrixOperator<M, MT, ST, GV, d, 1>::bind(m, matrix_id);
      Dune::GDT::bindings::MatrixOperator<M, MT, ST, GV, 1, d>::bind(m, matrix_id);
      Dune::GDT::bindings::MatrixOperator<M, MT, ST, GV, d, d>::bind(m, matrix_id);
    }
    // add your extra dimensions here
    // ...
    MatrixOperator_for_all_grids<M, MT, ST, typename GridTypes::tail_type>::bind(m, matrix_id);
  }
};

template <class M, class MT, class ST>
struct MatrixOperator_for_all_grids<M, MT, ST, boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/, const std::string& /*matrix_id*/) {}
};


PYBIND11_MODULE(_operators_matrix_based, m)
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
  py::module::import("dune.gdt._operators_interfaces_common");
  py::module::import("dune.gdt._operators_interfaces_eigen");
  py::module::import("dune.gdt._operators_interfaces_istl");
  py::module::import("dune.gdt._spaces_interface");

  //  MatrixOperator_for_all_grids<LA::CommonDenseMatrix<double>,
  //                               LA::bindings::Common,
  //                               void,
  //                               XT::Grid::AvailableGridTypes>::bind(m, "common_dense");
  // Generic linear solver missing for CommonSparseMatrix!
  //  MatrixOperator_for_all_grids<LA::CommonSparseMatrix<double>, LA::bindings::Common, LA::bindings::Sparse,
  //  XT::Grid::AvailableGridTypes>::bind(m, "common_sparse");
  //#if HAVE_EIGEN
  //  MatrixOperator_for_all_grids<LA::EigenDenseMatrix<double>,
  //                               LA::bindings::Eigen,
  //                               LA::bindings::Dense,
  //                               XT::Grid::AvailableGridTypes>::bind(m, "eigen_dense");
  //  MatrixOperator_for_all_grids<LA::EigenRowMajorSparseMatrix<double>,
  //                               LA::bindings::Eigen,
  //                               LA::bindings::Sparse,
  //                               XT::Grid::AvailableGridTypes>::bind(m, "eigen_sparse");
  //#endif
  MatrixOperator_for_all_grids<LA::IstlRowMajorSparseMatrix<double>,
                               LA::bindings::Istl,
                               void,
                               XT::Grid::AvailableGridTypes>::bind(m, "istl_sparse");
}
