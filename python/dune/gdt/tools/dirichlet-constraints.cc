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

#include <dune/gdt/tools/dirichlet-constraints.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/parameter.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/la/traits.hh>


namespace Dune {
namespace GDT {
namespace bindings {


template <class GV, size_t r = 1>
class DirichletConstraints
{
  using G = std::decay_t<XT::Grid::extract_grid_t<GV>>;
  static const size_t d = G::dimension;
  using I = XT::Grid::extract_intersection_t<GV>;
  using S = SpaceInterface<GV, r>;

public:
  using type = GDT::DirichletConstraints<I, S>;
  using base_type = XT::Grid::ElementFunctor<GV>;
  using bound_type = pybind11::class_<type, base_type>;

private:
  template <class M>
  static void addbind_apply_matrix(bound_type& c)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    c.def(
        "apply",
        [](type& self, M& matrix, const bool only_clear, const bool ensure_symmetry) {
          self.apply(matrix, only_clear, ensure_symmetry);
        },
        "matrix"_a,
        "only_clear"_a = false,
        "ensure_symmetry"_a = true,
        py::call_guard<py::gil_scoped_release>());
  } // ... addbind_apply_matrix(...)

  template <class V>
  static void addbind_apply_vector(bound_type& c)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    c.def(
        "apply",
        [](type& self, V& vector) { self.apply(vector); },
        "vector"_a,
        py::call_guard<py::gil_scoped_release>());
  } // ... addbind_apply_vector(...)

  template <class M, class V>
  static void addbind_apply_matrix_and_vector(bound_type& c)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    c.def(
        "apply",
        [](type& self, M& matrix, V& vector, const bool only_clear, const bool ensure_symmetry) {
          self.apply(matrix, vector, only_clear, ensure_symmetry);
        },
        "matrix"_a,
        "vector"_a,
        "only_clear"_a = false,
        "ensure_symmetry"_a = true,
        py::call_guard<py::gil_scoped_release>());
  } // ... addbind_apply_matrix_and_vector(...)

public:
  static bound_type bind(pybind11::module& m,
                         const std::string& grid_id,
                         const std::string& layer_id = "",
                         const std::string& class_id = "dirichlet_constraints")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id;
    class_name += "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    class_name += "_" + XT::Common::to_string(r);
    class_name += "d_space";
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def(py::init<const XT::Grid::BoundaryInfo<I>&, const S&>(), "boundary_info"_a, "space"_a, py::keep_alive<1, 2>());
    c.def_property_readonly("boundary_info", &type::boundary_info);
    c.def_property_readonly("dirichlet_DoFs", &type::dirichlet_DoFs);
    //    addbind_apply_matrix<XT::LA::CommonDenseMatrix<double>>(c);
    //    addbind_apply_vector<XT::LA::CommonDenseVector<double>>(c);
    //    addbind_apply_matrix_and_vector<XT::LA::CommonDenseMatrix<double>, XT::LA::CommonDenseVector<double>>(c);
    //    //    addbind_apply_matrix<XT::LA::CommonSparseMatrix<double>>(c);
    //    //    addbind_apply_matrix_and_vector<XT::LA::CommonSparseMatrix<double>,
    //    XT::LA::CommonDenseVector<double>>(c);
    //#if HAVE_EIGEN
    //    addbind_apply_matrix<XT::LA::EigenDenseMatrix<double>>(c);
    //    addbind_apply_vector<XT::LA::EigenDenseVector<double>>(c);
    //    addbind_apply_matrix_and_vector<XT::LA::EigenDenseMatrix<double>, XT::LA::IstlDenseVector<double>>(c);
    //    addbind_apply_matrix<XT::LA::EigenRowMajorSparseMatrix<double>>(c);
    //    addbind_apply_matrix_and_vector<XT::LA::EigenRowMajorSparseMatrix<double>,
    //    XT::LA::EigenDenseVector<double>>(c);
    //#endif
    addbind_apply_matrix<XT::LA::IstlRowMajorSparseMatrix<double>>(c);
    addbind_apply_vector<XT::LA::IstlDenseVector<double>>(c);
    addbind_apply_matrix_and_vector<XT::LA::IstlRowMajorSparseMatrix<double>, XT::LA::IstlDenseVector<double>>(c);

    // factories
    const auto FactoryName = XT::Common::to_camel_case(class_id);
    m.def(
        FactoryName.c_str(),
        [](const XT::Grid::BoundaryInfo<I>& boundary_info, const S& space) { return new type(boundary_info, space); },
        "boundary_info"_a,
        "source_space"_a,
        py::keep_alive<0, 1>());

    return c;
  } // ... bind(...)
}; // class DirichletConstraints


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct DirichletConstraints_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    using Dune::GDT::bindings::DirichletConstraints;
    using Dune::XT::Grid::bindings::grid_name;

    DirichletConstraints<GV>::bind(m, grid_name<G>::value());
    if (d > 1)
      DirichletConstraints<GV, d>::bind(m, grid_name<G>::value());
    // add your extra dimensions here
    // ...
    DirichletConstraints_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct DirichletConstraints_for_all_grids<boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_tools_dirichlet_constraints, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");

  py::module::import("dune.gdt._spaces_interface");

  DirichletConstraints_for_all_grids<>::bind(m);
}
