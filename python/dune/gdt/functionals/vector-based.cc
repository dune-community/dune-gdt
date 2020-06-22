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

#include <dune/gdt/functionals/vector-based.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/parameter.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/walker.hh>
#include <python/dune/xt/la/traits.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <class V, class VectorTag, class GV, size_t s_r = 1>
class VectorBasedFunctional
{
  using G = std::decay_t<XT::Grid::extract_grid_t<GV>>;
  static const size_t d = G::dimension;
  using GP = XT::Grid::GridProvider<G>;

public:
  using type = GDT::VectorBasedFunctional<V, GV, s_r>;
  using base_functional_type = GDT::FunctionalInterface<V, GV, s_r>;
  using base_walker_type = XT::Grid::Walker<GV>;
  using bound_type = pybind11::class_<type, base_functional_type, base_walker_type>;

private:
  using SS = typename type::SourceSpaceType;
  using E = typename type::E;
  using F = typename type::F;

public:
  static bound_type bind(pybind11::module& m,
                         const std::string& vector_id,
                         const std::string& class_id = "vector_functional",
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value(),
                         const std::string& layer_id = "")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id;
    class_name += "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    class_name += "_" + XT::Common::to_string(s_r) + "d_source_space";
    class_name += "_" + vector_id + "_vector";
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());

    c.def(py::init([](GP& grid, const SS& source_space, V& vector) {
            return new type(grid.leaf_view(), source_space, vector);
          }),
          "grid"_a,
          "source_space"_a,
          "vector"_a,
          py::keep_alive<1, 2>(),
          py::keep_alive<1, 3>(),
          py::keep_alive<1, 4>());
    c.def(py::init([](GP& grid, const SS& source_space) { return new type(grid.leaf_view(), source_space); }),
          "grid"_a,
          "source_space"_a,
          py::keep_alive<1, 2>(),
          py::keep_alive<1, 3>());

    // doing this so complicated to get an actual reference instead of a copy
    c.def_property("vector", (const V& (type::*)() const) & type::vector, (V & (type::*)()) & type::vector);

    // methods from walker base, to allow for overloads
    XT::Grid::bindings::Walker<G>::addbind_methods(c);

    // methods from functional base, to allow for overloads
    bindings::FunctionalInterface<V, G, s_r>::addbind_methods(c);

    // additional methods
    c.def("clear", [](type& self) { self.clear(); });
    c.def("append",
          [](type& self,
             const LocalElementFunctionalInterface<E, s_r, 1, F, typename type::DofFieldType>& local_functional,
             const XT::Common::Parameter& param,
             const XT::Grid::ElementFilter<GV>& filter) { self.append(local_functional, param, filter); },
          "local_element_functional"_a,
          "param"_a = XT::Common::Parameter(),
          "element_filter"_a = XT::Grid::ApplyOn::AllElements<GV>());
    c.def("__iadd__", // function ptr signature required for the right return type
          (type
           & (type::*)(const LocalElementFunctionalInterface<E, s_r, 1, F, typename type::DofFieldType>&,
                       const XT::Common::Parameter&,
                       const XT::Grid::ElementFilter<GV>&))
              & type::append,
          "local_element_functional"_a,
          "param"_a = XT::Common::Parameter(),
          "element_filter"_a = XT::Grid::ApplyOn::AllElements<GV>(),
          py::is_operator());
    c.def("assemble",
          [](type& self, const bool use_tbb) { self.assemble(use_tbb); },
          "parallel"_a = false,
          py::call_guard<py::gil_scoped_release>());

    // factories
    const auto FactoryName = XT::Common::to_camel_case(class_id);
    m.def(FactoryName.c_str(),
          [](GP& grid, const SS& source_space, V& vector) { return type(grid.leaf_view(), source_space, vector); },
          "grid"_a,
          "source_space"_a,
          "vector"_a,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>());
    m.def(FactoryName.c_str(),
          [](GP& grid, const SS& source_space, const VectorTag&) { return type(grid.leaf_view(), source_space); },
          "grid"_a,
          "source_space"_a,
          "vector_type"_a,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());
    m.def(XT::Common::to_camel_case(vector_id + "_" + class_id).c_str(),
          [](GP& grid, const SS& source_space) { return type(grid.leaf_view(), source_space); },
          "grid"_a,
          "source_space"_a,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());

    return c;
  } // ... bind(...)
}; // class VectorBasedFunctional


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class V, class VT, class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct VectorBasedFunctional_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m, const std::string& vector_id)
  {
    Dune::GDT::bindings::VectorBasedFunctional<V, VT, GV>::bind(m, vector_id);
    if (d > 1)
      Dune::GDT::bindings::VectorBasedFunctional<V, VT, GV, d>::bind(m, vector_id);
    // add your extra dimensions here
    // ...
    VectorBasedFunctional_for_all_grids<V, VT, typename GridTypes::tail_type>::bind(m, vector_id);
  }
};

template <class V, class VT>
struct VectorBasedFunctional_for_all_grids<V, VT, boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/, const std::string& /*vector_id*/) {}
};


PYBIND11_MODULE(_functionals_vector_based, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");

  //  py::module::import("dune.gdt._local_functionals_element_interface");
  py::module::import("dune.gdt._functionals_interfaces_common");
  py::module::import("dune.gdt._functionals_interfaces_eigen");
  py::module::import("dune.gdt._functionals_interfaces_istl");
  py::module::import("dune.gdt._spaces_interface");

  //  VectorBasedFunctional_for_all_grids<LA::CommonDenseVector<double>,
  //                               LA::bindings::Common,
  //                               XT::Grid::AvailableGridTypes>::bind(m, "common_dense");
  //#if HAVE_EIGEN
  //  VectorBasedFunctional_for_all_grids<LA::EigenDenseVector<double>,
  //                               LA::bindings::Eigen,
  //                               XT::Grid::AvailableGridTypes>::bind(m, "eigen_dense");
  //#endif
  VectorBasedFunctional_for_all_grids<LA::IstlDenseVector<double>, LA::bindings::Istl, XT::Grid::AvailableGridTypes>::
      bind(m, "istl");
}
