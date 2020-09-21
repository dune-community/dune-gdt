// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#include "config.h"

#include <dune/xt/grid/grids.hh>
#include <python/dune/xt/la/traits.hh>

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/python.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/la/container.hh>
#include <python/dune/xt/common/exceptions.bindings.hh>
#include <python/dune/xt/functions/interfaces/grid-function.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/traits.hh>
#include <python/dune/xt/la/container.bindings.hh>

#include <dune/gdt/discretefunction/dof-vector.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class V, class GV>
class DofVector
{
  using type = GDT::DofVector<V, GV>;
  using G = typename GV::Grid;

public:
  using bound_type = pybind11::class_<type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& layer_id = "",
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value(),
                         const std::string& class_id = "dof_vector")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id + "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    class_name += "_" + XT::LA::bindings::container_name<V>::value();
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());

    // doing this so complicated to get an actual reference instead of a copy
    c.def_property("vector", (const V& (type::*)() const) & type::vector, (V & (type::*)()) & type::vector);

    c.def("resize_after_adapt", [](type& self) { self.resize_after_adapt(); });

    return c;
  } // ... bind(...)
}; // class DofVector


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class V, class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct DofVector_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GV = typename G::LeafGridView;

  static void bind(pybind11::module& m)
  {
    Dune::GDT::bindings::DofVector<V, GV>::bind(m);

    DofVector_for_all_grids<V, Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <class V>
struct DofVector_for_all_grids<V, Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_discretefunction_dof_vector, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");

  DofVector_for_all_grids<LA::CommonDenseVector<double>>::bind(m);
#if HAVE_EIGEN
  DofVector_for_all_grids<LA::EigenDenseVector<double>>::bind(m);
#endif
  DofVector_for_all_grids<LA::IstlDenseVector<double>>::bind(m);
}
