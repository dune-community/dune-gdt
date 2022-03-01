// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2021)

#include "config.h"

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/grid/grids.hh>
#include <python/dune/xt/la/traits.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

#include <dune/xt/common/python.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/la/container.hh>
#include <python/dune/xt/functions/interfaces/grid-function.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/traits.hh>
#include <python/dune/xt/la/traits.hh>
#include <python/dune/xt/la/container.bindings.hh>

#include <dune/gdt/discretefunction/bochner.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class V, class MatchingVectorTag, class GV, size_t r = 1, size_t rC = 1, class R = double>
class DiscreteBochnerFunction
{
  using type = GDT::DiscreteBochnerFunction<V, GV, r, rC, R>;
  using G = typename GV::Grid;
  using BS = BochnerSpace<GV, r, rC, R>;

public:
  using bound_type = pybind11::class_<type>;

  static V& get_vec_ref(pybind11::handle list_element)
  {
    try {
      return list_element.cast<V&>();
    } catch (...) {
    }
    // if we came that far the above did not work, try
    try {
      return list_element.attr("impl").cast<V&>();
    } catch (...) {
    }
    // if we came that far the above did not work, give up
    DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
               "Could convert Python object of type [1] to C++ type [2] when trying to convert list_of_vectors:"
                   << "\n\n"
                   << "  - [1]: " << pybind11::str(list_element.get_type()).cast<std::string>() << "\n"
                   << "  - [2]: " << XT::Common::Typename<V>::value() << "&");
  } // ... get_vec_ref(...)

  static bound_type bind(pybind11::module& m,
                         const std::string& layer_id = "",
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value(),
                         const std::string& class_id = "discrete_bochner_function")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    std::string class_name = class_id + "_" + grid_id;
    if (!layer_id.empty())
      class_name += "_" + layer_id;
    if (r > 1)
      class_name += "_to_" + XT::Common::to_string(r) + "d";
    if (rC > 1)
      class_name += "x" + XT::Common::to_string(rC) + "d";
    class_name += "_" + XT::LA::bindings::container_name<V>::value();
    const auto ClassName = XT::Common::to_camel_case(class_name);
    const std::string default_name = "dune.gdt.discretebochnerfunction";
    bound_type c(m, ClassName.c_str(), ClassName.c_str());

    c.def_property_readonly("space", &type::space);
    c.def_property_readonly("name", &type::name);

    c.def("evaluate", &type::evaluate, "time"_a);
    c.def(
        "visualize",
        [](type& self, const std::string& filename_prefix) {
          return visualize(self, filename_prefix, true, VTK::appendedraw);
        },
        "filename_prefix"_a);

    m.def(
        XT::Common::to_camel_case(class_id).c_str(),
        [](const BS& bochner_space, py::list list_of_vectors, const MatchingVectorTag&, const std::string& name) {
          // try to interprete list_of_vectors as a ListVectorArray of correct length
          DUNE_THROW_IF(list_of_vectors.size() != bochner_space.temporal_space().mapper().size(),
                        XT::Common::Exceptions::wrong_input_given,
                        "list_of_vectors does not match bochner_space:"
                            << "\n"
                            << "  bochner_space.spatial_space().mapper().size() = "
                            << bochner_space.spatial_space().mapper().size() << "\n"
                            << "  list_of_vectors.size() = " << list_of_vectors.size());
          auto vector_array =
              std::make_unique<XT::LA::ListVectorArray<V>>(/*dim=*/bochner_space.spatial_space().mapper().size(),
                                                           /*length=*/0,
                                                           /*reserve=*/bochner_space.temporal_space().mapper().size());
          size_t counter = 0;
          auto time_points = bochner_space.time_points();
          for (py::handle list_element : list_of_vectors) {
            auto& vec = get_vec_ref(list_element);
            vector_array->append(vec, {{"_t", {time_points[counter]}}});
          }
          return new type(bochner_space, vector_array.release(), name);
        },
        "space"_a,
        "list_of_vectors"_a,
        "vector_type"_a,
        "name"_a = default_name,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>());
    if (std::is_same<MatchingVectorTag, XT::LA::bindings::Istl>::value)
      m.def(
          XT::Common::to_camel_case(class_id).c_str(),
          [](const BS& bochner_space, py::list list_of_vectors, const std::string& name, const MatchingVectorTag&) {
            // try to interprete list_of_vectors as a ListVectorArray of correct length
            DUNE_THROW_IF(list_of_vectors.size() != bochner_space.temporal_space().mapper().size(),
                          XT::Common::Exceptions::wrong_input_given,
                          "list_of_vectors does not match bochner_space:"
                              << "\n"
                              << "  bochner_space.spatial_space().mapper().size() = "
                              << bochner_space.spatial_space().mapper().size() << "\n"
                              << "  list_of_vectors.size() = " << list_of_vectors.size());
            auto vector_array = std::make_unique<XT::LA::ListVectorArray<V>>(
                /*dim=*/bochner_space.spatial_space().mapper().size(),
                /*length=*/0,
                /*reserve=*/bochner_space.temporal_space().mapper().size());
            size_t counter = 0;
            auto time_points = bochner_space.time_points();
            for (py::handle list_element : list_of_vectors) {
              auto& vec = get_vec_ref(list_element);
              vector_array->append(vec, {{"_t", {time_points[counter]}}});
            }
            return new type(bochner_space, vector_array.release(), name);
          },
          "space"_a,
          "list_of_vectors"_a,
          "name"_a = default_name,
          "vector_type"_a = XT::LA::bindings::Istl(),
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());

    m.def(
        XT::Common::to_camel_case(class_id).c_str(),
        [](const BS& space, const MatchingVectorTag&, const std::string& name) {
          return make_discrete_bochner_function<V>(space, name);
        },
        "space"_a,
        "vector_type"_a,
        "name"_a = default_name,
        py::keep_alive<0, 1>());
    if (std::is_same<MatchingVectorTag, XT::LA::bindings::Istl>::value)
      m.def(
          XT::Common::to_camel_case(class_id).c_str(),
          [](const BS& space, const std::string& name, const MatchingVectorTag&) {
            return make_discrete_bochner_function<V>(space, name);
          },
          "space"_a,
          "name"_a = default_name,
          "vector_type"_a = XT::LA::bindings::Istl(),
          py::keep_alive<0, 1>());
    return c;
  } // ... bind(...)
}; // class DiscreteBochnerFunction


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class V, class VT, class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct DiscreteBochnerFunction_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GV = typename G::LeafGridView;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    using Dune::GDT::bindings::DiscreteBochnerFunction;

    DiscreteBochnerFunction<V, VT, GV>::bind(m);
    if (d > 1)
      DiscreteBochnerFunction<V, VT, GV, d>::bind(m);
    // add your extra dimensions here
    // ...
    DiscreteBochnerFunction_for_all_grids<V, VT, Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <class V, class VT>
struct DiscreteBochnerFunction_for_all_grids<V, VT, Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_discretefunction_bochner, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");

  py::module::import("dune.gdt._discretefunction_discretefunction");
  py::module::import("dune.gdt._spaces_bochner");

  DiscreteBochnerFunction_for_all_grids<LA::CommonDenseVector<double>,
                                        LA::bindings::Common,
                                        XT::Grid::bindings::AvailableGridTypes>::bind(m);
#if HAVE_EIGEN
  DiscreteBochnerFunction_for_all_grids<LA::EigenDenseVector<double>,
                                        LA::bindings::Eigen,
                                        XT::Grid::bindings::AvailableGridTypes>::bind(m);
#endif
  DiscreteBochnerFunction_for_all_grids<LA::IstlDenseVector<double>,
                                        LA::bindings::Istl,
                                        XT::Grid::bindings::AvailableGridTypes>::bind(m);
}
