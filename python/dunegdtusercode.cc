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
#include <dune/pybindxi/stl.h>
#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/parameter.hh>
#include <dune/xt/la/container/istl.hh>
#include <dune/xt/common/print.hh>
#include <dune/xt/grid/type_traits.hh>

// new for partition of unity
#include <dune/gdt/interpolations.hh>
#include <dune/gdt/spaces/l2/discontinuous-lagrange.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
#include <dune/xt/grid/dd/glued.hh>


namespace py = pybind11;
using namespace Dune;


template <class V>
V& get_vec_ref(py::handle list_element)
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
             "Cannot convert Python object of type [1] to C++ type [2]:"
                 << "\n\n"
                 << "  - [1]: " << list_element.get_type() << "\n"
                 << "  - [2]: " << XT::Common::Typename<V>::value() << "&");
} // ... get_vec_ref(...)

using G = YASP_2D_EQUIDISTANT_OFFSET;
using MG = YASP_2D_EQUIDISTANT_OFFSET;
using GridGlueType = XT::Grid::DD::Glued<MG, G, XT::Grid::Layers::leaf>;

PYBIND11_MODULE(dunegdtusercode, m)
{
  using namespace pybind11::literals;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");
  py::module::import("dune.gdt");

  m.def(
      "list_test",
      [](py::list list_of_vectors) {
        using V = XT::LA::IstlDenseVector<double>;
        std::cout << "list has " << list_of_vectors.size() << " elements:" << std::endl;
        for (py::handle obj : list_of_vectors) { // iterators!
          std::cout << "  - " << obj.attr("__str__")().cast<std::string>() << std::endl;
          auto& vec = get_vec_ref<V>(obj);
          std::cout << "  = " << XT::Common::print(vec) << std::endl;
          vec[0] = 100;
          std::cout << "  = " << XT::Common::print(vec) << std::endl;
        }
      },
      "param"_a);
  
  m.def("compute_partition_of_unity",
        [](GridGlueType& dd_grid, const size_t I /*, const std::string space_type */) {
          const auto& coarse_grid_view = dd_grid.macro_grid_view();
          auto coarse_space = GDT::make_discontinuous_lagrange_space(coarse_grid_view, /*order=*/1);
          auto coarse_basis = coarse_space.basis().localize();
          std::vector<XT::LA::IstlDenseVector<double>> interpolated_basis;
          for (auto&& macro_element : elements(coarse_grid_view)) {
            if (dd_grid.subdomain(macro_element) == I) {
              // this is the subdomain we are interested in, create space
               auto subdomain_grid_view = dd_grid.local_grid(macro_element).leaf_view();
               auto subdomain_space =  GDT::make_continuous_lagrange_space(subdomain_grid_view, /*order=*/1);
//              auto subdomain_space = dd_grid.local_space(macro_element);
              coarse_basis->bind(macro_element);
              for (size_t ii = 0; ii < coarse_basis->size(); ++ii)
                interpolated_basis.push_back(Dune::GDT::default_interpolation<XT::LA::IstlDenseVector<double>>(
                                                 coarse_basis->order(),
                                                 [&](const auto& point_in_physical_coordinates, const auto&) {
                                                   const auto point_macro_reference_element =
                                                       macro_element.geometry().local(point_in_physical_coordinates);
                                                   return coarse_basis->evaluate_set(point_macro_reference_element)[ii];
                                                 },
                                                 subdomain_space)
                                                 .dofs()
                                                 .vector());
              break;
            }
          }
          DUNE_THROW_IF(interpolated_basis.size() == 0, InvalidStateException, "This should not happen, I = " << I);
          return interpolated_basis;
        },
        py::call_guard<py::gil_scoped_release>(),
        "dd_grid"_a,
        "I"_a //,
        /*"space_type"_a = "discontinuous_lagrange"*/);
} // PYBIND11_MODULE(dunegdtusercode, ...)
