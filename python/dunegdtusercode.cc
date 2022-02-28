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

namespace py = pybind11;
using namespace Dune;


template <class V>
V& get_vec_ref(py::handle list_element, const bool recurse = true)
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
} // PYBIND11_MODULE(dunegdtusercode, ...)
