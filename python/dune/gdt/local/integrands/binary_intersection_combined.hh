// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_GDT_LOCAL_INTEGRANDS_BINARY_INTERSECTION_COMBINED_HH
#define PYTHON_DUNE_GDT_LOCAL_INTEGRANDS_BINARY_INTERSECTION_COMBINED_HH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/gridprovider/provider.hh>

#include <dune/gdt/local/integrands/combined.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

#include "binary_intersection_interface.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <class I,
          size_t t_r = 1,
          size_t t_rC = 1,
          class TF = double,
          class F = double,
          size_t a_r = t_r,
          size_t a_rC = t_rC,
          class AF = TF>
class LocalBinaryIntersectionIntegrandSum
  : protected LocalBinaryIntersectionIntegrandInterface<I, t_r, t_rC, TF, F, a_r, a_rC, AF>
{
  using BaseType = LocalBinaryIntersectionIntegrandInterface<I, t_r, t_rC, TF, F, a_r, a_rC, AF>;
  using typename BaseType::G;

public:
  using type = GDT::LocalBinaryIntersectionIntegrandSum<I, t_r, t_rC, TF, F, a_r, a_rC, AF>;
  using bound_type = pybind11::class_<type, typename BaseType::type>;

  static bound_type bind(pybind11::module& m,
                         const std::string& class_id = "local_binary_intersection_integrand_sum",
                         const std::string& grid_id = XT::Grid::bindings::grid_name<G>::value(),
                         const std::string& layer_id = "")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto class_name = class_id + BaseType::id(grid_id, layer_id);
    const auto ClassName = XT::Common::to_camel_case(class_name);
    bound_type c(m, ClassName.c_str(), ClassName.c_str());

    return c;
  } // ... bind(...)
}; // class LocalBinaryIntersectionIntegrandSum


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // PYTHON_DUNE_GDT_LOCAL_INTEGRANDS_BINARY_INTERSECTION_COMBINED_HH
