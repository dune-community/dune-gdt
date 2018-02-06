// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_LOCAL_DIFFUSIVE_FLUX_ESTIMATION_OPERATOR_BINDINGS_HH
#define DUNE_GDT_LOCAL_DIFFUSIVE_FLUX_ESTIMATION_OPERATOR_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/functions/interfaces/localizable-function.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/local/integrands/ESV2007.hh>
#include <dune/gdt/local/operators/integrals.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class G>
class LocalDiffusiveFluxEstimationOperator
{
  static_assert(XT::Grid::is_grid<G>::value, "");
  typedef typename G::template Codim<0>::Entity E;
  typedef typename G::ctype D;
  static const constexpr size_t d = G::dimension;
  typedef double R;

  typedef XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1> ScalarFunction;
  typedef XT::Functions::LocalizableFunctionInterface<E, D, d, R, d> VectorFunction;
  typedef XT::Functions::LocalizableFunctionInterface<E, D, d, R, d, d> TensorFunction;
  typedef typename ScalarFunction::LocalfunctionType ScalarLocalFunctionType;
  typedef GDT::LocalVolumeTwoFormInterface<ScalarLocalFunctionType> InterfaceType;

public:
  typedef LocalVolumeIntegralOperator<LocalDiffusiveFluxEstimateESV2007Integrand<ScalarFunction,
                                                                                 VectorFunction,
                                                                                 TensorFunction>,
                                      ScalarLocalFunctionType>
      type;
  typedef pybind11::class_<type, InterfaceType> bound_type;

  static bound_type bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;
    using XT::Common::to_string;

    // bind interface, guard since we might not be the first to do so for this combination
    try {
      const auto InterfaceName = XT::Common::to_camel_case(
          "local_volume_two_form_interface_" + XT::Grid::bindings::grid_name<G>::value() + "_to_1x1");
      py::class_<InterfaceType>(m, InterfaceName.c_str(), "LocalVolumeTwoFormInterface");
    } catch (std::runtime_error&) {
    }

    const std::string class_name = "local_diffusive_flux_estimation_esv2007_operator";
    const auto ClassName =
        XT::Common::to_camel_case(class_name + "_" + XT::Grid::bindings::grid_name<G>::value() + "_to_1x1");

    bound_type c(m, ClassName.c_str(), "LocalVolumeIntegralOperator<LocalDiffusiveFluxEstimateESV2007Integrand<...>>");

    m.def(("make_" + class_name + "_to_1x1").c_str(),
          [](const ScalarFunction& diffusion_factor,
             const TensorFunction& diffusion_tensor,
             const VectorFunction& diffusive_flux,
             const size_t over_integrate) {
            return type(over_integrate, diffusion_factor, diffusion_tensor, diffusive_flux);
          },
          "diffusion_factor"_a,
          "diffusion_tensor"_a,
          "diffusive_flux"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>());

    return c;
  } // ... bind(...)
}; // class LocalDiffusiveFluxEstimationOperator


} // namespace bindings
} // namespace GDT
} // namespace Dune


#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_LOCAL_DIFFUSIVE_FLUX_ESTIMATION_OPERATOR_BINDINGS_HH
