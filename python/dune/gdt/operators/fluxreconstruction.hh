// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef PYTHON_DUNE_GDT_OPERATORS_FLUXRECONSTRUCTION_BINDINGS_HH
#define PYTHON_DUNE_GDT_OPERATORS_FLUXRECONSTRUCTION_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/la/container.hh>

#include <dune/gdt/spaces.hh>

#include <dune/gdt/operators/fluxreconstruction.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class G,
          SpaceType space_type,
          Backends space_backend,
          XT::Grid::Layers space_layer_type,
          int p,
          class R,
          size_t r,
          XT::LA::Backends la_backend /*,
          XT::Grid::Layers layer_type = space_layer_type,
          XT::Grid::Backends layer_backend =
              SpaceProvider<G, space_layer_type, space_type, space_backend, p, R, r>::layer_backend*/>
class DiffusiveFluxReconstructionOperator
{
  typedef typename SpaceProvider<G, space_layer_type, space_type, space_backend, p, R, r>::type S;
  typedef typename S::GridLayerType GL;
  typedef typename XT::LA::Container<R, la_backend>::VectorType V;
  typedef typename S::EntityType E;
  typedef typename S::DomainFieldType D;
  static const size_t d = S::dimDomain;
  typedef XT::Functions::GridFunctionInterface<E, D, d, R, 1> ScalarFunctionType;
  typedef XT::Functions::GridFunctionInterface<E, D, d, R, d, d> TensorFunctionType;

public:
  static void bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    //    m.def("apply_diffusive_flux_reconstruction_operator",
    //          [](const ScalarFunctionType& diffusion_factor,
    //             const ScalarFunctionType& source,
    //             DiscreteFunction<S, V>& range) {
    //            py::gil_scoped_release DUNE_UNUSED(release);
    //            GDT::DiffusiveFluxReconstructionOperator<GL,
    //                                                     ScalarFunctionType,
    //                                                     void,
    //                                                     LocalEllipticIpdgIntegrands::Method::swipdg_affine_factor>(
    //                range.space().grid_layer(), diffusion_factor)
    //                .apply(source, range);
    //          },
    //          "diffusion_tensor"_a,
    //          "source"_a,
    //          "range"_a);
    m.def("apply_diffusive_flux_reconstruction_operator",
          [](const ScalarFunctionType& diffusion_factor,
             const TensorFunctionType& diffusion_tensor,
             const ScalarFunctionType& source,
             DiscreteFunction<S, V>& range) {
            py::gil_scoped_release DUNE_UNUSED(release);
            GDT::DiffusiveFluxReconstructionOperator<GL,
                                                     ScalarFunctionType,
                                                     TensorFunctionType,
                                                     LocalEllipticIpdgIntegrands::Method::swipdg_affine_factor>(
                range.space().grid_layer(), diffusion_factor, diffusion_tensor)
                .apply(source, range);
          },
          "diffusion_factor"_a,
          "diffusion_tensor"_a,
          "source"_a,
          "range"_a);
  }
}; // class DiffusiveFluxReconstructionOperator


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
#endif // PYTHON_DUNE_GDT_OPERATORS_FLUXRECONSTRUCTION_BINDINGS_HH
