// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_FUNCTIONALS_L2_BINDINGS_HH
#define DUNE_GDT_FUNCTIONALS_L2_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/grid/grids.bindings.hh>
#include <dune/xt/la/container.bindings.hh>

#include <dune/gdt/spaces/interface.bindings.hh>
#include <dune/gdt/type_traits.hh>

#include "l2.hh"
#include "base.bindings.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <class F, class SP, class V /*= typename XT::LA::Container<typename SP::type::RangeFieldType>::VectorType,
          class GridView = typename Space::GridLayerType,
          class Field = typename Space::RangeFieldType*/>
class L2VolumeVectorFunctional
{
  typedef typename SP::type S;
  static_assert(is_space<S>::value, "");

public:
  typedef GDT::L2VolumeVectorFunctional<F, S, V> type;
  typedef pybind11::class_<type> bound_type;

  static bound_type bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = XT::Common::to_camel_case("l2_volume_vector_functional_" + space_name<SP>::value() + "_"
                                                     + XT::LA::bindings::container_name<V>::value());

    auto c = VectorFunctionalBase<type>::bind(m, ClassName);

    m.def(std::string("make_l2_volume_vector_functional_" + XT::LA::bindings::container_name<V>::value()).c_str(),
          [](const F& function, const S& space, const size_t over_integrate) {
            return make_l2_volume_vector_functional<V>(function, space, over_integrate).release(); // <-            b.c.
          }, //    L2VolumeVectorFunctional is not movable, returning the raw pointer lets pybind11 correctly manage the
          "function"_a, //                                                                                        memory
          "space"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());

    m.def("make_l2_volume_vector_functional",
          [](const F& function, V& vector, const S& space, const size_t over_integrate) {
            return make_l2_volume_vector_functional(function, vector, space, over_integrate).release(); //       <- s.a.
          },
          "function"_a,
          "vector"_a,
          "space"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>());

    return c;
  } // ... bind(...)
}; // class L2VolumeVectorFunctional


template <class F,
          class SP,
          class V /*= typename XT::LA::Container<typename SP::type::RangeFieldType>::VectorType*/,
          class GL /*= typename SP::type::GridLayerType,
          class Field = typename SP::type::RangeFieldType*/>
class L2FaceVectorFunctional
{
  typedef typename SP::type S;
  static_assert(is_space<S>::value, "");

public:
  typedef GDT::L2FaceVectorFunctional<F, S, V> type;
  typedef pybind11::class_<type> bound_type;

  static bound_type bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const std::string class_name = "l2_face_vector_functional";
    const auto ClassName = XT::Common::to_camel_case(class_name + "_" + space_name<SP>::value() + "_"
                                                     + XT::LA::bindings::container_name<V>::value());

    auto c = VectorFunctionalBase<type>::bind(m, ClassName);

    m.def(std::string("make_" + class_name + "_" + XT::LA::bindings::container_name<V>::value()).c_str(),
          [](const F& function, const S& space, const size_t over_integrate) {
            return make_l2_face_vector_functional<V>(function, space, over_integrate).release(); //              <- b.c.
          }, //      L2FaceVectorFunctional is not movable, returning the raw pointer lets pybind11 correctly manage the
          "function"_a, //                                                                                        memory
          "space"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());
    m.def(std::string("make_" + class_name + "_" + XT::LA::bindings::container_name<V>::value()).c_str(),
          [](const F& function,
             const S& space,
             const XT::Grid::ApplyOn::WhichIntersection<GL>& which_intersections,
             const size_t over_integrate) {
            return make_l2_face_vector_functional<V>(function, space, over_integrate, which_intersections.copy())
                .release(); //                                                                                   <- s.a.
          },
          "function"_a,
          "space"_a,
          "which_intersections"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());

    m.def(std::string("make_" + class_name).c_str(),
          [](const F& function, V& vector, const S& space, const size_t over_integrate) {
            return make_l2_face_vector_functional(function, vector, space, over_integrate).release(); //         <- s.a.
          },
          "function"_a,
          "vector"_a,
          "space"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>());
    m.def(std::string("make_" + class_name).c_str(),
          [](const F& function,
             V& vector,
             const S& space,
             const XT::Grid::ApplyOn::WhichIntersection<GL>& which_intersections,
             const size_t over_integrate) {
            return make_l2_face_vector_functional(function, vector, space, over_integrate, which_intersections.copy())
                .release(); //                                                                                   <- s.a.
          },
          "function"_a,
          "vector"_a,
          "space"_a,
          "which_intersections"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>());

    return c;
  } // ... bind(...)
}; // XT::LA::bindings::container_name<V>::value()


} // namespace bindings
} // namespace GDT
} // namespace Dune


#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_FUNCTIONALS_L2_BINDINGS_HH
