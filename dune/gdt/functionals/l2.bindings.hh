// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_FUNCTIONALS_L2_PBH
#define DUNE_GDT_FUNCTIONALS_L2_PBH
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include "l2.hh"
#include "base.pbh"

namespace Dune {
namespace GDT {


template <class FunctionType,
          class Space,
          class Vector = typename XT::LA::Container<typename Space::RangeFieldType>::VectorType,
          class GridView = typename Space::GridViewType,
          class Field = typename Space::RangeFieldType>
pybind11::class_<L2VolumeVectorFunctional<FunctionType, Space, Vector, GridView, Field>>
bind_l2_volume_vector_functional(pybind11::module& m, const std::string& space_id, const std::string& la_id)
{
  static_assert(std::is_same<GridView, typename Space::GridViewType>::value, "Not tested yet!");

  namespace py = pybind11;
  using namespace pybind11::literals;

  typedef L2VolumeVectorFunctional<FunctionType, Space, Vector, GridView, Field> C;

  auto c = bind_vector_functional<C>(m, "L2VolumeVectorFunctional__" + la_id + "__" + space_id);

  m.def(std::string("make_l2_volume_vector_functional__" + la_id).c_str(),
        [](const FunctionType& function, const Space& space, const size_t over_integrate) {
          return make_l2_volume_vector_functional<Vector>(function, space, over_integrate).release(); // b.c.
        }, //      L2VolumeVectorFunctional is not movable, returning the raw pointer lets pybind11 correctly
        "function"_a, //                                                                    manage the memory
        "space"_a,
        "over_integrate"_a = 0,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>());

  m.def("make_l2_volume_vector_functional",
        [](const FunctionType& function, Vector& vector, const Space& space, const size_t over_integrate) {
          return make_l2_volume_vector_functional(function, vector, space, over_integrate).release(); // s.a.
        },
        "function"_a,
        "vector"_a,
        "space"_a,
        "over_integrate"_a = 0,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>(),
        py::keep_alive<0, 3>());

  return c;

} // ... bind_l2_volume_vector_functional(...)


template <class FunctionType,
          class Space,
          class Vector = typename XT::LA::Container<typename Space::RangeFieldType>::VectorType,
          class GridView = typename Space::GridViewType,
          class Field = typename Space::RangeFieldType>
pybind11::class_<L2FaceVectorFunctional<FunctionType, Space, Vector, GridView, Field>>
bind_l2_face_vector_functional(pybind11::module& m, const std::string& space_id, const std::string& la_id)
{
  static_assert(std::is_same<GridView, typename Space::GridViewType>::value, "Not tested yet!");

  namespace py = pybind11;
  using namespace pybind11::literals;

  typedef L2FaceVectorFunctional<FunctionType, Space, Vector, GridView, Field> C;

  auto c = bind_vector_functional<C>(m, "L2FaceVectorFunctional__" + la_id + "__" + space_id);

  m.def(std::string("make_l2_face_vector_functional__" + la_id).c_str(),
        [](const FunctionType& function, const Space& space, const size_t over_integrate) {
          return make_l2_face_vector_functional<Vector>(function, space, over_integrate).release(); // b.c.
        }, //      L2FaceVectorFunctional is not movable, returning the raw pointer lets pybind11 correctly
        "function"_a, //                                                                    manage the memory
        "space"_a,
        "over_integrate"_a = 0,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>());
  m.def(std::string("make_l2_face_vector_functional__" + la_id).c_str(),
        [](const FunctionType& function,
           const Space& space,
           const XT::Grid::ApplyOn::WhichIntersection<GridView>& which_intersections,
           const size_t over_integrate) {
          return make_l2_face_vector_functional<Vector>(function, space, over_integrate, which_intersections.copy())
              .release(); // b.c.
        }, //      L2FaceVectorFunctional is not movable, returning the raw pointer lets pybind11 correctly
        "function"_a, //                                                                    manage the memory
        "space"_a,
        "which_intersections"_a,
        "over_integrate"_a = 0,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>());

  m.def("make_l2_face_vector_functional",
        [](const FunctionType& function, Vector& vector, const Space& space, const size_t over_integrate) {
          return make_l2_face_vector_functional(function, vector, space, over_integrate).release(); // s.a.
        },
        "function"_a,
        "vector"_a,
        "space"_a,
        "over_integrate"_a = 0,
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>(),
        py::keep_alive<0, 3>());
  m.def("make_l2_face_vector_functional",
        [](const FunctionType& function,
           Vector& vector,
           const Space& space,
           const XT::Grid::ApplyOn::WhichIntersection<GridView>& which_intersections,
           const size_t over_integrate) {
          return make_l2_face_vector_functional(function, vector, space, over_integrate, which_intersections.copy())
              .release(); // s.a.
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
} // ... bind_l2_face_vector_functional(...)


} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_FUNCTIONALS_L2_PBH
