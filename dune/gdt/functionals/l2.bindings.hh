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

#include <dune/gdt/spaces.bindings.hh>
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


// begin: this is what we need for the lib

#define _DUNE_GDT_FUNCTIONALS_L2_BIND_LIB(_prefix, _d, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la)        \
  _prefix class Dune::GDT::bindings::                                                                                  \
      L2VolumeVectorFunctional<Dune::XT::Functions::                                                                   \
                                   LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<                      \
                                                                    typename Dune::XT::Grid::                          \
                                                                        Layer<_GRID,                                   \
                                                                              Dune::XT::Grid::Layers::_layer,          \
                                                                              Dune::XT::Grid::Backends::_g_backend,    \
                                                                              Dune::XT::Grid::DD::                     \
                                                                                  SubdomainGrid<_GRID>>::type>,        \
                                                                double,                                                \
                                                                _d,                                                    \
                                                                double,                                                \
                                                                1,                                                     \
                                                                1>,                                                    \
                               Dune::GDT::SpaceProvider<_GRID,                                                         \
                                                        Dune::XT::Grid::Layers::_layer,                                \
                                                        Dune::GDT::SpaceType::_s_type,                                 \
                                                        Dune::GDT::Backends::_s_backend,                               \
                                                        _p,                                                            \
                                                        double,                                                        \
                                                        1,                                                             \
                                                        1>,                                                            \
                               typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::VectorType>;     \
  _prefix class Dune::GDT::bindings::                                                                                  \
      L2FaceVectorFunctional<Dune::XT::Functions::                                                                     \
                                 LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<                        \
                                                                  typename Dune::XT::Grid::                            \
                                                                      Layer<_GRID,                                     \
                                                                            Dune::XT::Grid::Layers::_layer,            \
                                                                            Dune::XT::Grid::Backends::_g_backend,      \
                                                                            Dune::XT::Grid::DD::                       \
                                                                                SubdomainGrid<_GRID>>::type>,          \
                                                              double,                                                  \
                                                              _d,                                                      \
                                                              double,                                                  \
                                                              1,                                                       \
                                                              1>,                                                      \
                             Dune::GDT::SpaceProvider<_GRID,                                                           \
                                                      Dune::XT::Grid::Layers::_layer,                                  \
                                                      Dune::GDT::SpaceType::_s_type,                                   \
                                                      Dune::GDT::Backends::_s_backend,                                 \
                                                      _p,                                                              \
                                                      double,                                                          \
                                                      1,                                                               \
                                                      1>,                                                              \
                             typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::VectorType,        \
                             typename Dune::XT::Grid::Layer<_GRID,                                                     \
                                                            Dune::XT::Grid::Layers::_layer,                            \
                                                            Dune::XT::Grid::Backends::_g_backend,                      \
                                                            Dune::XT::Grid::DD::SubdomainGrid<_GRID>>::type>


#if HAVE_DUNE_ALUGRID
#define DUNE_GDT_FUNCTIONALS_L2_BIND_LIB_ALU(_prefix, _layer, _g_backend, _s_type, _s_backend, _p, _la)                \
  _DUNE_GDT_FUNCTIONALS_L2_BIND_LIB(                                                                                   \
      _prefix, 2, ALU_2D_SIMPLEX_CONFORMING, _layer, _g_backend, _s_type, _s_backend, _p, _la)
#else
#define DUNE_GDT_FUNCTIONALS_L2_BIND_LIB_ALU(_prefix, _layer, _g_backend, _s_type, _s_backend, _p, _la)
#endif

#define DUNE_GDT_FUNCTIONALS_L2_BIND_LIB_YASP(_prefix, _layer, _g_backend, _s_type, _s_backend, _p, _la)               \
  _DUNE_GDT_FUNCTIONALS_L2_BIND_LIB(                                                                                   \
      _prefix, 1, YASP_1D_EQUIDISTANT_OFFSET, _layer, _g_backend, _s_type, _s_backend, _p, _la);                       \
  _DUNE_GDT_FUNCTIONALS_L2_BIND_LIB(                                                                                   \
      _prefix, 2, YASP_2D_EQUIDISTANT_OFFSET, _layer, _g_backend, _s_type, _s_backend, _p, _la)

// alu_fem_istl.cc
#if HAVE_DUNE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_ISTL
DUNE_GDT_FUNCTIONALS_L2_BIND_LIB_ALU(extern template, leaf, part, cg, fem, 1, istl_sparse);
DUNE_GDT_FUNCTIONALS_L2_BIND_LIB_ALU(extern template, level, part, cg, fem, 1, istl_sparse);
DUNE_GDT_FUNCTIONALS_L2_BIND_LIB_ALU(extern template, dd_subdomain, part, cg, fem, 1, istl_sparse);
#endif

// yasp_fem_istl.cc
#if HAVE_DUNE_FEM && HAVE_DUNE_ISTL
// DUNE_GDT_FUNCTIONALS_L2_BIND_LIB_YASP(extern template, leaf, part, cg, fem, 1, istl_sparse);
// DUNE_GDT_FUNCTIONALS_L2_BIND_LIB_YASP(extern template, level, part, cg, fem, 1, istl_sparse);
// DUNE_GDT_FUNCTIONALS_L2_BIND_LIB_YASP(extern template, dd_subdomain, part, cg, fem, 1, istl_sparse);
#endif

// end: this is what we need for the lib


// begin: this is what we need for the .so

#define _DUNE_GDT_FUNCTIONALS_L2_BIND(_m, _d, _GRID, _layer, _g_backend, _s_type, _s_backend, _p, _la)                 \
  Dune::GDT::bindings::                                                                                                \
      L2VolumeVectorFunctional<Dune::XT::Functions::                                                                   \
                                   LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<                      \
                                                                    typename Dune::XT::Grid::                          \
                                                                        Layer<_GRID,                                   \
                                                                              Dune::XT::Grid::Layers::_layer,          \
                                                                              Dune::XT::Grid::Backends::_g_backend,    \
                                                                              Dune::XT::Grid::DD::                     \
                                                                                  SubdomainGrid<_GRID>>::type>,        \
                                                                double,                                                \
                                                                _d,                                                    \
                                                                double,                                                \
                                                                1,                                                     \
                                                                1>,                                                    \
                               Dune::GDT::SpaceProvider<_GRID,                                                         \
                                                        Dune::XT::Grid::Layers::_layer,                                \
                                                        Dune::GDT::SpaceType::_s_type,                                 \
                                                        Dune::GDT::Backends::_s_backend,                               \
                                                        _p,                                                            \
                                                        double,                                                        \
                                                        1,                                                             \
                                                        1>,                                                            \
                               typename Dune::XT::LA::Container<double,                                                \
                                                                Dune::XT::LA::Backends::_la>::VectorType>::bind(_m);   \
  Dune::GDT::bindings::                                                                                                \
      L2FaceVectorFunctional<Dune::XT::Functions::                                                                     \
                                 LocalizableFunctionInterface<Dune::XT::Grid::extract_entity_t<                        \
                                                                  typename Dune::XT::Grid::                            \
                                                                      Layer<_GRID,                                     \
                                                                            Dune::XT::Grid::Layers::_layer,            \
                                                                            Dune::XT::Grid::Backends::_g_backend,      \
                                                                            Dune::XT::Grid::DD::                       \
                                                                                SubdomainGrid<_GRID>>::type>,          \
                                                              double,                                                  \
                                                              _d,                                                      \
                                                              double,                                                  \
                                                              1,                                                       \
                                                              1>,                                                      \
                             Dune::GDT::SpaceProvider<_GRID,                                                           \
                                                      Dune::XT::Grid::Layers::_layer,                                  \
                                                      Dune::GDT::SpaceType::_s_type,                                   \
                                                      Dune::GDT::Backends::_s_backend,                                 \
                                                      _p,                                                              \
                                                      double,                                                          \
                                                      1,                                                               \
                                                      1>,                                                              \
                             typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::VectorType,        \
                             typename Dune::XT::Grid::Layer<_GRID,                                                     \
                                                            Dune::XT::Grid::Layers::_layer,                            \
                                                            Dune::XT::Grid::Backends::_g_backend,                      \
                                                            Dune::XT::Grid::DD::SubdomainGrid<_GRID>>::type>::bind(_m)


#if HAVE_DUNE_ALUGRID
#define DUNE_GDT_FUNCTIONALS_L2_BIND_ALU(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la)                         \
  _DUNE_GDT_FUNCTIONALS_L2_BIND(_m, 2, ALU_2D_SIMPLEX_CONFORMING, _layer, _g_backend, _s_type, _s_backend, _p, _la)
#else
#define DUNE_GDT_FUNCTIONALS_L2_BIND_ALU(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la)
#endif

#define DUNE_GDT_FUNCTIONALS_L2_BIND_YASP(_m, _layer, _g_backend, _s_type, _s_backend, _p, _la)                        \
  _DUNE_GDT_FUNCTIONALS_L2_BIND(_m, 1, YASP_1D_EQUIDISTANT_OFFSET, _layer, _g_backend, _s_type, _s_backend, _p, _la);  \
  _DUNE_GDT_FUNCTIONALS_L2_BIND(_m, 2, YASP_2D_EQUIDISTANT_OFFSET, _layer, _g_backend, _s_type, _s_backend, _p, _la)


// end: this is what we need for the .so


#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_FUNCTIONALS_L2_BINDINGS_HH
