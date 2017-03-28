// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BINDINGS_HH
#define DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/grid/grids.bindings.hh>
#include <dune/xt/la/type_traits.hh>
#include <dune/xt/la/container.bindings.hh>

#include <dune/gdt/type_traits.hh>
#include <dune/gdt/spaces.bindings.hh>

#include "default.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <class SP, class V>
class ConstDiscreteFunction
{
  typedef typename SP::type S;
  static_assert(is_space<S>::value, "");
  static_assert(XT::LA::is_vector<V>::value, "");

public:
  typedef GDT::ConstDiscreteFunction<S, V> type;

private:
  typedef XT::Functions::LocalizableFunctionInterface<typename S::EntityType,
                                                      typename S::DomainFieldType,
                                                      S::dimDomain,
                                                      typename S::RangeFieldType,
                                                      S::dimRange,
                                                      S::dimRangeCols>
      BaseType;

public:
  typedef pybind11::class_<type, BaseType> bound_type;

  static bound_type bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = XT::Common::to_camel_case("const_discrete_function_" + space_name<SP>::value() + "_"
                                                     + XT::LA::bindings::container_name<V>::value());

    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def(py::init<const S&, V&, const std::string>(),
          "space"_a,
          "vector"_a,
          "name"_a = "gdt.constdiscretefunction",
          py::keep_alive<1, 2>(),
          py::keep_alive<1, 3>());
    c.def("space", [](type& self) { return self.space(); });
    c.def("vector", [](type& self) { return self.vector(); });
    c.def("visualize",
          [](type& self, const std::string filename, const bool subsampling) {
            return self.visualize(filename, subsampling);
          },
          "filename"_a,
          "subsampling"_a = (S::polOrder > 1));

    m.def(std::string("make_const_discrete_function").c_str(),
          [](const S& space, V& vector, const std::string& name) {
            return make_const_discrete_function(space, vector, name);
          },
          "space"_a,
          "vector"_a,
          "name"_a = "gdt.constdiscretefunction",
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());

    return c;
  } // ... bind(...)
}; // class ConstDiscreteFunction


template <class SP, class V>
class DiscreteFunction
{
  typedef typename SP::type S;
  static_assert(is_space<S>::value, "");
  static_assert(XT::LA::is_vector<V>::value, "");

public:
  typedef GDT::DiscreteFunction<S, V> type;

private:
  typedef XT::Functions::LocalizableFunctionInterface<typename S::EntityType,
                                                      typename S::DomainFieldType,
                                                      S::dimDomain,
                                                      typename S::RangeFieldType,
                                                      S::dimRange,
                                                      S::dimRangeCols>
      BaseType;

public:
  typedef pybind11::class_<type, BaseType> bound_type;

  static bound_type bind(pybind11::module& m, const std::string& space_id, const std::string& la_id)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = XT::Common::to_camel_case("discrete_function_" + space_name<SP>::value() + "_"
                                                     + XT::LA::bindings::container_name<V>::value());

    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    c.def(py::init<const S&, V&, const std::string>(),
          "space"_a,
          "vector"_a,
          "name"_a = "gdt.discretefunction",
          py::keep_alive<1, 2>(),
          py::keep_alive<1, 3>());
    c.def("space", [](type& self) { return self.space(); });
    c.def("vector", [](type& self) { return self.vector(); });
    c.def("visualize",
          [](type& self, const std::string filename, const bool subsampling) {
            return self.visualize(filename, subsampling);
          },
          "filename"_a,
          "subsampling"_a = (S::polOrder > 1));

    m.def(
        std::string("make_discrete_function").c_str(),
        [](const S& space, V& vector, const std::string& name) { return make_discrete_function(space, vector, name); },
        "space"_a,
        "vector"_a,
        "name"_a = "gdt.discretefunction",
        py::keep_alive<0, 1>(),
        py::keep_alive<0, 2>());

    return c;
  } // ... bind(...)

}; // class DiscreteFunction


} // namespace bindings
} // namespace GDT
} // namespace Dune


// begin: this is what we need for the .so

#if HAVE_DUNE_FEM

// * fem backend
//   - cg spaces

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_CG(_m, _GRID, _layer, _p, _r, _rC, _la)                            \
  Dune::GDT::bindings::ConstDiscreteFunction<Dune::GDT::CgSpaceProvider<_GRID,                                         \
                                                                        Dune::XT::Grid::Layers::_layer,                \
                                                                        Dune::GDT::ChooseSpaceBackend::fem,            \
                                                                        _p,                                            \
                                                                        double,                                        \
                                                                        _r,                                            \
                                                                        _rC>,                                          \
                                             typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::   \
                                                 VectorType>::bind(_m)

#if HAVE_ALBERTA
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_CG_ALBERTA(_m, _layer, _p, _la)                                    \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_CG(_m, ALBERTA_2D, _layer, _p, 1, 1, _la)
#else
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_CG_ALBERTA(_m, _p, _layer, _la)
#endif

#if HAVE_DUNE_ALUGRID
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_CG_ALU(_m, _layer, _p, _la)                                        \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_CG(_m, ALU_2D_SIMPLEX_CONFORMING, _layer, _p, 1, 1, _la)
#else
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_CG_ALU(_m, _layer, _p, _la)
#endif

#if HAVE_DUNE_UGGRID || HAVE_UG
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_CG_UG(_m, _layer, _p, _la)                                         \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_CG(_m, UG_2D, _layer, _p, 1, 1, _la)
#else
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_CG_UG(_m, _layer, _p, _la)
#endif

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_CG_YASP(_m, _layer, _p, _la)                                       \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_CG(_m, YASP_2D_EQUIDISTANT_OFFSET, _layer, _p, 1, 1, _la)

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_CG_ALL(_m, _layer, _p, _la)                                        \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_CG_ALBERTA(_m, _layer, _p, _la);                                         \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_CG_ALU(_m, _layer, _p, _la);                                             \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_CG_UG(_m, _layer, _p, _la);                                              \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_CG_YASP(_m, _layer, _p, _la)

//   - dg spaces

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG(_m, _GRID, _layer, _p, _r, _rC, _la)                            \
  Dune::GDT::bindings::ConstDiscreteFunction<Dune::GDT::DgSpaceProvider<_GRID,                                         \
                                                                        Dune::XT::Grid::Layers::_layer,                \
                                                                        Dune::GDT::ChooseSpaceBackend::fem,            \
                                                                        _p,                                            \
                                                                        double,                                        \
                                                                        _r,                                            \
                                                                        _rC>,                                          \
                                             typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::   \
                                                 VectorType>::bind(_m)

#if HAVE_ALBERTA
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG_ALBERTA(_m, _layer, _p, _la)                                    \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG(_m, ALBERTA_2D, _layer, _p, 1, 1, _la)
#else
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG_ALBERTA(_m, _p, _layer, _la)
#endif

#if HAVE_DUNE_ALUGRID
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG_ALU(_m, _layer, _p, _la)                                        \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG(_m, ALU_2D_SIMPLEX_CONFORMING, _layer, _p, 1, 1, _la)
#else
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG_ALU(_m, _layer, _p, _la)
#endif

//#if HAVE_DUNE_UGGRID || HAVE_UG // <- does not work
//#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG_UG(_m, _layer, _p, _la)
//  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG(_m, UG_2D, _layer, _p, 1, 1, _la)
//#else
//#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG_UG(_m, _layer, _p, _la)
//#endif

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG_YASP(_m, _layer, _p, _la)                                       \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG(_m, YASP_2D_EQUIDISTANT_OFFSET, _layer, _p, 1, 1, _la)

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG_ALL(_m, _layer, _p, _la)                                        \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG_ALBERTA(_m, _layer, _p, _la);                                         \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG_ALU(_m, _layer, _p, _la);                                             \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG_YASP(_m, _layer, _p, _la)
//  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG_UG(_m, _layer, _p, _la); // <- does not work

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM(_m, _la)                                                           \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_CG_ALL(_m, dd_subdomain, 1, _la);                                        \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_CG_ALL(_m, leaf, 1, _la);                                                \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_CG_ALL(_m, level, 1, _la);                                               \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG_ALL(_m, dd_subdomain, 1, _la);                                        \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG_ALL(_m, leaf, 1, _la);                                                \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG_ALL(_m, level, 1, _la);                                               \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG_ALL(_m, dd_subdomain, 2, _la);                                        \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG_ALL(_m, leaf, 2, _la);                                                \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM_DG_ALL(_m, level, 2, _la)

#else // HAVE_DUNE_FEM
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM(_m, _la)
#endif

// * gdt backend
//   - fv spaces

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT_FV(_m, _GRID, _layer, _r, _rC, _la)                                \
  Dune::GDT::bindings::ConstDiscreteFunction<Dune::GDT::FvSpaceProvider<_GRID,                                         \
                                                                        Dune::XT::Grid::Layers::_layer,                \
                                                                        Dune::GDT::ChooseSpaceBackend::gdt,            \
                                                                        double,                                        \
                                                                        _r,                                            \
                                                                        _rC>,                                          \
                                             typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::   \
                                                 VectorType>::bind(_m)

#if HAVE_ALBERTA
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT_FV_ALBERTA(_m, _layer, _la)                                        \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT_FV(_m, ALBERTA_2D, _layer, 1, 1, _la)
#else
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT_FV_ALBERTA(_m, _layer, _la)
#endif

#if HAVE_DUNE_ALUGRID
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT_FV_ALU(_m, _layer, _la)                                            \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT_FV(_m, ALU_2D_SIMPLEX_CONFORMING, _layer, 1, 1, _la)
#else
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT_FV_ALU(_m, _layer, _la)
#endif

#if HAVE_DUNE_UGGRID || HAVE_UG
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT_FV_UG(_m, _layer, _la)                                             \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT_FV(_m, UG_2D, _layer, 1, 1, _la)
#else
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT_FV_UG(_m, _layer, _la)
#endif

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT_FV_YASP(_m, _layer, _la)                                           \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT_FV(_m, YASP_2D_EQUIDISTANT_OFFSET, _layer, 1, 1, _la)

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT_FV_ALL(_m, _layer, _la)                                            \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT_FV_ALBERTA(_m, _layer, _la);                                             \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT_FV_ALU(_m, _layer, _la);                                                 \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT_FV_UG(_m, _layer, _la);                                                  \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT_FV_YASP(_m, _layer, _la)

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT(_m, _la)                                                           \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT_FV_ALL(_m, leaf, _la);                                                   \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT_FV_ALL(_m, level, _la)


#if HAVE_DUNE_PDELAB

// * pdelab backend
//   - cg spaces

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB_CG(_m, _GRID, _layer, _p, _r, _rC, _la)                         \
  Dune::GDT::bindings::ConstDiscreteFunction<Dune::GDT::CgSpaceProvider<_GRID,                                         \
                                                                        Dune::XT::Grid::Layers::_layer,                \
                                                                        Dune::GDT::ChooseSpaceBackend::pdelab,         \
                                                                        _p,                                            \
                                                                        double,                                        \
                                                                        _r,                                            \
                                                                        _rC>,                                          \
                                             typename Dune::XT::LA::Container<double, Dune::XT::LA::Backends::_la>::   \
                                                 VectorType>::bind(_m)

#if HAVE_ALBERTA
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB_CG_ALBERTA(_m, _layer, _p, _la)                                 \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB_CG(_m, ALBERTA_2D, _layer, _p, 1, 1, _la)
#else
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB_CG_ALBERTA(_m, _p, _layer, _la)
#endif

#if HAVE_DUNE_ALUGRID
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB_CG_ALU(_m, _layer, _p, _la)                                     \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB_CG(_m, ALU_2D_SIMPLEX_CONFORMING, _layer, _p, 1, 1, _la)
#else
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB_CG_ALU(_m, _layer, _p, _la)
#endif

//#if HAVE_DUNE_UGGRID || HAVE_UG // <- does not work
//#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB_CG_UG(_m, _layer, _p, _la)
//  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB_CG(_m, UG_2D, _layer, _p, 1, 1, _la)
//#else
//#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB_CG_UG(_m, _layer, _p, _la)
//#endif

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB_CG_YASP(_m, _layer, _p, _la)                                    \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB_CG(_m, YASP_2D_EQUIDISTANT_OFFSET, _layer, _p, 1, 1, _la)

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB_CG_ALL(_m, _layer, _p, _la)                                     \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB_CG_ALBERTA(_m, _layer, _p, _la);                                      \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB_CG_ALU(_m, _layer, _p, _la);                                          \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB_CG_YASP(_m, _layer, _p, _la)
//_DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB_CG_UG(_m, _layer, _p, _la); // <- does not work

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB(_m, _la)                                                        \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB_CG_ALL(_m, leaf, 1, _la);                                             \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB_CG_ALL(_m, level, 1, _la);

#else // HAVE_DUNE_PDELAB
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB(_m, _la)
#endif

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_LA(_m, _la)                                                            \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM(_m, _la);                                                                \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT(_m, _la);                                                                \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB(_m, _la)

#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_COMMON(_m) _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_LA(_m, common_dense)

#if HAVE_EIGEN
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_EIGEN(_m) _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_LA(_m, eigen_dense)
#else
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_EIGEN(_m)
#endif

#if HAVE_DUNE_ISTL
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ISTL(_m) _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_LA(_m, istl_dense)
#else
#define _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ISTL(_m)
#endif

#define DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND(_m)                                                                     \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_COMMON(_m);                                                                  \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_EIGEN(_m);                                                                   \
  _DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_ISTL(_m)

// end: this is what we need for the .so


#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BINDINGS_HH
