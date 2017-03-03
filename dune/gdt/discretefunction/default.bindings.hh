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


template <class S, class V>
class ConstDiscreteFunction
{
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

  static bound_type bind(pybind11::module& m, const std::string& space_id, const std::string& la_id)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    bound_type c(m, std::string("ConstDiscreteFunction__" + space_id + "__" + la_id).c_str());
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


template <class S, class V>
class DiscreteFunction
{
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

    bound_type c(m, std::string("DiscreteFunction__" + space_id + "__" + la_id).c_str());
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


#define DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT(_prefix, _GRID, _LA)                                                \
  _prefix class ConstDiscreteFunction<FV_SPACE(_GRID, leaf, gdt, 1, 1), _LA>;                                          \
  _prefix class ConstDiscreteFunction<FV_SPACE(_GRID, level, gdt, 1, 1), _LA>;                                         \
  _prefix class DiscreteFunction<FV_SPACE(_GRID, leaf, gdt, 1, 1), _LA>;                                               \
  _prefix class DiscreteFunction<FV_SPACE(_GRID, level, gdt, 1, 1), _LA>

#define DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM(_prefix, _GRID, _LA)                                                \
  _prefix class ConstDiscreteFunction<CG_SPACE(_GRID, leaf, fem, 1, 1, 1), _LA>;                                       \
  _prefix class ConstDiscreteFunction<CG_SPACE(_GRID, level, fem, 1, 1, 1), _LA>;                                      \
  _prefix class ConstDiscreteFunction<DG_SPACE(_GRID, leaf, fem, 1, 1, 1), _LA>;                                       \
  _prefix class ConstDiscreteFunction<DG_SPACE(_GRID, level, fem, 1, 1, 1), _LA>;                                      \
  _prefix class ConstDiscreteFunction<DG_SPACE(_GRID, leaf, fem, 2, 1, 1), _LA>;                                       \
  _prefix class ConstDiscreteFunction<DG_SPACE(_GRID, level, fem, 2, 1, 1), _LA>;                                      \
  _prefix class DiscreteFunction<CG_SPACE(_GRID, leaf, fem, 1, 1, 1), _LA>;                                            \
  _prefix class DiscreteFunction<CG_SPACE(_GRID, level, fem, 1, 1, 1), _LA>;                                           \
  _prefix class DiscreteFunction<DG_SPACE(_GRID, leaf, fem, 1, 1, 1), _LA>;                                            \
  _prefix class DiscreteFunction<DG_SPACE(_GRID, level, fem, 1, 1, 1), _LA>;                                           \
  _prefix class DiscreteFunction<DG_SPACE(_GRID, leaf, fem, 2, 1, 1), _LA>;                                            \
  _prefix class DiscreteFunction<DG_SPACE(_GRID, level, fem, 2, 1, 1), _LA>


#define DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB(_prefix, _GRID, _LA)                                             \
  _prefix class ConstDiscreteFunction<CG_SPACE(_GRID, leaf, pdelab, 1, 1, 1), _LA>;                                    \
  _prefix class ConstDiscreteFunction<CG_SPACE(_GRID, level, pdelab, 1, 1, 1), _LA>;                                   \
  _prefix class DiscreteFunction<CG_SPACE(_GRID, leaf, pdelab, 1, 1, 1), _LA>;                                         \
  _prefix class DiscreteFunction<CG_SPACE(_GRID, level, pdelab, 1, 1, 1), _LA>


// these lines have to match the corresponding ones in the .cc source file
DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT(extern template, YASP_2D_EQUIDISTANT_OFFSET, COMMON_DENSE_VECTOR);
#if HAVE_EIGEN
DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT(extern template, YASP_2D_EQUIDISTANT_OFFSET, EIGEN_DENSE_VECTOR);
#endif
#if HAVE_DUNE_ISTL
DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT(extern template, YASP_2D_EQUIDISTANT_OFFSET, ISTL_DENSE_VECTOR);
#endif
#if HAVE_DUNE_FEM
DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM(extern template, YASP_2D_EQUIDISTANT_OFFSET, COMMON_DENSE_VECTOR);
#if HAVE_EIGEN
DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM(extern template, YASP_2D_EQUIDISTANT_OFFSET, EIGEN_DENSE_VECTOR);
#endif
#if HAVE_DUNE_ISTL
DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM(extern template, YASP_2D_EQUIDISTANT_OFFSET, ISTL_DENSE_VECTOR);
#endif
#endif // HAVE_DUNE_FEM
#if HAVE_DUNE_PDELAB
DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB(extern template, YASP_2D_EQUIDISTANT_OFFSET, COMMON_DENSE_VECTOR);
#if HAVE_EIGEN
DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB(extern template, YASP_2D_EQUIDISTANT_OFFSET, EIGEN_DENSE_VECTOR);
#endif
#if HAVE_DUNE_ISTL
DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB(extern template, YASP_2D_EQUIDISTANT_OFFSET, ISTL_DENSE_VECTOR);
#endif
#endif // HAVE_DUNE_PDELAB

#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID
DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT(extern template, ALU_2D_SIMPLEX_CONFORMING, COMMON_DENSE_VECTOR);
#if HAVE_EIGEN
DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT(extern template, ALU_2D_SIMPLEX_CONFORMING, EIGEN_DENSE_VECTOR);
#endif
#if HAVE_DUNE_ISTL
DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_GDT(extern template, ALU_2D_SIMPLEX_CONFORMING, ISTL_DENSE_VECTOR);
#endif
#if HAVE_DUNE_FEM
DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM(extern template, ALU_2D_SIMPLEX_CONFORMING, COMMON_DENSE_VECTOR);
#if HAVE_EIGEN
DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM(extern template, ALU_2D_SIMPLEX_CONFORMING, EIGEN_DENSE_VECTOR);
#endif
#if HAVE_DUNE_ISTL
DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_FEM(extern template, ALU_2D_SIMPLEX_CONFORMING, ISTL_DENSE_VECTOR);
#endif
#endif // HAVE_DUNE_FEM
#if HAVE_DUNE_PDELAB
DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB(extern template, ALU_2D_SIMPLEX_CONFORMING, COMMON_DENSE_VECTOR);
#if HAVE_EIGEN
DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB(extern template, ALU_2D_SIMPLEX_CONFORMING, EIGEN_DENSE_VECTOR);
#endif
#if HAVE_DUNE_ISTL
DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BIND_PDELAB(extern template, ALU_2D_SIMPLEX_CONFORMING, ISTL_DENSE_VECTOR);
#endif
#endif // HAVE_DUNE_PDELAB
#endif // HAVE_ALUGRID || HAVE_DUNE_ALUGRID


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_DISCRETEFUNCTION_DEFAULT_BINDINGS_HH
