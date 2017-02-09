// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_SPACES_CONSTRAINTS_BINDINGS_HH
#define DUNE_GDT_SPACES_CONSTRAINTS_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <boost/numeric/conversion/cast.hpp>

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/intersection.hh>
#include <dune/xt/grid/grids.bindings.hh>

#include <dune/xt/la/container.hh>
#include <dune/xt/la/container.bindings.hh>

#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/spaces/cg.bindings.hh>

#include "constraints.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <class I>
class DirichletConstraints
{
public:
  typedef GDT::DirichletConstraints<I> type;
  typedef pybind11::class_<type> bound_type;

  static bound_type bind(pybind11::module& m, const std::string& intersection_id)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    bound_type c(m, std::string("DirichletConstraints__" + intersection_id).c_str());
    c.def("__init__",
          [](type& self, const XT::Grid::BoundaryInfo<I>& boundary_info, const ssize_t size, const bool set) {
            try {
              new (&self) type(boundary_info, boost::numeric_cast<size_t>(size), set);
            } catch (boost::bad_numeric_cast& ee) {
              DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
                         "Given size has to be positive!\n\n The error in boost while converting '"
                             << size
                             << "' to '"
                             << XT::Common::Typename<size_t>::value()
                             << "' was: "
                             << ee.what());
            }
          },
          "boundary_info"_a,
          "size"_a,
          "set_diagonal_entries"_a = true,
          py::keep_alive<1, 2>());
    c.def("boundary_info", &type::boundary_info);
    c.def("size", &type::size);

    m.def(std::string("make_dirichlet_constraints").c_str(),
          [](const XT::Grid::BoundaryInfo<I>& boundary_info, const ssize_t size, const bool set) {
            size_t size__as_size_t = 0;
            try {
              size__as_size_t = boost::numeric_cast<size_t>(size);
            } catch (boost::bad_numeric_cast& ee) {
              DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
                         "Given size has to be positive!\n\n The error in boost while converting '"
                             << size
                             << "' to '"
                             << XT::Common::Typename<size_t>::value()
                             << "' was: "
                             << ee.what());
            }
            return type(boundary_info, size__as_size_t, set);
          },
          "boundary_info"_a,
          "size"_a,
          "set_diagonal_entries"_a,
          py::keep_alive<0, 1>());

    return c;
  } // ... bind(...)

  template <XT::LA::Backends la_backend, class R = double>
  static void addbind(bound_type& c)
  {
    typedef typename XT::LA::Container<R, la_backend>::MatrixType M;
    typedef typename XT::LA::Container<R, la_backend>::VectorType V;
    using namespace pybind11::literals;

    c.def("apply", [](const type& self, M& matrix) { self.apply(matrix); }, "matrix"_a);
    c.def("apply", [](const type& self, V& vector) { self.apply(vector); }, "vector"_a);
    c.def("apply", [](const type& self, M& matrix, V& vector) { self.apply(matrix, vector); }, "matrix"_a, "vector"_a);
  }

  template <class T, class GV, class A>
  static void addbind(pybind11::class_<GDT::SystemAssembler<T, GV, A>>& bound_system_assembler)
  {
    using namespace pybind11::literals;

    bound_system_assembler.def("append",
                               [](GDT::SystemAssembler<T, GV, A>& self,
                                  GDT::DirichletConstraints<typename XT::Grid::Intersection<GV>::type>& constraints) {
                                 self.append(constraints);
                               },
                               "dirichlet_constraints"_a);
  } // ... addbind_to_SystemAssembler(...)
}; // class addbind_to_SystemAssembler


#define DUNE_GDT_SPACES_CONSTRAINTS_BIND(_prefix, _GRID, _layer, _backend)                                             \
  _prefix class DirichletConstraints<typename XT::Grid::Intersection<                                                  \
      typename XT::Grid::Layer<_GRID, XT::Grid::Layers::_layer, XT::Grid::Backends::_backend>::type>::type>

#define DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(_prefix, _GRID, _layer, _backend, _la)                                  \
  _prefix void DirichletConstraints<typename XT::Grid::Intersection<                                                   \
      typename XT::Grid::Layer<_GRID, XT::Grid::Layers::_layer, XT::Grid::Backends::_backend>::type>::type>::          \
      addbind<XT::LA::Backends::_la>(                                                                                  \
          typename DirichletConstraints<typename XT::Grid::Intersection<                                               \
              typename XT::Grid::Layer<_GRID, XT::Grid::Layers::_layer, XT::Grid::Backends::_backend>::type>::type>::  \
              bound_type&)

#define DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_ASSEMBLER_FEM(_prefix, _GRID)                                              \
  DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_ASSEMBLER_(_prefix, CG_SPACE(_GRID, leaf, fem, 1, 1, 1));                        \
  DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_ASSEMBLER_(_prefix, CG_SPACE(_GRID, level, fem, 1, 1, 1))

#define DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_ASSEMBLER_PDELAB(_prefix, _GRID)                                           \
  DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_ASSEMBLER_(_prefix, CG_SPACE(_GRID, leaf, pdelab, 1, 1, 1));                     \
  DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_ASSEMBLER_(_prefix, CG_SPACE(_GRID, level, pdelab, 1, 1, 1))

#define DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_ASSEMBLER_(_prefix, _SPACE)                                                \
  _prefix void DirichletConstraints<typename XT::Grid::Intersection<typename _SPACE::GridViewType>::type>::            \
      addbind<_SPACE, typename _SPACE::GridViewType, _SPACE>(                                                          \
          pybind11::class_<GDT::SystemAssembler<_SPACE, typename _SPACE::GridViewType, _SPACE>>&)


// these lines have to match the corresponding ones in the .hh header file
DUNE_GDT_SPACES_CONSTRAINTS_BIND(extern template, YASP_2D_EQUIDISTANT_OFFSET, leaf, view);
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(extern template, YASP_2D_EQUIDISTANT_OFFSET, leaf, view, common_dense);
#if HAVE_EIGEN
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(extern template, YASP_2D_EQUIDISTANT_OFFSET, leaf, view, eigen_dense);
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(extern template, YASP_2D_EQUIDISTANT_OFFSET, leaf, view, eigen_sparse);
#endif
#if HAVE_DUNE_ISTL
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(extern template, YASP_2D_EQUIDISTANT_OFFSET, leaf, view, istl_sparse);
#endif
#if HAVE_DUNE_FEM
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_ASSEMBLER_FEM(extern template, YASP_2D_EQUIDISTANT_OFFSET);
#endif
#if HAVE_DUNE_PDELAB
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_ASSEMBLER_PDELAB(extern template, YASP_2D_EQUIDISTANT_OFFSET);
#endif

#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID
DUNE_GDT_SPACES_CONSTRAINTS_BIND(extern template, ALU_2D_SIMPLEX_CONFORMING, leaf, view);
DUNE_GDT_SPACES_CONSTRAINTS_BIND(extern template, ALU_2D_SIMPLEX_CONFORMING, level, view);
DUNE_GDT_SPACES_CONSTRAINTS_BIND(extern template, ALU_2D_SIMPLEX_CONFORMING, leaf, view);
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(extern template, ALU_2D_SIMPLEX_CONFORMING, leaf, view, common_dense);
#if HAVE_EIGEN
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(extern template, ALU_2D_SIMPLEX_CONFORMING, leaf, view, eigen_dense);
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(extern template, ALU_2D_SIMPLEX_CONFORMING, leaf, view, eigen_sparse);
#endif
#if HAVE_DUNE_ISTL
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(extern template, ALU_2D_SIMPLEX_CONFORMING, leaf, view, istl_sparse);
#endif
#if HAVE_DUNE_FEM
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_ASSEMBLER_FEM(extern template, ALU_2D_SIMPLEX_CONFORMING);
#endif
#if HAVE_DUNE_PDELAB
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_ASSEMBLER_PDELAB(extern template, ALU_2D_SIMPLEX_CONFORMING);
#endif
#endif


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_SPACES_CONSTRAINTS_BINDINGS_HH
