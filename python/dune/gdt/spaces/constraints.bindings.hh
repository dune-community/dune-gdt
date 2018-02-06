// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
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
#include <python/dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/la/container.hh>
#include <python/dune/xt/la/container.bindings.hh>

#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/spaces/cg/interface.hh>
#include <dune/gdt/spaces/cg.bindings.hh>

#include <dune/gdt/spaces/constraints.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class I, class G>
class DirichletConstraints
{
public:
  typedef GDT::DirichletConstraints<I> type;
  typedef pybind11::class_<type> bound_type;

  static bound_type bind(pybind11::module& m, const std::string& layer_name)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto grid_name = XT::Grid::bindings::grid_name<G>::value();
    const auto ClassName = XT::Common::to_camel_case("DirichletConstraints_" + layer_name + "_" + grid_name);

    bound_type c(m, ClassName.c_str(), ClassName.c_str());
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
  } // ... addbind(...)

private:
  template <class T, bool is_cg = is_cg_space<T>::value>
  struct addbind_assembler
  {
    template <class GL, class A>
    void operator()(pybind11::class_<GDT::SystemAssembler<T, GL, A>, XT::Grid::Walker<GL>>& bound_system_assembler)
    {
      using namespace pybind11::literals;

      bound_system_assembler.def("append",
                                 [](GDT::SystemAssembler<T, GL, A>& self,
                                    GDT::DirichletConstraints<XT::Grid::extract_intersection_t<GL>>& constraints) {
                                   self.append(constraints);
                                 },
                                 "dirichlet_constraints"_a);
    } // ... addbind(...)
  }; // struct addbind_assembler

  template <class T>
  struct addbind_assembler<T, false>
  {
    template <class GL, class A>
    void operator()(pybind11::class_<GDT::SystemAssembler<T, GL, A>, XT::Grid::Walker<GL>>& /*bound_system_assembler*/)
    {
    }
  };

public:
  template <class T, class GL, class A>
  static void addbind(pybind11::class_<GDT::SystemAssembler<T, GL, A>, XT::Grid::Walker<GL>>& bound_system_assembler)
  {
    addbind_assembler<T>()(bound_system_assembler);
  }
}; // class DirichletConstraints


} // namespace bindings
} // namespace GDT
} // namespace Dune


// begin: this is what we need for the .so

#define _DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(_c, _GRID, _layer, _backend, _la)                                      \
  Dune::GDT::bindings::                                                                                                \
      DirichletConstraints<Dune::XT::Grid::extract_intersection_t<                                                     \
                               typename Dune::XT::Grid::Layer<_GRID,                                                   \
                                                              Dune::XT::Grid::Layers::_layer,                          \
                                                              Dune::XT::Grid::Backends::_backend,                      \
                                                              Dune::XT::Grid::DD::SubdomainGrid<_GRID>>::type>,        \
                           _GRID>::addbind<Dune::XT::LA::Backends::_la>(_c)

#define _DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA_COMMON(_c, _GRID, _layer, _backend)
//  _DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(_c, _GRID, _layer, _backend, common_dense)

//#if HAVE_EIGEN
//#define _DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA_EIGEN(_c, _GRID, _layer, _backend)                                   \
//  _DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(_c, _GRID, _layer, _backend, eigen_sparse)
//#else
#define _DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA_EIGEN(_c, _GRID, _layer, _backend)
//#endif

#if HAVE_DUNE_ISTL
#define _DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA_ISTL(_c, _GRID, _layer, _backend)                                      \
  _DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(_c, _GRID, _layer, _backend, istl_sparse)
#else
#define _DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA_ISTL(_c, _GRID, _layer, _backend)
#endif

#define _DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA_ALL(_c, _GRID, _layer, _backend)                                       \
  _DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA_COMMON(_c, _GRID, _layer, _backend);                                         \
  _DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA_EIGEN(_c, _GRID, _layer, _backend);                                          \
  _DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA_ISTL(_c, _GRID, _layer, _backend)

#define _DUNE_GDT_SPACES_CONSTRAINTS_BIND(_m, _GRID, _layer, _backend, _layer_name)                                    \
  auto dirichlet_constraints_##_GRID##_##_layer##_##_backend = Dune::GDT::bindings::                                   \
      DirichletConstraints<Dune::XT::Grid::extract_intersection_t<                                                     \
                               typename Dune::XT::Grid::Layer<_GRID,                                                   \
                                                              Dune::XT::Grid::Layers::_layer,                          \
                                                              Dune::XT::Grid::Backends::_backend,                      \
                                                              Dune::XT::Grid::DD::SubdomainGrid<_GRID>>::type>,        \
                           _GRID>::bind(_m, _layer_name);                                                              \
  _DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA_ALL(                                                                         \
      dirichlet_constraints_##_GRID##_##_layer##_##_backend, _GRID, _layer, _backend)


//#if HAVE_ALBERTA
//#define _DUNE_GDT_SPACES_CONSTRAINTS_BIND_ALBERTA(_m)                                                                \
//  _DUNE_GDT_SPACES_CONSTRAINTS_BIND(_m, ALBERTA_2D, leaf, view, "");                                                 \
//  _DUNE_GDT_SPACES_CONSTRAINTS_BIND(_m, ALBERTA_2D, dd_subdomain, part, "dd_subdomain")
//#else
#define _DUNE_GDT_SPACES_CONSTRAINTS_BIND_ALBERTA(_m)
//#endif

#if HAVE_DUNE_ALUGRID
#define _DUNE_GDT_SPACES_CONSTRAINTS_BIND_ALU(_m)                                                                      \
  _DUNE_GDT_SPACES_CONSTRAINTS_BIND(_m, ALU_2D_SIMPLEX_CONFORMING, leaf, view, "leaf");                                \
  _DUNE_GDT_SPACES_CONSTRAINTS_BIND(_m, ALU_2D_SIMPLEX_CONFORMING, level, view, "level");                              \
  _DUNE_GDT_SPACES_CONSTRAINTS_BIND(_m, ALU_2D_SIMPLEX_CONFORMING, dd_subdomain, view, "dd_subdomain")
#else
#define _DUNE_GDT_SPACES_CONSTRAINTS_BIND_ALU(_m)
#endif

//#if HAVE_DUNE_UGGRID
//#define _DUNE_GDT_SPACES_CONSTRAINTS_BIND_UG(_m)                                                                     \
//  _DUNE_GDT_SPACES_CONSTRAINTS_BIND(_m, UG_2D, leaf, view, "leaf");                                                  \
//  _DUNE_GDT_SPACES_CONSTRAINTS_BIND(_m, UG_2D, level, view, "level");                                                \
//  _DUNE_GDT_SPACES_CONSTRAINTS_BIND(_m, UG_2D, dd_subdomain, part, "dd_subdomain")
//#else
#define _DUNE_GDT_SPACES_CONSTRAINTS_BIND_UG(_m)
//#endif

#define _DUNE_GDT_SPACES_CONSTRAINTS_BIND_YASP(_m)                                                                     \
  _DUNE_GDT_SPACES_CONSTRAINTS_BIND(_m, YASP_1D_EQUIDISTANT_OFFSET, leaf, view, "");                                   \
  _DUNE_GDT_SPACES_CONSTRAINTS_BIND(_m, YASP_1D_EQUIDISTANT_OFFSET, dd_subdomain, view, "dd_subdomain");               \
  _DUNE_GDT_SPACES_CONSTRAINTS_BIND(_m, YASP_2D_EQUIDISTANT_OFFSET, leaf, view, "");                                   \
  _DUNE_GDT_SPACES_CONSTRAINTS_BIND(_m, YASP_2D_EQUIDISTANT_OFFSET, dd_subdomain, view, "dd_subdomain")

#define DUNE_GDT_SPACES_CONSTRAINTS_BIND(_m)                                                                           \
  _DUNE_GDT_SPACES_CONSTRAINTS_BIND_ALBERTA(_m);                                                                       \
  _DUNE_GDT_SPACES_CONSTRAINTS_BIND_ALU(_m);                                                                           \
  _DUNE_GDT_SPACES_CONSTRAINTS_BIND_UG(_m);                                                                            \
  _DUNE_GDT_SPACES_CONSTRAINTS_BIND_YASP(_m)

// end: this is what we need for the .so


#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_SPACES_CONSTRAINTS_BINDINGS_HH
