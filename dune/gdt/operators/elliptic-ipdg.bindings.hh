// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BINDINGS_HH
#define DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/grids.bindings.hh>

#include <dune/gdt/spaces.bindings.hh>

#include "elliptic-ipdg.hh"
#include "base.bindings.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <class DF,
          typename DT, // may be void
          class R,
          LocalEllipticIpdgIntegrands::Method method,
          class M = typename XT::LA::Container<typename R::RangeFieldType>::MatrixType /*,
          class GV = typename R::GridViewType,
          class S = R,
          class F = typename R::RangeFieldType*/>
class EllipticIpdgMatrixOperator
{
public:
  typedef GDT::EllipticIpdgMatrixOperator<DF, DT, R, method, M /*, GV, S, F*/> type;
  typedef pybind11::class_<type> bound_type;

private:
  template <bool single_diffusion = std::is_same<DT, void>::value, bool anything = false>
  struct diffusion_switch
  {
    static std::string suffix()
    {
      return "diffusion_factor_and_tensor";
    }

    template <class C>
    static void addbind_factory_methods(pybind11::module& m, const std::string& method_id, const std::string& la_id)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      m.def(
          std::string(method_id + "__" + la_id).c_str(),
          [](const DF& diffusion_factor,
             const DT& diffusion_tensor,
             const XT::Grid::BoundaryInfo<typename R::GridViewType::Intersection>& boundary_info,
             const R& space,
             const size_t over_integrate) {
            return make_elliptic_ipdg_matrix_operator<M, method>(
                       diffusion_factor, diffusion_tensor, boundary_info, space, over_integrate)
                .release(); //         <- b.c. EllipticIpdgMatrixOperator is not movable, returning the raw pointer lets
          }, //                                                                     pybind11 correctly manage the memory
          "diffusion_factor"_a,
          "diffusion_tensor"_a,
          "boundary_info"_a,
          "space"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>(),
          py::keep_alive<0, 4>());

      m.def(
          std::string(method_id).c_str(),
          [](const DF& diffusion_factor,
             const DT& diffusion_tensor,
             const XT::Grid::BoundaryInfo<typename R::GridViewType::Intersection>& boundary_info,
             M& matrix,
             const R& space,
             const size_t over_integrate) {
            return make_elliptic_ipdg_matrix_operator<method>(
                       diffusion_factor, diffusion_tensor, boundary_info, matrix, space, over_integrate)
                .release(); //                                                                     <- s.a. for release()
          },
          "diffusion_factor"_a,
          "diffusion_tensor"_a,
          "boundary_info"_a,
          "matrix"_a,
          "space"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>(),
          py::keep_alive<0, 4>(),
          py::keep_alive<0, 5>());
    } // ... addbind_factory_methods(...)

  }; // struct diffusion_switch

  template <bool anything>
  struct diffusion_switch<true, anything>
  {
    static std::string suffix()
    {
      return "single_diffusion";
    }

    template <class C>
    static void addbind_factory_methods(pybind11::module& m, const std::string& method_id, const std::string& la_id)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      m.def(
          std::string(method_id + "__" + la_id).c_str(),
          [](const DF& diffusion,
             const XT::Grid::BoundaryInfo<typename R::GridViewType::Intersection>& boundary_info,
             const R& space,
             const size_t over_integrate) {
            return make_elliptic_ipdg_matrix_operator<M, method>(diffusion, boundary_info, space, over_integrate)
                .release(); //                                                                     <- s.a. for release()
          },
          "diffusion"_a,
          "boundary_info"_a,
          "space"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>());

      m.def(
          std::string(method_id).c_str(),
          [](const DF& diffusion,
             const XT::Grid::BoundaryInfo<typename R::GridViewType::Intersection>& boundary_info,
             M& matrix,
             const R& space,
             const size_t over_integrate) {
            return make_elliptic_ipdg_matrix_operator<method>(diffusion, boundary_info, matrix, space, over_integrate)
                .release(); //                                                                     <- s.a. for release()
          },
          "diffusion"_a,
          "boundary_info"_a,
          "matrix"_a,
          "space"_a,
          "over_integrate"_a = 0,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>(),
          py::keep_alive<0, 4>());
    } // ... addbind_factory_methods(...)
  }; // struct diffusion_switch<..., void>

public:
  static bound_type
  bind(pybind11::module& m, const std::string& space_id, const std::string& la_id, const std::string& method_id)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const std::string suffix = la_id + "__" + space_id + "_" + diffusion_switch<>::suffix();

    auto c = MatrixOperatorBase<type>::bind(m, "Elliptic" + method_id + "MatrixOperator__" + suffix);

    diffusion_switch<>::template addbind_factory_methods<type>(
        m, "make_elliptic_" + XT::Common::to_lower(method_id) + "_matrix_operator", la_id);

    return c;
  } // ... bind(...)

}; // EllipticIpdgMatrixOperator


// If everyone just had enough memory, we could just use a single line
//     DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND(template);
// in a source file to instantiate everything (together with
//     DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND(extern template);
// in this header. Alas, we can use the latter in this header, but need to distribute the load over several sources by
// using the specialized macros below...
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic" // because of the extra ; in some places
#endif

#define DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND(_PRE)                                                                    \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_YASPGRID(_PRE);                                                                \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_ALUGRID(_PRE)

#define DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_YASPGRID(_PRE)                                                           \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_LA(_PRE, YASP_2D_EQUIDISTANT_OFFSET)

#if HAVE_DUNE_ALUGRID
#define DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_ALUGRID(_PRE)                                                            \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_LA(_PRE, ALU_2D_SIMPLEX_CONFORMING)
#else
#define DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_ALUGRID(_PRE)
#endif // HAVE_DUNE_ALUGRID

#define DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_LA(_PRE, _GRID)                                                          \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_COMMON(_PRE, _GRID);                                                           \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_EIGEN(_PRE, _GRID);                                                            \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_ISTL(_PRE, _GRID)

#define DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_COMMON(_PRE, _GRID)                                                      \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_FEM(_PRE, _GRID, common_dense);                                                \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_GDT(_PRE, _GRID, common_dense);                                                \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_PDELAB(_PRE, _GRID, common_dense)

#if HAVE_EIGEN
#define DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_EIGEN(_PRE, _GRID)                                                       \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_FEM(_PRE, _GRID, eigen_dense);                                                 \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_FEM(_PRE, _GRID, eigen_sparse);                                                \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_GDT(_PRE, _GRID, eigen_dense);                                                 \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_GDT(_PRE, _GRID, eigen_sparse);                                                \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_PDELAB(_PRE, _GRID, eigen_dense);                                              \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_PDELAB(_PRE, _GRID, eigen_sparse)
#else
#define DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_EIGEN(_PRE, _GRID)
#endif

#if HAVE_DUNE_ISTL
#define DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_ISTL(_PRE, _GRID)                                                        \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_FEM(_PRE, _GRID, istl_dense);                                                  \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_GDT(_PRE, _GRID, istl_dense);                                                  \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_PDELAB(_PRE, _GRID, istl_dense)
#else
#define DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_ISTL(_PRE, _GRID)
#endif

#if HAVE_DUNE_FEM
#define DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_FEM(_PRE, _GRID, _LA)                                                    \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_METHODS(_PRE, CG_SPACE(_GRID, leaf, fem, 1, 1, 1), _LA);                       \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_METHODS(_PRE, CG_SPACE(_GRID, level, fem, 1, 1, 1), _LA);                      \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_METHODS(_PRE, DG_SPACE(_GRID, leaf, fem, 1, 1, 1), _LA);                       \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_METHODS(_PRE, DG_SPACE(_GRID, level, fem, 1, 1, 1), _LA);                      \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_METHODS(_PRE, DG_SPACE(_GRID, leaf, fem, 2, 1, 1), _LA);                       \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_METHODS(_PRE, DG_SPACE(_GRID, level, fem, 2, 1, 1), _LA)
#else // HAVE_DUNE_FEM
#define DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_FEM(_PRE, _GRID)
#endif

#define DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_GDT(_PRE, _GRID, _LA)                                                    \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_METHODS(_PRE, FV_SPACE(_GRID, leaf, gdt, 1, 1), _LA);                          \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_METHODS(_PRE, FV_SPACE(_GRID, level, gdt, 1, 1), _LA)

#if HAVE_DUNE_PDELAB
#define DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_PDELAB(_PRE, _GRID, _LA)                                                 \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_METHODS(_PRE, CG_SPACE(_GRID, leaf, pdelab, 1, 1, 1), _LA);                    \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_METHODS(_PRE, CG_SPACE(_GRID, level, pdelab, 1, 1, 1), _LA)
#else
#define DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_PDELAB(_PRE, _GRID, _LA)
#endif

#define DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_METHODS(_PRE, _SPACE, _LA)                                               \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_(_PRE, _SPACE, sipdg, _LA);                                                    \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_(_PRE, _SPACE, swipdg, _LA);                                                   \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_(_PRE, _SPACE, swipdg_affine_factor, _LA);                                     \
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_(_PRE, _SPACE, swipdg_affine_tensor, _LA)

#define DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_(_PRE, _SPACE, _METHOD, _LA)                                             \
  _PRE class EllipticIpdgMatrixOperator<XT::Functions::LocalizableFunctionInterface<typename _SPACE::EntityType,       \
                                                                                    typename _SPACE::DomainFieldType,  \
                                                                                    _SPACE::dimDomain,                 \
                                                                                    typename _SPACE::RangeFieldType,   \
                                                                                    1,                                 \
                                                                                    1>,                                \
                                        XT::Functions::LocalizableFunctionInterface<typename _SPACE::EntityType,       \
                                                                                    typename _SPACE::DomainFieldType,  \
                                                                                    _SPACE::dimDomain,                 \
                                                                                    typename _SPACE::RangeFieldType,   \
                                                                                    _SPACE::dimDomain,                 \
                                                                                    _SPACE::dimDomain>,                \
                                        _SPACE,                                                                        \
                                        LocalEllipticIpdgIntegrands::Method::_METHOD,                                  \
                                        typename XT::LA::Container<typename _SPACE::RangeFieldType,                    \
                                                                   XT::LA::Backends::_LA>::MatrixType>;                \
  _PRE class EllipticIpdgMatrixOperator<XT::Functions::LocalizableFunctionInterface<typename _SPACE::EntityType,       \
                                                                                    typename _SPACE::DomainFieldType,  \
                                                                                    _SPACE::dimDomain,                 \
                                                                                    typename _SPACE::RangeFieldType,   \
                                                                                    1,                                 \
                                                                                    1>,                                \
                                        void,                                                                          \
                                        _SPACE,                                                                        \
                                        LocalEllipticIpdgIntegrands::Method::_METHOD,                                  \
                                        typename XT::LA::Container<typename _SPACE::RangeFieldType,                    \
                                                                   XT::LA::Backends::_LA>::MatrixType>;                \
  _PRE class EllipticIpdgMatrixOperator<XT::Functions::LocalizableFunctionInterface<typename _SPACE::EntityType,       \
                                                                                    typename _SPACE::DomainFieldType,  \
                                                                                    _SPACE::dimDomain,                 \
                                                                                    typename _SPACE::RangeFieldType,   \
                                                                                    _SPACE::dimDomain,                 \
                                                                                    _SPACE::dimDomain>,                \
                                        void,                                                                          \
                                        _SPACE,                                                                        \
                                        LocalEllipticIpdgIntegrands::Method::_METHOD,                                  \
                                        typename XT::LA::Container<typename _SPACE::RangeFieldType,                    \
                                                                   XT::LA::Backends::_LA>::MatrixType>


DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND(extern template);


#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif


} // naemspace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BINDINGS_HH
