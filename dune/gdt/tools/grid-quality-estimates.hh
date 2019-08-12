// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_TOOLS_GRID_QUALITY_ESTIMATES_HH
#define DUNE_GDT_TOOLS_GRID_QUALITY_ESTIMATES_HH

#include <limits>
#include <vector>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/conversion.hh>
#include <dune/xt/la/generalized-eigen-solver.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/intersection.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/integrands/laplace.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {


template <class GV, size_t r>
double estimate_inverse_inequality_constant(const SpaceInterface<GV, r>& space)
{
  DUNE_THROW_IF(!XT::Common::Lapacke::available(), XT::Common::Exceptions::dependency_missing, "lapacke");
  using E = XT::Grid::extract_entity_t<GV>;
  double result = std::numeric_limits<double>::min();
  auto basis = space.basis().localize();
  for (auto&& element : elements(space.grid_view())) {
    basis->bind(element);
    const double h = XT::Grid::diameter(element);
    auto H1_product_matrix = XT::LA::convert_to<XT::LA::CommonDenseMatrix<double>>(
        LocalElementIntegralBilinearForm<E, r>(LocalLaplaceIntegrand<E, r>()).apply2(*basis, *basis));
    auto L2_product_matrix = XT::LA::convert_to<XT::LA::CommonDenseMatrix<double>>(
        LocalElementIntegralBilinearForm<E, r>(LocalElementProductIntegrand<E, r>()).apply2(*basis, *basis));
    auto evs =
        XT::LA::make_generalized_eigen_solver(H1_product_matrix,
                                              L2_product_matrix,
                                              {{"type", XT::LA::generalized_eigen_solver_types(H1_product_matrix)[0]},
                                               {"compute_eigenvectors", "false"},
                                               {"assert_real_eigenvalues", "1e-15"}})
            .real_eigenvalues();
    double min_ev = std::numeric_limits<double>::max();
    for (auto&& ev : evs)
      if (std::abs(ev) > 1e-7) // TODO: find a better tolerance here!
        min_ev = std::min(min_ev, ev);
    // the smalles nonzero eigenvalue is (C_I / h)^2
    result = std::max(result, h * std::sqrt(min_ev));
  }
  return result;
} // ... estimate_inverse_inequality_constant(...)


template <class GV, size_t r>
double estimate_combined_inverse_trace_inequality_constant(const SpaceInterface<GV, r>& space)
{
  DUNE_THROW_IF(!XT::Common::Lapacke::available(), XT::Common::Exceptions::dependency_missing, "lapacke");
  using E = XT::Grid::extract_entity_t<GV>;
  using I = XT::Grid::extract_intersection_t<GV>;
  double result = std::numeric_limits<double>::min();
  auto basis = space.basis().localize();
  for (auto&& element : elements(space.grid_view())) {
    basis->bind(element);
    const double h = XT::Grid::diameter(element);
    XT::LA::CommonDenseMatrix<double> L2_face_product_matrix(basis->size(), basis->size(), 0.);
    DynamicMatrix<double> tmp_L2_face_product_matrix(basis->size(), basis->size(), 0.);
    for (auto&& intersection : intersections(space.grid_view(), element)) {
      LocalIntersectionIntegralBilinearForm<I, r>(LocalIntersectionProductIntegrand<I, r>(1.))
          .apply2(intersection, *basis, *basis, tmp_L2_face_product_matrix);
      for (size_t ii = 0; ii < basis->size(); ++ii)
        for (size_t jj = 0; jj < basis->size(); ++jj)
          L2_face_product_matrix.add_to_entry(ii, jj, tmp_L2_face_product_matrix[ii][jj]);
    }
    auto L2_element_product_matrix = XT::LA::convert_to<XT::LA::CommonDenseMatrix<double>>(
        LocalElementIntegralBilinearForm<E, r>(LocalElementProductIntegrand<E, r>(1.)).apply2(*basis, *basis));
    auto evs = XT::LA::make_generalized_eigen_solver(
                   L2_face_product_matrix,
                   L2_element_product_matrix,
                   {{"type", XT::LA::generalized_eigen_solver_types(L2_face_product_matrix)[0]},
                    {"compute_eigenvectors", "false"},
                    {"assert_real_eigenvalues", "1e-15"}})
                   .real_eigenvalues();
    double min_ev = std::numeric_limits<double>::max();
    for (auto&& ev : evs)
      if (std::abs(ev) > 1e-7) // TODO: find a better tolerance here!
        min_ev = std::min(min_ev, ev);
    // the smalles nonzero eigenvalue is (C_M (1 + C_I)) / h
    result = std::max(result, h * min_ev);
  }
  return result;
} // ... estimate_combined_inverse_trace_inequality_constant(...)


template <class GV>
double estimate_element_to_intersection_equivalence_constant(
    const GridView<GV>& grid_view,
    const std::function<double(const XT::Grid::extract_intersection_t<GridView<GV>>&)>& intersection_diameter =
        [](const auto& intersection) {
          if (GridView<GV>::dimension == 1) {
            if (intersection.neighbor())
              return 0.5 * (XT::Grid::diameter(intersection.inside()) + XT::Grid::diameter(intersection.outside()));
            else
              return XT::Grid::diameter(intersection.inside());
          } else
            return XT::Grid::diameter(intersection);
        })
{
  auto result = std::numeric_limits<double>::min();
  for (auto&& element : elements(grid_view)) {
    const double h = XT::Grid::diameter(element);
    for (auto&& intersection : intersections(grid_view, element))
      result = std::max(result, intersection_diameter(intersection) / h);
    return result;
  }
} // ... estimate_element_to_intersection_equivalence_constant(...)


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TOOLS_GRID_QUALITY_ESTIMATES_HH
