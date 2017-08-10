// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_OPERATORS_FV_REALIZABILITY_HH
#define DUNE_GDT_OPERATORS_FV_REALIZABILITY_HH

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/common/fvector.hh>

#include <dune/xt/grid/walker/functors.hh>

#if HAVE_QHULL
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacetList.h>
#endif // HAVE_QHULL

#if HAVE_LPSOLVE
namespace lpsolve {
#include <lpsolve/lp_lib.h>
}
#endif // HAVE_LPSOLVE

namespace Dune {
namespace GDT {


template <class EntityType>
class NonLimitingRealizabilityLimiter : public std::function<void(const EntityType&)>
{
  typedef std::function<void(const EntityType&)> BaseType;

public:
  NonLimitingRealizabilityLimiter()
    : BaseType([](const EntityType&) {})
  {
  }

  template <class BasisFunctionType, class QuadratureType, class RangeFieldType = double>
  NonLimitingRealizabilityLimiter(const BasisFunctionType& /*basis_functions*/,
                                  const QuadratureType& /*quadrature*/,
                                  const RangeFieldType /*epsilon*/ = 1e-14)
    : BaseType([](const EntityType&) {})
  {
  }

  template <class DiscreteFunctionType>
  void set_source(const DiscreteFunctionType*)
  {
  }

  template <class ReconstructedValuesType>
  void set_reconstructed_values(const ReconstructedValuesType*)
  {
  }
}; // class NonLimitingRealizabilityLimiter


#if HAVE_QHULL

template <class DiscreteFunctionType, class BasisFunctionType, size_t dimDomain, size_t dimRange>
class ConvexHullLocalRealizabilityLimiter
    : public XT::Grid::Functor::Codim0<typename DiscreteFunctionType::SpaceType::GridLayerType>
{
  typedef typename DiscreteFunctionType::SpaceType::GridLayerType GridLayerType;
  typedef typename DiscreteFunctionType::EntityType EntityType;
  typedef typename GridLayerType::template Codim<0>::Geometry::LocalCoordinate DomainType;
  typedef typename DiscreteFunctionType::RangeType RangeType;
  typedef typename GridLayerType::IndexSet IndexSetType;
  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef typename Dune::QuadratureRule<RangeFieldType, dimDomain> QuadratureType;
  typedef typename std::vector<std::pair<RangeType, RangeFieldType>> PlaneCoefficientsType;

public:
  // cell averages includes left and right boundary values as the two last indices in each dimension
  explicit ConvexHullLocalRealizabilityLimiter(const BasisFunctionType& basis_functions,
                                               const QuadratureType& quadrature,
                                               RangeFieldType epsilon = 1e-14)
    : basis_functions_(basis_functions)
    , quadrature_(quadrature)
    , epsilon_(epsilon)
  {
    if (is_instantiated_)
      DUNE_THROW(InvalidStateException,
                 "This class uses several static variables to save its state between time "
                 "steps, so using several instances at the same time may result in undefined "
                 "behavior!");
    is_instantiated_ = true;
    if (!plane_coefficients_)
      plane_coefficients_ = calculate_plane_coefficients();
  }

  ~ConvexHullLocalRealizabilityLimiter()
  {
    is_instantiated_ = false;
  }

  void set_source(const DiscreteFunctionType* source)
  {
    source_ = source;
    index_set_ = &(source_->space().grid_layer().indexSet());
  }

  void set_reconstructed_values(
      std::vector<std::map<DomainType, RangeType, XT::Common::FieldVectorLess>>* reconstructed_values)
  {
    reconstructed_values_ = reconstructed_values;
  }

  void apply_local(const EntityType& entity)
  {
    auto& local_reconstructed_values = (*reconstructed_values_)[index_set_->index(entity)];

    // get cell average
    const RangeType& u_bar =
        source_->local_function(entity)->evaluate(entity.geometry().local(entity.geometry().center()));

    // vector to store thetas for each local reconstructed value
    std::vector<RangeFieldType> thetas(local_reconstructed_values.size(), -epsilon_);

    size_t ll = -1;
    for (const auto& pair : local_reconstructed_values) {
      ++ll;
      // rescale u_l, u_bar
      auto u_l = pair.second;
      auto u_bar_minus_u_l = u_bar - u_l;
      const auto factor = basis_functions_.realizability_limiter_max(u_l, u_bar);
      u_l /= factor;
      u_bar_minus_u_l /= factor;

      for (const auto& coeffs : *plane_coefficients_) {
        const RangeType& a = coeffs.first;
        const RangeFieldType& b = coeffs.second;
        RangeFieldType theta_li = (b - a * u_l) / (a * u_bar_minus_u_l);
        if (XT::Common::FloatCmp::ge(theta_li, -epsilon_) && XT::Common::FloatCmp::le(theta_li, 1.))
          thetas[ll] = std::max(thetas[ll], theta_li);
      } // coeffs
    } // ll
    for (auto& theta : thetas)
      theta = std::min(epsilon_ + theta, 1.);

    auto theta_entity = *std::max_element(thetas.begin(), thetas.end());
    if (XT::Common::FloatCmp::ne(theta_entity, 0.)) {
      for (auto& pair : local_reconstructed_values) {
        auto& u = pair.second;
        auto u_scaled = u;
        u_scaled *= (1 - theta_entity);
        auto u_bar_scaled = u_bar;
        u_bar_scaled *= theta_entity;
        u = u_scaled + u_bar_scaled;
      }
    }
  } // void apply_local(...)

private:
  // calculate half space representation of realizable set
  std::shared_ptr<const PlaneCoefficientsType> calculate_plane_coefficients()
  {
    using orgQhull::Qhull;
    Qhull qhull;
    std::vector<FieldVector<RangeFieldType, dimRange>> points(quadrature_.size() + 1);
    points[0] = FieldVector<RangeFieldType, dimRange>(0);
    size_t ii = 1;
    for (const auto& quad_point : quadrature_)
      points[ii++] = basis_functions_.evaluate(quad_point.position());

    qhull.runQhull("Realizable set", int(dimRange), int(points.size()), &(points[0][0]), "Qt T1");
    //    qhull.outputQhull("n");
    const auto facet_end = qhull.endFacet();
    std::vector<std::pair<RangeType, RangeFieldType>> plane_coefficients(qhull.facetList().count());
    ii = 0;
    for (auto facet = qhull.beginFacet(); facet != facet_end; facet = facet.next(), ++ii) {
      for (size_t jj = 0; jj < dimRange; ++jj)
        plane_coefficients[ii].first[jj] = *(facet.hyperplane().coordinates() + jj);
      plane_coefficients[ii].second = -facet.hyperplane().offset();
    }
    return std::make_shared<const PlaneCoefficientsType>(plane_coefficients);
  }

  const DiscreteFunctionType* source_;
  const IndexSetType* index_set_;
  std::vector<std::map<DomainType, RangeType, XT::Common::FieldVectorLess>>* reconstructed_values_;
  const BasisFunctionType& basis_functions_;
  const QuadratureType& quadrature_;
  const RangeFieldType epsilon_;
  static bool is_instantiated_;
  static std::shared_ptr<const PlaneCoefficientsType> plane_coefficients_;
}; // class ConvexHullLocalRealizabilityLimiter

template <class DiscreteFunctionType, class BasisFunctionType, size_t dimDomain, size_t dimRange>
bool ConvexHullLocalRealizabilityLimiter<DiscreteFunctionType, BasisFunctionType, dimDomain, dimRange>::
    is_instantiated_ = false;

template <class DiscreteFunctionType, class BasisFunctionType, size_t dimDomain, size_t dimRange>
std::shared_ptr<const typename ConvexHullLocalRealizabilityLimiter<DiscreteFunctionType,
                                                                   BasisFunctionType,
                                                                   dimDomain,
                                                                   dimRange>::PlaneCoefficientsType>
    ConvexHullLocalRealizabilityLimiter<DiscreteFunctionType, BasisFunctionType, dimDomain, dimRange>::
        plane_coefficients_;

#else // HAVE_QHULL

template <class DiscreteFunctionType, class BasisFunctionType, size_t dimDomain, size_t dimRange>
class ConvexHullLocalRealizabilityLimiter
{
  static_assert(Dune::AlwaysFalse<DiscreteFunctionType>::value, "You are missing Qhull!");
};

#endif // HAVE_QHULL

#if HAVE_LPSOLVE

template <class DiscreteFunctionType, class BasisFunctionType, size_t dimDomain, size_t dimRange>
class LPLocalRealizabilityLimiter
    : public XT::Grid::Functor::Codim0<typename DiscreteFunctionType::SpaceType::GridLayerType>
{
  typedef typename DiscreteFunctionType::SpaceType::GridLayerType GridLayerType;
  typedef typename DiscreteFunctionType::EntityType EntityType;
  typedef typename GridLayerType::template Codim<0>::Geometry::LocalCoordinate DomainType;
  typedef typename DiscreteFunctionType::RangeType RangeType;
  typedef typename GridLayerType::IndexSet IndexSetType;
  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef typename Dune::QuadratureRule<RangeFieldType, dimDomain> QuadratureType;
  using BasisValuesMatrixType = std::vector<RangeType>;

public:
  // cell averages includes left and right boundary values as the two last indices in each dimension
  explicit LPLocalRealizabilityLimiter(const BasisFunctionType& basis_functions,
                                       const QuadratureType& quadrature,
                                       RangeFieldType epsilon = 1e-14)
    : basis_functions_(basis_functions)
    , quadrature_(quadrature)
    , num_quad_points_(quadrature_.size())
    , epsilon_(epsilon)
    , M_(quadrature_.size())
  {
    for (size_t ii = 0; ii < quadrature_.size(); ++ii)
      M_[ii] = basis_functions_.evaluate(quadrature_[ii].position());
  }

  void set_source(const DiscreteFunctionType* source)
  {
    source_ = source;
    index_set_ = &(source_->space().grid_layer().indexSet());
  }

  void set_reconstructed_values(
      std::vector<std::map<DomainType, RangeType, XT::Common::FieldVectorLess>>* reconstructed_values)
  {
    reconstructed_values_ = reconstructed_values;
  }

  void apply_local(const EntityType& entity)
  {
    auto& local_reconstructed_values = (*reconstructed_values_)[index_set_->index(entity)];

    // get cell average
    const RangeType& u_bar =
        source_->local_function(entity)->evaluate(entity.geometry().local(entity.geometry().center()));

    // vector to store thetas for each local reconstructed value
    std::vector<RangeFieldType> thetas(local_reconstructed_values.size(), -epsilon_);

    size_t ll = -1;
    for (const auto& pair : local_reconstructed_values) {
      ++ll;
      auto u_l = pair.second;
      auto u_l_minus_u_bar = u_l - u_bar;

      //            const auto max_val = *std::max_element(u_l_minus_u_bar.begin(), u_l_minus_u_bar.end());
      //            const auto min_val = *std::min_element(u_l_minus_u_bar.begin(), u_l_minus_u_bar.end());
      //            const auto abs_max = std::max(std::abs(max_val), std::abs(min_val));
      //      if (XT::Common::FloatCmp::le(abs_max, 1e-10)) {
      //        thetas[ll] = 1.;
      //      } else {
      // solve LP:
      // min \theta s.t.
      // (\sum x_i v_i) + \theta (u - \bar{u}) = u
      // x_i, \theta >= 0
      auto theta = solve_linear_program(u_l, u_l_minus_u_bar);
      assert(XT::Common::FloatCmp::ge(theta, epsilon_));
      if (theta < 1 + epsilon_)
        theta = theta - epsilon_;
      else
        theta = 1.;
      thetas[ll] = theta;
      //      }
    } // ll

    auto theta_entity = *std::min_element(thetas.begin(), thetas.end());
    if (XT::Common::FloatCmp::ne(theta_entity, 1.)) {
      for (auto& pair : local_reconstructed_values) {
        auto& u = pair.second;
        auto u_scaled = u;
        u_scaled *= theta_entity;
        auto u_bar_scaled = u_bar;
        u_bar_scaled *= 1 - theta_entity;
        u = u_scaled + u_bar_scaled;
      }
    }
  } // void apply_local(...)

private:
  RangeFieldType solve_linear_program(RangeType u_l, RangeType u_l_minus_u_bar)
  {
    typename lpsolve::lprec* lp;

    //     scale
    //    const auto max_val = *std::max_element(u_l_minus_u_bar.begin(), u_l_minus_u_bar.end());
    //    const auto min_val = *std::min_element(u_l_minus_u_bar.begin(), u_l_minus_u_bar.end());
    //    const auto abs_max = std::max(std::abs(max_val), std::abs(min_val));
    //    u_l_minus_u_bar /= abs_max;

    // We start with creating a model with dimRange+1 rows and num_quad_points+1 columns */
    constexpr int num_rows = int(dimRange + 1);
    int num_cols = int(quadrature_.size() + 1); /* variables are x_1, ..., x_{num_quad_points}, \theta */
    lp = lpsolve::make_lp(num_rows, num_cols);
    if (!lp)
      DUNE_THROW(Dune::MathError, "Couldn't construct linear program");

    /* let us name our variables. Not required, but can be useful for debugging */
    std::vector<char> name;
    for (int ii = 0; ii < num_cols; ++ii) {
      auto name_string = "x" + XT::Common::to_string(ii + 1);
      name.resize(name_string.size());
      std::copy(name_string.begin(), name_string.end(), name.begin());
      lpsolve::set_col_name(lp, ii + 1, name.data());
    }
    name.resize(5);
    name = {'t', 'h', 'e', 't', 'a'};
    lpsolve::set_col_name(lp, num_cols, name.data());

    // In the call to set_column, the first entry (row 0) is the value of the objective function
    // (c_i in the objective function c^T x), the other entries are the entries of the i-th column
    // in the constraints matrix. The 0-th column is the rhs vector. The entry (0, 0) corresponds
    // to the initial value of the objective function.
    std::array<REAL, num_rows + 1> column;
    // set rhs (column 0)
    column[0] = 0.;
    std::copy(u_l.begin(), u_l.end(), column.begin() + 1);
    column[dimRange + 1] = 1.;
    lpsolve::set_rh_vec(lp, column.data());
    lpsolve::set_rh(lp, 0, column[0]);
    // set columns for quadrature points
    column[0] = 0.;
    for (int ii = 0; ii < int(M_.size()); ++ii) {
      const auto& v_i = M_[ii];
      std::copy(v_i.begin(), v_i.end(), column.begin() + 1);
      column[dimRange + 1] = 1.;
      lpsolve::set_column(lp, ii + 1, column.data());
    }
    // set last column
    column[0] = 1.;
    std::copy(u_l_minus_u_bar.begin(), u_l_minus_u_bar.end(), column.begin() + 1);
    column[dimRange + 1] = 0.;
    std::cout << "theta col" << XT::Common::to_string(column, 15) << std::endl;
    lpsolve::set_column(lp, num_cols, column.data());
    for (int ii = 1; ii <= num_rows; ++ii)
      lpsolve::set_constr_type(lp, ii, EQ);

    // set bounds for all variables. This should not be necessary, as 0 <= x <= inf is
    // the default for all variables.
    for (int ii = 1; ii <= num_cols; ++ii)
      lpsolve::set_bounds(lp, ii, 0., lpsolve::get_infinite(lp));

    // set starting point for iteration. We can only set the variable to its lower or upper bound.
    // The variable is set to its lower bound if the value in initial_basis is negative and to its upper_
    // bound if the value is positive. We want to set all variables to its lower bound.
    //    std::array<REAL, num_rows + 1> initial_basis;
    //    const auto basisret = set_basis(lp, initial_basis, FALSE);
    //    if (basisret == FALSE)
    //      DUNE_THROW(Dune::MathError, "Failed to set initial basis!");


    /* make LP maximize instead of minimize */
    //    set_maxim(lp);

    /* print LP */
    /* this only works if this is a console application. If not, use write_lp and a filename */
    lpsolve::write_LP(lp, stdout);

    /* I only want to see important messages on screen while solving */
    lpsolve::set_verbose(lp, FULL);

    /* Now let lpsolve calculate a solution */
    const auto solve_status = lpsolve::solve(lp);
    if (solve_status != OPTIMAL) {
      std::cout << solve_status << std::endl;
      DUNE_THROW(Dune::MathError, "An unexpected error occured while solving the linear program");
    }
    RangeFieldType lambda = lpsolve::get_objective(lp); // * abs_max;

    /* create space large enough for one row */
    REAL* row = (REAL*)malloc(num_cols * sizeof(REAL));
    /* variable values */
    lpsolve::get_variables(lp, row);
    for (size_t ii = 0; ii <= size_t(num_cols); ++ii)
      std::cout << row[ii] << " ";
    //    assert(XT::Common::FloatCmp::eq(lambda, row[num_cols - 1]));

    /* free allocated memory */
    if (row)
      free(row);
    if (lp)
      lpsolve::delete_lp(lp);

    return lambda;
  }

  const DiscreteFunctionType* source_;
  const IndexSetType* index_set_;
  std::vector<std::map<DomainType, RangeType, XT::Common::FieldVectorLess>>* reconstructed_values_;
  const BasisFunctionType& basis_functions_;
  const QuadratureType& quadrature_;
  const size_t num_quad_points_;
  const RangeFieldType epsilon_;
  BasisValuesMatrixType M_;
}; // class LPLocalRealizabilityLimiter

#else // HAVE_LPSOLVE

template <class DiscreteFunctionType, class BasisFunctionType, size_t dimDomain, size_t dimRange>
class LPLocalRealizabilityLimiter
{
  static_assert(Dune::AlwaysFalse<DiscreteFunctionType>::value, "You are missing LPSolve!");
};

#endif // HAVE_LPSOLVE


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_REALIZABILITY_HH
