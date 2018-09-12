// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_HATFUNCTIONS_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_HATFUNCTIONS_HH

#include <vector>
#include <string>

#include <dune/xt/common/fmatrix.hh>

#include <dune/gdt/test/hyperbolic/problems/momentmodels/triangulation.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


template <class DomainFieldType,
          size_t dimDomain,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols = 1,
          size_t dimFlux = dimDomain>
class HatFunctionMomentBasis
{
  //  static_assert(false, "Not implemented for this dimension!");
};

template <class DomainFieldType, class RangeFieldType, size_t rangeDim, size_t rangeDimCols, size_t fluxDim>
class HatFunctionMomentBasis<DomainFieldType, 1, RangeFieldType, rangeDim, rangeDimCols, fluxDim>
    : public BasisfunctionsInterface<DomainFieldType, 1, RangeFieldType, rangeDim, rangeDimCols, fluxDim>
{
public:
  static const size_t dimDomain = 1;
  static const size_t dimRange = rangeDim;
  static const size_t dimRangeCols = rangeDimCols;

private:
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols, fluxDim> BaseType;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::MatrixType;
  using typename BaseType::StringifierType;
  using typename BaseType::QuadraturesType;
  using TriangulationType = typename BaseType::Triangulation1dType;
  template <class DiscreteFunctionType>
  using VisualizerType = typename BaseType::template VisualizerType<DiscreteFunctionType>;

  static std::string static_id()
  {
    return "hatfunctions";
  }

  HatFunctionMomentBasis(const QuadraturesType& quadratures)
    : BaseType(quadratures)
    , triangulation_(BaseType::create_1d_triangulation(dimRange - 1))
  {
  }

  HatFunctionMomentBasis(const size_t quad_order = 15, const size_t DXTC_DEBUG_ONLY(quad_refinements) = 0)
    : BaseType(BaseType::gauss_lobatto_quadratures(dimRange - 1, quad_order))
    , triangulation_(BaseType::create_1d_triangulation(dimRange - 1))
  {
    assert(quad_refinements == 0 && "Refinement of the quadrature intervals not implemented for this basis!");
  }

  virtual RangeType evaluate(const DomainType& v) const override final
  {
    RangeType ret(0);
    for (size_t ii = 0; ii < dimRange; ++ii) {
      if (ii < dimRange - 1 && XT::Common::FloatCmp::ge(v[0], triangulation_[ii])
          && XT::Common::FloatCmp::le(v[0], triangulation_[ii + 1]))
        ret[ii] = (v - triangulation_[ii + 1]) / (triangulation_[ii] - triangulation_[ii + 1]);
      if (ii > 0 && XT::Common::FloatCmp::ge(v[0], triangulation_[ii - 1])
          && XT::Common::FloatCmp::le(v[0], triangulation_[ii]))
        ret[ii] = (v - triangulation_[ii - 1]) / (triangulation_[ii] - triangulation_[ii - 1]);
    }
    return ret;
  } // ... evaluate(...)

  virtual RangeType evaluate(const DomainType& v, const size_t interval_index) const override final
  {
    RangeType ret(0);
    ret[interval_index] = (v - triangulation_[interval_index + 1])
                          / (triangulation_[interval_index] - triangulation_[interval_index + 1]);
    ret[interval_index + 1] =
        (v - triangulation_[interval_index]) / (triangulation_[interval_index + 1] - triangulation_[interval_index]);
    return ret;
  } // ... evaluate(...)

  virtual RangeType integrated() const override final
  {
    RangeType ret(0);
    ret[0] = triangulation_[1] - triangulation_[0];
    for (size_t ii = 1; ii < dimRange - 1; ++ii)
      ret[ii] = triangulation_[ii + 1] - triangulation_[ii - 1];
    ret[dimRange - 1] = triangulation_[dimRange - 1] - triangulation_[dimRange - 2];
    ret *= 0.5;
    return ret;
  }

  // returns matrix with entries <h_i h_j>
  virtual MatrixType mass_matrix() const override final
  {
    MatrixType ret(dimRange, dimRange, 0);
    ret[0][0] = (triangulation_[1] - triangulation_[0]) / 3.;
    for (size_t rr = 0; rr < dimRange; ++rr) {
      if (rr > 0 && rr < dimRange - 1)
        ret[rr][rr] = (triangulation_[rr + 1] - triangulation_[rr - 1]) / 3.;
      if (rr > 0)
        ret[rr][rr - 1] = (triangulation_[rr] - triangulation_[rr - 1]) / 6.;
      if (rr < dimRange - 1)
        ret[rr][rr + 1] = (triangulation_[rr + 1] - triangulation_[rr]) / 6.;
    }
    ret[dimRange - 1][dimRange - 1] = (triangulation_[dimRange - 1] - triangulation_[dimRange - 2]) / 3.;
    return ret;
  }

  virtual MatrixType mass_matrix_inverse() const override final
  {
    return tridiagonal_matrix_inverse<RangeFieldType, dimRange>(mass_matrix());
  }

  // returns matrix with entries <v h_i h_j>
  virtual FieldVector<MatrixType, 1> mass_matrix_with_v() const override final
  {
    MatrixType ret(dimRange, dimRange, 0.);
    ret[0][0] = (triangulation_[1] * triangulation_[1] + 2 * triangulation_[1] * triangulation_[0]
                 - 3 * triangulation_[0] * triangulation_[0])
                / 12.;
    for (size_t rr = 0; rr < dimRange; ++rr) {
      if (rr > 0 && rr < dimRange - 1)
        ret[rr][rr] = (triangulation_[rr + 1] * triangulation_[rr + 1] + 2 * triangulation_[rr + 1] * triangulation_[rr]
                       - 2 * triangulation_[rr] * triangulation_[rr - 1]
                       - triangulation_[rr - 1] * triangulation_[rr - 1])
                      / 12.;
      if (rr > 0)
        ret[rr][rr - 1] =
            (triangulation_[rr] * triangulation_[rr] - triangulation_[rr - 1] * triangulation_[rr - 1]) / 12.;
      if (rr < dimRange - 1)
        ret[rr][rr + 1] =
            (triangulation_[rr + 1] * triangulation_[rr + 1] - triangulation_[rr] * triangulation_[rr]) / 12.;
    }
    ret[dimRange - 1][dimRange - 1] = (3 * triangulation_[dimRange - 1] * triangulation_[dimRange - 1]
                                       - 2 * triangulation_[dimRange - 1] * triangulation_[dimRange - 2]
                                       - triangulation_[dimRange - 2] * triangulation_[dimRange - 2])
                                      / 12.;
    return ret;
  }

  // returns matrices with entries <v h_i h_j>_- and <v h_i h_j>_+
  virtual FieldVector<FieldVector<MatrixType, 2>, 1> kinetic_flux_matrices() const override final
  {
    FieldVector<FieldVector<MatrixType, 2>, 1> ret(FieldVector<MatrixType, 2>(MatrixType(dimRange, dimRange, 0.)));
    auto mm_with_v = mass_matrix_with_v();
    auto& ret_neg = ret[0][0];
    auto& ret_pos = ret[0][1];
    size_t N = dimRange;
    for (size_t nn = 0; nn < N; ++nn) {
      for (size_t mm = (nn > 0 ? nn - 1 : 0); mm <= (nn < N - 1 ? nn + 1 : nn); ++mm) {
        if (N % 2) {
          if (nn < N / 2 || mm < N / 2)
            ret_neg[nn][mm] = mm_with_v[0][nn][mm];
          else if (nn > N / 2 || mm > N / 2)
            ret_pos[nn][mm] = mm_with_v[0][nn][mm];
          else { // nn == mm == N/2
            ret_neg[nn][mm] = -std::pow(triangulation_[mm - 1], 2) / 12.;
            ret_pos[nn][mm] = std::pow(triangulation_[mm + 1], 2) / 12.;
          }
        } else {
          if (nn < N / 2 - 1 || mm < N / 2 - 1)
            ret_neg[nn][mm] = mm_with_v[0][nn][mm];
          else if (nn > N / 2 || mm > N / 2)
            ret_pos[nn][mm] = mm_with_v[0][nn][mm];
          else if (nn == N / 2 && mm == nn) {
            ret_neg[nn][mm] = -std::pow(triangulation_[mm], 2) / 48.;
            ret_pos[nn][mm] = (5 * std::pow(triangulation_[mm], 2) + 8 * triangulation_[mm] * triangulation_[mm + 1]
                               + 4 * std::pow(triangulation_[mm + 1], 2))
                              / 48.;
          } else if (nn == N / 2 - 1 && mm == nn) {
            ret_neg[nn][mm] = (-5 * std::pow(triangulation_[mm], 2) - 8 * triangulation_[mm] * triangulation_[mm - 1]
                               - 4 * std::pow(triangulation_[mm - 1], 2))
                              / 48.;
            ret_pos[nn][mm] = std::pow(triangulation_[mm], 2) / 48.;
          } else { // (((mm == N / 2 && nn == N / 2 - 1) || (mm == N / 2 - 1 && nn == N / 2))) {
            ret_neg[nn][mm] = -std::pow(triangulation_[mm], 2) / 16.;
            ret_pos[nn][mm] = std::pow(triangulation_[mm], 2) / 16.;
          }
        } // else (N % 2)
      } // mm
    } // nn
    return ret;
  }

  virtual MatrixType reflection_matrix(const DomainType& n) const override final
  {
    MatrixType ret(dimRange, dimRange, 0);
    for (size_t ii = 0; ii < dimDomain; ++ii)
      if (XT::Common::FloatCmp::ne(n[ii], 0.))
        if (XT::Common::FloatCmp::ne(std::abs(n[ii]), 1.))
          DUNE_THROW(NotImplemented, "Implemented only for +-e_i where e_i is the i-th canonical basis vector!");
    const auto mass_mat = mass_matrix();
    for (size_t ii = 0; ii < dimRange; ++ii)
      for (size_t jj = 0; jj < dimRange; ++jj)
        ret[ii][jj] = mass_mat[ii][dimRange - 1 - jj];
    ret.rightmultiply(mass_matrix_inverse());
#ifndef NDEBUG
    for (size_t ii = 0; ii < dimRange; ++ii)
      for (size_t jj = 0; jj < dimRange; ++jj)
        if (std::isnan(ret[ii][jj]) || std::isinf(ret[ii][jj]))
          DUNE_THROW(Dune::MathError,
                     "Calculation of reflection matrix failed for normal n = " + XT::Common::to_string(n));
#endif
    return ret;
  }

  template <class DiscreteFunctionType>
  VisualizerType<DiscreteFunctionType> visualizer() const
  {
    return [](const DiscreteFunctionType& u_n, const std::string& filename_prefix, const size_t ii) {
      sum_visualizer(u_n, filename_prefix, ii);
    };
  }

  static StringifierType stringifier()
  {
    return [](const RangeType& val) {
      RangeFieldType psi(0);
      for (const auto& entry : val)
        psi += entry;
      return XT::Common::to_string(psi, 15);
    };
  } // ... stringifier()

  const TriangulationType& triangulation() const
  {
    return triangulation_;
  }

  virtual RangeType alpha_iso() const override final
  {
    return RangeType(1.);
  }

  virtual RangeFieldType density(const RangeType& u) const override final
  {
    return std::accumulate(u.begin(), u.end(), RangeFieldType(0));
  }

  virtual RangeFieldType density_min(const RangeType& u) const override final
  {
    return *std::min_element(u.begin(), u.end());
  }

  virtual size_t density_factor() const override final
  {
    return dimRange;
  }

  virtual std::string short_id() const override final
  {
    return "1dhf";
  }

  // get indices of all faces that contain point v
  std::vector<size_t> get_face_indices(const DomainType& v) const
  {
    std::vector<size_t> face_indices;
    for (size_t jj = 0; jj < triangulation_.size() - 1; ++jj)
      if (XT::Common::FloatCmp::ge(v[0], triangulation_[jj]) && XT::Common::FloatCmp::le(v[0], triangulation_[jj + 1]))
        face_indices.push_back(jj);
    assert(face_indices.size());
    return face_indices;
  }

private:
  const TriangulationType triangulation_;
}; // class HatFunctionMomentBasis<DomainFieldType, 1, ...>

template <class DomainFieldType, class RangeFieldType, size_t refinements, size_t fluxDim>
class HatFunctionMomentBasis<DomainFieldType, 3, RangeFieldType, refinements, 1, fluxDim>
    : public BasisfunctionsInterface<DomainFieldType,
                                     3,
                                     RangeFieldType,
                                     OctaederStatistics<refinements>::num_vertices(),
                                     1,
                                     fluxDim>
{
public:
  static constexpr size_t dimDomain = 3;
  static constexpr size_t dimRange = OctaederStatistics<refinements>::num_vertices();
  static constexpr size_t dimRangeCols = 1;
  static constexpr size_t dimFlux = fluxDim;

private:
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols, dimFlux> BaseType;

public:
  typedef SphericalTriangulation<DomainFieldType> TriangulationType;
  using typename BaseType::DomainType;
  using typename BaseType::MatrixType;
  using typename BaseType::QuadraturesType;
  using typename BaseType::RangeType;
  using typename BaseType::StringifierType;
  template <class DiscreteFunctionType>
  using VisualizerType = typename BaseType::template VisualizerType<DiscreteFunctionType>;

  using BaseType::barycentre_rule;

  HatFunctionMomentBasis(const QuadraturesType& quadratures)
    : BaseType(refinements, quadratures)
  {
    assert(triangulation_.vertices().size() == dimRange);
  }

  HatFunctionMomentBasis(const size_t quad_refinements,
                         const QuadratureRule<RangeFieldType, 2>& reference_quadrature_rule)
    : BaseType(refinements)
  {
    quadratures_ = triangulation_.quadrature_rules(quad_refinements, reference_quadrature_rule);
    assert(triangulation_.vertices().size() == dimRange);
  }

  // This constructor is here for compatibility with the one-dimensional basis to simplify testing
  HatFunctionMomentBasis(const size_t fekete_rule_num = 3,
                         const size_t quad_refinements =
#if HAVE_FEKETE
                             0
#else
                             7
#endif
                         )
    : BaseType(refinements)
  {
#if HAVE_FEKETE
    const QuadratureRule<RangeFieldType, 2> reference_quadrature_rule =
        FeketeQuadrature<DomainFieldType>::get(fekete_rule_num);
#else
    DUNE_UNUSED_PARAMETER(fekete_rule_num);
    const QuadratureRule<RangeFieldType, 2> reference_quadrature_rule = barycentre_rule();
#endif
    quadratures_ = triangulation_.quadrature_rules(quad_refinements, reference_quadrature_rule);
    assert(triangulation_.vertices().size() == dimRange);
  }

  virtual RangeType evaluate(const DomainType& v) const override
  {
    RangeType ret(0);
    bool success = false;
    // walk over faces
    for (const auto& face : triangulation_.faces()) {
      const auto& vertices = face->vertices();
      DomainType barycentric_coords(0);
      success = calculate_barycentric_coordinates(v, vertices, barycentric_coords);
      if (success) {
        for (size_t ii = 0; ii < 3; ++ii)
          ret[vertices[ii]->index()] = barycentric_coords[ii];
        break;
      }
    } // faces
    assert(success);
    return ret;
  } // ... evaluate(...)

  virtual RangeType evaluate(const DomainType& v, const size_t face_index) const override final
  {
    RangeType ret(0);
    const auto& face = triangulation_.faces()[face_index];
    const auto& vertices = face->vertices();
    DomainType barycentric_coords(0);
    bool success = calculate_barycentric_coordinates(v, vertices, barycentric_coords);
    assert(success);
#ifdef NDEBUG
    static_cast<void>(success);
#endif
    for (size_t ii = 0; ii < 3; ++ii)
      ret[vertices[ii]->index()] = barycentric_coords[ii];
    return ret;
  } // ... evaluate(...)

  template <class DiscreteFunctionType>
  VisualizerType<DiscreteFunctionType> visualizer() const
  {
    return [](const DiscreteFunctionType& u_n, const std::string& filename_prefix, const size_t ii) {
      sum_visualizer(u_n, filename_prefix, ii);
    };
  }

  static StringifierType stringifier()
  {
    return [](const RangeType& val) {
      RangeFieldType psi(0);
      for (const auto& entry : val)
        psi += entry;
      return XT::Common::to_string(psi, 15);
    };
  } // ... stringifier()

  const TriangulationType& triangulation() const
  {
    return triangulation_;
  }

  virtual RangeFieldType unit_ball_volume() const override final
  {
    return BaseType::unit_ball_volume_quad();
  }

  virtual RangeType alpha_iso() const override final
  {
    return RangeType(1.);
  }

  virtual RangeFieldType density(const RangeType& u) const override final
  {
    return std::accumulate(u.begin(), u.end(), RangeFieldType(0));
  }

  virtual RangeFieldType density_min(const RangeType& u) const override final
  {
    return *std::min_element(u.begin(), u.end());
  }

  virtual size_t density_factor() const override final
  {
    return dimRange;
  }

  virtual std::string short_id() const override final
  {
    return "3dhf";
  }

  // get indices of all faces that contain point v
  std::vector<size_t> get_face_indices(const DomainType& v) const
  {
    return triangulation_.get_face_indices(v);
  }

  // calculates <b(v) dirac(v-dirac_position)>
  RangeType integrate_dirac_at(const DomainType& dirac_position) const
  {
    return evaluate(dirac_position);
  }

protected:
  using BaseType::parallel_quadrature;
  using BaseType::integrated_initializer;

  template <class VertexVectorType>
  bool calculate_barycentric_coordinates(const DomainType& v, const VertexVectorType& vertices, DomainType& ret) const
  {
    if (XT::Common::FloatCmp::ne(v.two_norm2(), 1.))
      DUNE_THROW(Dune::MathError, "Wrong input given!");
    Dune::FieldMatrix<RangeFieldType, 3, 3> gradients(0);
    for (size_t ii = 0; ii < 3; ++ii) {
      // copy vertices to gradients
      gradients[ii] = vertices[ii]->position();
      const auto scalar_prod = v * gradients[ii];
      // if v is not on the same half space of the sphere as the vertices, return false
      // assumes the triangulation is fine enough that vertices[ii]*vertices[jj] >= 0 for all triangles
      if (XT::Common::FloatCmp::lt(scalar_prod, 0.))
        return false;
      else if (XT::Common::FloatCmp::ge(scalar_prod, 1.)) {
        ret *= 0.;
        ret[ii] = 1.;
        return true;
      }
      auto v_scaled = v;
      v_scaled *= scalar_prod;
      gradients[ii] -= v_scaled;
      // scale with factor
      auto denominator = std::sqrt(1. - std::pow(scalar_prod, 2));
      if (std::isnan(denominator))
        DUNE_THROW(Dune::MathError, "NaN in evaluation!");
      if (std::isnan(std::acos(scalar_prod)))
        DUNE_THROW(Dune::MathError, "wrong value of scalar_prod!");
      gradients[ii] *= XT::Common::FloatCmp::eq(denominator, 0.) ? 0. : std::acos(scalar_prod) / denominator;
    } // ii
    // Calculate barycentric coordinates for 0 w.r.t to the points g_i = gradients[i]
    // For that purpose, solve the overdetermined system  A (h0 h1)^T = b
    // for the matrix A = (g_0-g_2 g_1-g_2) and the right-hand side b = -g_2.
    // The solution is (A^T A)^{-1} A^T b.
    // The third coordinate is calculated from the condition h0+h1+h2=1.
    Dune::XT::Common::FieldMatrix<RangeFieldType, 3, 2> A;
    Dune::XT::Common::FieldMatrix<RangeFieldType, 2, 3> AT;
    Dune::XT::Common::FieldVector<RangeFieldType, 2> solution;
    AT[0] = gradients[0];
    AT[1] = gradients[1];
    AT[0] -= gradients[2];
    AT[1] -= gradients[2];
    for (size_t ii = 0; ii < 3; ++ii)
      for (size_t jj = 0; jj < 2; ++jj)
        A[ii][jj] = AT[jj][ii];
    Dune::XT::Common::FieldMatrix<RangeFieldType, 2, 2> AT_A = AT.rightmultiplyany(A);
    gradients[2] *= -1;
    Dune::XT::Common::FieldVector<RangeFieldType, 2> AT_b;
    AT.mv(gradients[2], AT_b);
    AT_A.solve(solution, AT_b);
    ret[0] = solution[0];
    ret[1] = solution[1];
    ret[2] = 1. - ret[0] - ret[1];
    if (XT::Common::FloatCmp::lt(ret[0], 0.) || XT::Common::FloatCmp::lt(ret[1], 0.)
        || XT::Common::FloatCmp::lt(ret[2], 0.))
      return false;
    return true;
  } // bool calculate_barycentric_coordinates(...)

  using BaseType::quadratures_;
  using BaseType::triangulation_;
}; // class HatFunctionMomentBasis<DomainFieldType, 3, ...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_HATFUNCTIONS_HH
