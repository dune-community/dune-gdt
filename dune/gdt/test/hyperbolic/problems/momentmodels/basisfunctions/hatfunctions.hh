// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_HATFUNCTIONS_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_HATFUNCTIONS_HH

#include <vector>
#include <string>

#include <dune/gdt/test/hyperbolic/problems/momentmodels/triangulation.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


template <class DomainFieldType,
          size_t dimDomain,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols = 1,
          size_t dimFlux = dimDomain>
class HatFunctions
{
  //  static_assert(false, "Not implemented for this dimension!");
};

template <class DomainFieldType, class RangeFieldType, size_t rangeDim, size_t rangeDimCols, size_t fluxDim>
class HatFunctions<DomainFieldType, 1, RangeFieldType, rangeDim, rangeDimCols, fluxDim>
    : public BasisfunctionsInterface<DomainFieldType, 1, RangeFieldType, rangeDim, rangeDimCols, fluxDim>
{
public:
  static const size_t dimDomain = 1;
  static const size_t dimRange = rangeDim;
  static const size_t dimRangeCols = rangeDimCols;

private:
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols, fluxDim> BaseType;

public:
  typedef typename Dune::QuadratureRule<DomainFieldType, dimDomain> QuadratureType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::MatrixType;
  template <class DiscreteFunctionType>
  using VisualizerType = typename BaseType::template VisualizerType<DiscreteFunctionType>;
  typedef RangeType TriangulationType;

  HatFunctions(const TriangulationType triangulation = create_triangulation(),
               const QuadratureType& /*quadrature*/ = QuadratureType())
    : triangulation_(triangulation)
  {
  }

  static TriangulationType create_triangulation()
  {
    RangeType ret;
    for (size_t ii = 0; ii < dimRange; ++ii)
      ret[ii] = -1. + 2. * ii / (dimRange - 1.);
    return ret;
  }

  virtual RangeType evaluate(const DomainType& v) const override
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

  virtual RangeType integrated() const override
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
  virtual MatrixType mass_matrix() const override
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

  virtual MatrixType mass_matrix_inverse() const override
  {
    return tridiagonal_matrix_inverse<RangeFieldType, dimRange>(mass_matrix());
  }

  // returns matrix with entries <v h_i h_j>
  virtual FieldVector<MatrixType, 1> mass_matrix_with_v() const override
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

  template <class DiscreteFunctionType>
  VisualizerType<DiscreteFunctionType> visualizer() const
  {
    return [](const DiscreteFunctionType& u_n, const std::string& filename_prefix, const size_t ii) {
      sum_visualizer<DiscreteFunctionType, dimRange>(u_n, filename_prefix, ii);
    };
  }

  std::pair<RangeType, RangeType> calculate_isotropic_distribution(const RangeType& u) const
  {
    RangeFieldType psi_iso(0);
    for (size_t ii = 0; ii < dimRange; ++ii)
      psi_iso += u[ii];
    psi_iso /= 2.;
    RangeType alpha_iso(std::log(psi_iso)), u_iso;
    u_iso = integrated();
    u_iso *= psi_iso / 2.;
    return std::make_pair(u_iso, alpha_iso);
  }

  const TriangulationType& triangulation() const
  {
    return triangulation_;
  }

  RangeFieldType realizability_limiter_max(const RangeType& u, const RangeType& u_bar) const
  {
    return 2 * std::max(std::accumulate(u.begin(), u.end(), RangeFieldType(0)),
                        std::accumulate(u_bar.begin(), u_bar.end(), RangeFieldType(0)));
  }

private:
  const TriangulationType triangulation_;
}; // class HatFunctions<DomainFieldType, 1, ...>

template <class DomainFieldType, class RangeFieldType, size_t rangeDim, size_t rangeDimCols, size_t fluxDim>
class HatFunctions<DomainFieldType, 2, RangeFieldType, rangeDim, rangeDimCols, fluxDim>
    : public BasisfunctionsInterface<DomainFieldType, 2, RangeFieldType, rangeDim, rangeDimCols, fluxDim>
{
public:
  static const size_t dimDomain = 2;
  static const size_t dimRange = rangeDim;
  static const size_t dimRangeCols = rangeDimCols;

private:
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, rangeDimCols, fluxDim> BaseType;

public:
  typedef typename Dune::QuadratureRule<DomainFieldType, dimDomain> QuadratureType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::MatrixType;
  template <class DiscreteFunctionType>
  using VisualizerType = typename BaseType::template VisualizerType<DiscreteFunctionType>;
  typedef RangeType TriangulationType;

  HatFunctions(const TriangulationType triangulation = create_triangulation(),
               const QuadratureType& /*quadrature*/ = QuadratureType())
    : triangulation_(triangulation)
  {
  }

  static TriangulationType create_triangulation()
  {
    RangeType ret;
    for (size_t ii = 0; ii < dimRange; ++ii)
      ret[ii] = -1. + 2. * ii / (dimRange - 1.);
    return ret;
  }

  virtual RangeType evaluate(const DomainType& v) const override
  {
    RangeType ret(0);
    return ret;
  } // ... evaluate(...)

  virtual RangeType integrated() const override
  {
    RangeType ret(0);
    return ret;
  }

  // returns matrix with entries <h_i h_j>
  virtual MatrixType mass_matrix() const override
  {
    MatrixType ret(dimRange, dimRange, 0);
    return ret;
  }

  virtual MatrixType mass_matrix_inverse() const override
  {
    MatrixType ret(dimRange, dimRange, 0);
    return ret;
  }

  // returns matrix with entries <v h_i h_j>
  virtual FieldVector<MatrixType, fluxDim> mass_matrix_with_v() const override
  {
    MatrixType ret(dimRange, dimRange, 0);
    return ret;
  }

  template <class DiscreteFunctionType>
  VisualizerType<DiscreteFunctionType> visualizer() const
  {
    return [](const DiscreteFunctionType& u_n, const std::string& filename_prefix, const size_t ii) {
      sum_visualizer<DiscreteFunctionType, dimRange>(u_n, filename_prefix, ii);
    };
  }

  std::pair<RangeType, RangeType> calculate_isotropic_distribution(const RangeType& u) const
  {
    RangeFieldType psi_iso(0);
    for (size_t ii = 0; ii < dimRange; ++ii)
      psi_iso += u[ii];
    psi_iso /= 2.;
    RangeType alpha_iso(std::log(psi_iso)), u_iso;
    u_iso = integrated();
    u_iso *= psi_iso / 2.;
    return std::make_pair(u_iso, alpha_iso);
  }

  const TriangulationType& triangulation() const
  {
    return triangulation_;
  }

  RangeFieldType realizability_limiter_max(const RangeType& u, const RangeType& u_bar) const
  {
    return 2 * std::max(std::accumulate(u.begin(), u.end(), RangeFieldType(0)),
                        std::accumulate(u_bar.begin(), u_bar.end(), RangeFieldType(0)));
  }

private:
  const TriangulationType triangulation_;
}; // class HatFunctions<DomainFieldType, 2, ...>

template <class DomainFieldType, class RangeFieldType, size_t rangeDim, size_t rangeDimCols, size_t fluxDim>
class HatFunctions<DomainFieldType, 3, RangeFieldType, rangeDim, rangeDimCols, fluxDim>
    : public BasisfunctionsInterface<DomainFieldType, 3, RangeFieldType, rangeDim, rangeDimCols, fluxDim>
{
public:
  static const size_t dimDomain = 3;
  static const size_t dimRange = rangeDim;
  static const size_t dimRangeCols = rangeDimCols;
  static const size_t dimFlux = fluxDim;

private:
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols, dimFlux> BaseType;

public:
  typedef typename Dune::QuadratureRule<DomainFieldType, dimDomain> QuadratureType;
  typedef SphericalTriangulation<DomainFieldType> TriangulationType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::MatrixType;
  template <class DiscreteFunctionType>
  using VisualizerType = typename BaseType::template VisualizerType<DiscreteFunctionType>;

  HatFunctions(const TriangulationType& triangulation, const QuadratureType& quadrature)
    : triangulation_(triangulation)
    , quadrature_(quadrature)
  {
    assert(triangulation_.vertices().size() == dimRange);
  }

  HatFunctions(const size_t refinements = 0,
               const size_t quadrature_refinements = 4,
               std::vector<Dune::XT::Common::FieldVector<DomainFieldType, dimDomain>> initial_points =
                   {{1., 0., 0.}, {-1., 0., 0.}, {0., 1., 0.}, {0., -1., 0.}, {0., 0., 1.}, {0., 0., -1.}})
    : triangulation_(initial_points, refinements)
    , quadrature_(triangulation_.quadrature_rule(quadrature_refinements))
  {
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

  // avoid recalculation of integral by using a static local variable that is initialized on first call
  virtual RangeType integrated() const override
  {
    static const RangeType ret = integrated_initializer(quadrature_);
    return ret;
  }

  virtual MatrixType mass_matrix() const override
  {
    MatrixType A(dimRange, dimRange, 0);
    parallel_quadrature(quadrature_, A, size_t(-1));
    return A;
  } // ... mass_matrix()

  virtual MatrixType mass_matrix_inverse() const override
  {
    auto ret = mass_matrix();
    ret.invert();
    return ret;
  }

  virtual FieldVector<MatrixType, dimFlux> mass_matrix_with_v() const override
  {
    FieldVector<MatrixType, dimFlux> B(MatrixType(dimRange, dimRange, 0.));
    for (size_t dd = 0; dd < dimFlux; ++dd)
      parallel_quadrature(quadrature_, B[dd], dd);
    return B;
  } // ... mass_matrix_with_v()

  template <class DiscreteFunctionType>
  VisualizerType<DiscreteFunctionType> visualizer() const
  {
    return [](const DiscreteFunctionType& u_n, const std::string& filename_prefix, const size_t ii) {
      sum_visualizer<DiscreteFunctionType, dimRange>(u_n, filename_prefix, ii);
    };
  }

  std::pair<RangeType, RangeType> calculate_isotropic_distribution(const RangeType& u) const
  {
    RangeFieldType psi_iso(0);
    for (size_t ii = 0; ii < dimRange; ++ii)
      psi_iso += u[ii];
    psi_iso /= 4. * M_PI;
    RangeType alpha_iso(std::log(psi_iso));
    auto u_iso = integrated();
    u_iso *= psi_iso;
    return std::make_pair(u_iso, alpha_iso);
  }

  RangeFieldType realizability_limiter_max(const RangeType& u, const RangeType& u_bar) const
  {
    return 2 * std::max(std::accumulate(u.begin(), u.end(), RangeFieldType(0)),
                        std::accumulate(u_bar.begin(), u_bar.end(), RangeFieldType(0)));
  }

  const QuadratureType& quadrature() const
  {
    return quadrature_;
  }

protected:
  using BaseType::parallel_quadrature;
  using BaseType::integrated_initializer;

  template <class VertexVectorType>
  bool calculate_barycentric_coordinates(const DomainType& v, const VertexVectorType& vertices, DomainType& ret) const
  {
    Dune::FieldMatrix<RangeFieldType, 3, 3> gradients(0);
    for (size_t ii = 0; ii < 3; ++ii) {
      // copy vertices to gradients
      gradients[ii] = vertices[ii]->position();
      const auto scalar_prod = v * gradients[ii];
      // if v is not on the same octant of the sphere as the vertices, return false
      // assumes the triangulation is fine enough that vertices[ii]*vertices[jj] >= 0 for all triangles
      if (XT::Common::FloatCmp::lt(scalar_prod, 0.))
        return false;
      auto v_scaled = v;
      v_scaled *= scalar_prod;
      gradients[ii] -= v_scaled;
      // scale with factor
      auto denominator = std::sqrt(1. - std::pow(scalar_prod, 2));
      gradients[ii] *= XT::Common::FloatCmp::eq(denominator, 0.) ? 0. : std::acos(scalar_prod) / denominator;
    } // ii
    // Calculate barycentric coordinates for 0 w.r.t to the points g_i = gradients[i]
    // For that purpose, solve the overdetermined system  A (h0 h1)^T = b
    // for the matrix A = (g_0-g_2 g_1-g_2) and the right-hand side b = -g_2.
    // The solution is (A^T A)^{-1} A^T b.
    // The third coordinate is calculated from the condition h0+h1+h2=1.
    Dune::FieldMatrix<RangeFieldType, 3, 2> A;
    Dune::FieldMatrix<RangeFieldType, 2, 3> AT;
    Dune::FieldVector<RangeFieldType, 2> solution;
    AT[0] = gradients[0];
    AT[1] = gradients[1];
    AT[0] -= gradients[2];
    AT[1] -= gradients[2];
    for (size_t ii = 0; ii < 3; ++ii)
      for (size_t jj = 0; jj < 2; ++jj)
        A[ii][jj] = AT[jj][ii];
    Dune::FieldMatrix<RangeFieldType, 2, 2> AT_A = AT.rightmultiplyany(A);
    gradients[2] *= -1;
    FieldVector<RangeFieldType, 2> AT_b;
    AT.mv(gradients[2], AT_b);
    AT_A.solve(solution, AT_b);
    ret[0] = solution[0];
    ret[1] = solution[1];
    ret[2] = 1. - ret[0] - ret[1];
    if (XT::Common::FloatCmp::lt(ret[0], 0.) || XT::Common::FloatCmp::lt(ret[1], 0.))
      return false;
    if (XT::Common::FloatCmp::lt(ret[2], 0.))
      return false;
    return true;
  } // bool calculate_barycentric_coordinates(...)

  const TriangulationType triangulation_;
  const QuadratureType quadrature_;
}; // class HatFunctions<DomainFieldType, 3, ...>


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_HATFUNCTIONS_HH
