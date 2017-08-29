// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_PIECEWISEMONOMIALS_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_PIECEWISEMONOMIALS_HH

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
          size_t dimFlux = dimDomain,
          size_t order = 1>
class PiecewiseMonomials
{
  //  static_assert(false, "Not implemented for this combination of dimension and order!");
};

template <class DomainFieldType, class RangeFieldType, size_t rangeDim, size_t rangeDimCols, size_t dimFlux>
class PiecewiseMonomials<DomainFieldType, 1, RangeFieldType, rangeDim, rangeDimCols, dimFlux, 1>
    : public BasisfunctionsInterface<DomainFieldType, 1, RangeFieldType, rangeDim, rangeDimCols, dimFlux>
{
public:
  static const size_t dimDomain = 1;
  static const size_t dimRange = rangeDim;
  static const size_t dimRangeCols = rangeDimCols;
  static_assert(!(dimRange % 2), "dimRange has to be even!");

private:
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols> BaseType;

public:
  typedef typename Dune::QuadratureRule<DomainFieldType, dimDomain> QuadratureType;
  typedef FieldVector<DomainFieldType, dimRange + 1> TriangulationType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::MatrixType;
  template <class DiscreteFunctionType>
  using VisualizerType = typename BaseType::template VisualizerType<DiscreteFunctionType>;

  PiecewiseMonomials(const TriangulationType& triangulation = create_triangulation(),
                     const QuadratureType& /*quadrature*/ = QuadratureType())
    : triangulation_(triangulation)
  {
  }

  static TriangulationType create_triangulation()
  {
    TriangulationType ret;
    for (size_t ii = 0; ii < dimRange / 2 + 1; ++ii)
      ret[ii] = -1. + 4. * ii / dimRange;
    return ret;
  }

  virtual RangeType evaluate(const DomainType& v) const override final
  {
    RangeType ret(0);
    for (size_t ii = 0; ii < dimRange / 2; ++ii) {
      if (XT::Common::FloatCmp::ge(v[0], triangulation_[ii])
          && XT::Common::FloatCmp::le(v[0], triangulation_[ii + 1])) {
        ret[2 * ii] = 1;
        ret[2 * ii + 1] = v[0];
      }
    }
    return ret;
  } // ... evaluate(...)

  virtual RangeType integrated() const override final
  {
    RangeType ret(0);
    for (size_t ii = 0; ii < dimRange / 2; ++ii) {
      ret[2 * ii] = triangulation_[ii + 1] - triangulation_[ii];
      ret[2 * ii + 1] = (std::pow(triangulation_[ii + 1], 2) - std::pow(triangulation_[ii], 2)) / 2.;
    }
    return ret;
  }

  // returns matrix with entries <h_i h_j>
  virtual MatrixType mass_matrix() const override
  {
    MatrixType M(dimRange, dimRange, 0.);
    for (size_t ii = 0; ii < dimRange / 2; ++ii) {
      M[2 * ii][2 * ii] = triangulation_[ii + 1] - triangulation_[ii];
      M[2 * ii + 1][2 * ii + 1] = (std::pow(triangulation_[ii + 1], 3) - std::pow(triangulation_[ii], 3)) / 3.;
      M[2 * ii][2 * ii + 1] = (std::pow(triangulation_[ii + 1], 2) - std::pow(triangulation_[ii], 2)) / 2.;
      M[2 * ii + 1][2 * ii] = M[2 * ii][2 * ii + 1];
    }
    return M;
  }

  virtual MatrixType mass_matrix_inverse() const override
  {
    return tridiagonal_matrix_inverse<RangeFieldType, dimRange>(mass_matrix());
  }

  // returns matrix with entries <v h_i h_j>
  virtual FieldVector<MatrixType, dimDomain> mass_matrix_with_v() const override
  {
    MatrixType B(dimRange, dimRange, 0.);
    for (size_t ii = 0; ii < dimRange / 2; ++ii) {
      B[2 * ii][2 * ii] = (std::pow(triangulation_[ii + 1], 2) - std::pow(triangulation_[ii], 2)) / 2.;
      B[2 * ii + 1][2 * ii + 1] = (std::pow(triangulation_[ii + 1], 4) - std::pow(triangulation_[ii], 4)) / 4.;
      B[2 * ii][2 * ii + 1] = (std::pow(triangulation_[ii + 1], 3) - std::pow(triangulation_[ii], 3)) / 3.;
      B[2 * ii + 1][2 * ii] = B[2 * ii][2 * ii + 1];
    }
    return FieldVector<MatrixType, dimDomain>(B);
  }

  template <class DiscreteFunctionType>
  VisualizerType<DiscreteFunctionType> visualizer() const
  {
    return [](const DiscreteFunctionType& u_n, const std::string& filename_prefix, const size_t ii) {
      sum_divisible_by_visualizer<DiscreteFunctionType, dimRange>(u_n, filename_prefix, ii, 2);
    };
  }

  std::pair<RangeType, RangeType> calculate_isotropic_distribution(const RangeType& u) const
  {
    RangeType alpha_iso(0);
    RangeFieldType psi_iso(0);
    for (size_t ii = 0; ii < dimRange; ii += 2) {
      psi_iso += u[ii];
      alpha_iso[ii] = 1;
    }
    psi_iso /= 2.;
    alpha_iso *= std::log(psi_iso);
    RangeType u_iso = integrated();
    u_iso *= psi_iso;
    return std::make_pair(u_iso, alpha_iso);
  }

  const TriangulationType& triangulation() const
  {
    return triangulation_;
  }

  RangeFieldType realizability_limiter_max(const RangeType& u, const RangeType& u_bar) const
  {
    RangeFieldType u_sum;
    auto u_bar_sum = u_sum;
    for (size_t ii = 0; ii < u.size(); ii += 4) {
      u_sum += u[ii];
      u_bar_sum += u_bar[ii];
    }
    return 2 * std::max(u_sum, u_bar_sum);
  }

private:
  const TriangulationType triangulation_;
}; // class PiecewiseMonomials<DomainFieldType, 1, ...>

template <class DomainFieldType, class RangeFieldType, size_t rangeDim, size_t rangeDimCols, size_t dimFlux>
class PiecewiseMonomials<DomainFieldType, 3, RangeFieldType, rangeDim, rangeDimCols, dimFlux, 1>
    : public BasisfunctionsInterface<DomainFieldType, 3, RangeFieldType, rangeDim, rangeDimCols, dimFlux>
{
public:
  static const size_t dimDomain = 3;
  static const size_t dimRange = rangeDim;
  static const size_t dimRangeCols = rangeDimCols;

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

  PiecewiseMonomials(const TriangulationType& triangulation, const QuadratureType& quadrature)
    : triangulation_(triangulation)
    , quadrature_(quadrature)
  {
    assert(4 * triangulation_.faces().size() == dimRange);
  }

  PiecewiseMonomials(const size_t refinements = 0,
                     const size_t quadrature_refinements = 4,
                     std::vector<Dune::XT::Common::FieldVector<DomainFieldType, dimDomain>> initial_points =
                         {{1., 0., 0.}, {-1., 0., 0.}, {0., 1., 0.}, {0., -1., 0.}, {0., 0., 1.}, {0., 0., -1.}})
    : triangulation_(initial_points, refinements)
    , quadrature_(triangulation_.quadrature_rule(quadrature_refinements))
  {
    assert(4 * triangulation_.faces().size() == dimRange);
  }

  PiecewiseMonomials(const size_t refinements,
                     const QuadratureType& quadrature,
                     std::vector<Dune::XT::Common::FieldVector<DomainFieldType, dimDomain>> initial_points =
                         {{1., 0., 0.}, {-1., 0., 0.}, {0., 1., 0.}, {0., -1., 0.}, {0., 0., 1.}, {0., 0., -1.}})
    : triangulation_(initial_points, refinements)
    , quadrature_(quadrature)
  {
    assert(4 * triangulation_.faces().size() == dimRange);
  }

  virtual RangeType evaluate(const DomainType& v) const override final
  {
    size_t dummy;
    return evaluate(v, true, dummy);
  } // ... evaluate(...)

  virtual RangeType evaluate(const DomainType& v, bool split_boundary, size_t& num_faces) const
  {
    RangeType ret(0);
    FieldMatrix<RangeFieldType, 3, 3> vertices_matrix;
    FieldMatrix<RangeFieldType, 3, 3> determinant_matrix;
    bool success = false;
    size_t num_adjacent_faces = 0;
    for (const auto& face : triangulation_.faces()) {
      bool v_in_this_facet = true;
      bool second_check = true;
      const auto& vertices = face->vertices();
      for (size_t ii = 0; ii < 3; ++ii) {
        // if v is not on the same octant of the sphere as the vertices, return false
        // assumes the triangulation is fine enough that vertices[ii]*vertices[jj] >= 0 for all triangles
        const auto scalar_prod = v * vertices[ii]->position();
        if (XT::Common::FloatCmp::lt(scalar_prod, 0.)) {
          v_in_this_facet = false;
          second_check = false;
          break;
        } else if (XT::Common::FloatCmp::eq(scalar_prod, 1.)) {
          ++num_adjacent_faces;
          second_check = false;
          break;
        }
        vertices_matrix[ii] = vertices[ii]->position();
      } // ii

      if (second_check) {
        // Vertices are ordered counterclockwise, so if the point is inside the spherical triangle,
        // the coordinate system formed by two adjacent vertices and v is always right-handed, i.e.
        // the triple product is positive.
        // The triple products that need to be positive are the determinants of the matrices (v1, v2, v), (v2, v3, v),
        // (v3, v1, v), where vi is the ith vertex. Swapping two columns changes the sign of det, the matrices used
        // below all have an even number of column swaps.
        // As we have checked before that v is in the same octant as the vertices, the determinant is 0 iff v is on
        // an edge of the facet. In that case, assign half of the basis function to this edge.
        for (size_t ii = 0; ii < 3; ++ii) {
          determinant_matrix = vertices_matrix;
          determinant_matrix[ii] = v;
          auto det = determinant_matrix.determinant();
          if (XT::Common::FloatCmp::eq(det, 0.)) {
            ++num_adjacent_faces;
            break;
          } else if (det < 0.) {
            v_in_this_facet = false;
            break;
          }
        } // ii
      } // if (second_check)
      if (v_in_this_facet) {
        const auto face_index = face->index();
        ret[4 * face_index] = 1.;
        for (size_t ii = 1; ii < 4; ++ii) {
          assert(4 * face_index + ii < ret.size());
          ret[4 * face_index + ii] = v[ii - 1];
        }
        success = true;
      }
    } // faces
    if (split_boundary && num_adjacent_faces > 0)
      ret /= RangeFieldType(num_adjacent_faces);
    num_faces = num_adjacent_faces > 0 ? num_adjacent_faces : 1;
    assert(success);
    return ret;
  } // ... evaluate(...)

  // returns <b>, where b is the basis functions vector
  virtual RangeType integrated() const override
  {
    static const RangeType ret = integrated_initializer(quadrature_);
    return ret;
  }

  virtual MatrixType mass_matrix() const override
  {
    MatrixType M(dimRange, dimRange, 0.);
    parallel_quadrature(quadrature_, M, size_t(-1));
    return M;
  } // ... mass_matrix()

  virtual MatrixType mass_matrix_inverse() const override
  {
    auto ret = mass_matrix();
    ret.invert();
    return ret;
  }

  virtual FieldVector<MatrixType, dimFlux> mass_matrix_with_v() const override
  {
    FieldVector<MatrixType, dimFlux> B(MatrixType(dimRange, dimRange, 0));
    for (size_t dd = 0; dd < dimFlux; ++dd)
      parallel_quadrature(quadrature_, B[dd], dd);
    return B;
  }

  template <class DiscreteFunctionType>
  VisualizerType<DiscreteFunctionType> visualizer() const
  {
    return [](const DiscreteFunctionType& u_n, const std::string& filename_prefix, const size_t ii) {
      sum_divisible_by_visualizer<DiscreteFunctionType, dimRange>(u_n, filename_prefix, ii, 4);
    };
  }

  std::pair<RangeType, RangeType> calculate_isotropic_distribution(const RangeType& u) const
  {
    RangeFieldType psi_iso(0);
    RangeType alpha_iso(0);
    for (size_t ii = 0; ii < dimRange; ii += 4) {
      psi_iso += u[ii];
      alpha_iso[ii] = 1.;
    }
    psi_iso /= 4. * M_PI;
    alpha_iso *= std::log(psi_iso);
    auto u_iso = integrated();
    u_iso *= psi_iso;
    return std::make_pair(u_iso, alpha_iso);
  }

  RangeFieldType realizability_limiter_max(const RangeType& u, const RangeType& u_bar) const
  {
    RangeFieldType u_sum(0);
    auto u_bar_sum = u_sum;
    for (size_t ii = 0; ii < u.size(); ii += 4) {
      u_sum += u[ii];
      u_bar_sum += u_bar[ii];
    }
    return 2 * std::max(u_sum, u_bar_sum);
  }

  const TriangulationType& triangulation() const
  {
    return triangulation_;
  }

  const QuadratureType& quadrature() const
  {
    return quadrature_;
  }

protected:
  using BaseType::parallel_quadrature;

  virtual void calculate_in_thread(const QuadratureType& quadrature,
                                   MatrixType& local_matrix,
                                   const size_t v_index,
                                   const std::vector<size_t>& indices) const override
  {
    for (const auto& jj : indices) {
      const auto& quad_point = quadrature[jj];
      const auto& v = quad_point.position();
      size_t num_adjacent_faces;
      const auto basis_evaluated = evaluate(v, false, num_adjacent_faces);
      const auto& weight = quad_point.weight();
      const auto factor = (v_index == size_t(-1)) ? 1. : v[v_index];
      for (size_t kk = 0; kk < local_matrix.N(); kk += 4)
        for (size_t nn = kk; nn < kk + 4; ++nn)
          for (size_t mm = kk; mm < kk + 4; ++mm)
            local_matrix[nn][mm] += basis_evaluated[nn] * basis_evaluated[mm] * factor * weight / num_adjacent_faces;
    } // ii
  } // void calculate_in_thread(...)

  using BaseType::integrated_initializer;

  const TriangulationType triangulation_;
  const QuadratureType quadrature_;
}; // class PiecewiseMonomials<DomainFieldType, 3, ...>


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_PIECEWISEMONOMIALS_HH
