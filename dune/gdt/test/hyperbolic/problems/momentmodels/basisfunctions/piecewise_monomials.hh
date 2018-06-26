// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_PIECEWISEMONOMIALS_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_PIECEWISEMONOMIALS_HH

#include "base.hh"

#include <dune/xt/common/fvector.hh>

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
  typedef FieldVector<DomainFieldType, dimRange / 2 + 1> TriangulationType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::MatrixType;
  using typename BaseType::StringifierType;
  template <class DiscreteFunctionType>
  using VisualizerType = typename BaseType::template VisualizerType<DiscreteFunctionType>;

  static std::string static_id()
  {
    return "pcw";
  }

  PiecewiseMonomials(const TriangulationType& triangulation = create_triangulation(),
                     const QuadratureType& /*quadrature*/ = QuadratureType())
    : triangulation_(triangulation)
  {
  }

  static TriangulationType create_triangulation()
  {
    TriangulationType ret;
    for (size_t ii = 0; ii < dimRange / 2 + 1; ++ii)
      ret[ii] = -1. + (4. * ii) / dimRange;
    return ret;
  }

  virtual RangeType evaluate(const DomainType& v) const override final
  {
    size_t dummy;
    return evaluate(v, false, dummy);
  } // ... evaluate(...)

  virtual RangeType evaluate(const DomainType& v, bool split_boundary, size_t& /*num_faces*/) const
  {
    RangeType ret(0);
    bool boundary = false;
    for (size_t ii = 0; ii < dimRange / 2; ++ii) {
      if (XT::Common::FloatCmp::eq(v[0], triangulation_[ii]))
        boundary = true;
      if (XT::Common::FloatCmp::ge(v[0], triangulation_[ii])
          && XT::Common::FloatCmp::le(v[0], triangulation_[ii + 1])) {
        ret[2 * ii] = 1;
        ret[2 * ii + 1] = v[0];
      }
    }
    if (XT::Common::FloatCmp::eq(v[0], triangulation_[dimRange / 2]))
      boundary = true;
    if (split_boundary && boundary)
      ret /= 2.;
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

  // returns matrices with entries <v h_i h_j>_- and <v h_i h_j>_+
  virtual FieldVector<FieldVector<MatrixType, 2>, 1> kinetic_flux_matrices() const
  {
    FieldVector<FieldVector<MatrixType, 2>, 1> ret(FieldVector<MatrixType, 2>(MatrixType(dimRange, dimRange, 0.)));
    auto mm_with_v = mass_matrix_with_v();
    auto& ret_neg = ret[0][0];
    auto& ret_pos = ret[0][1];
    for (size_t ii = 0; ii < dimRange / 2; ++ii) {
      if (dimRange / 2 % 2) {
        // if there is an odd number of intervals, the middle interval is split in 2 parts
        // copy other intervals
        for (size_t nn = 0; nn < dimRange / 2 - 1; ++nn)
          for (size_t mm = 0; mm < dimRange / 2 - 1; ++mm)
            ret_neg[nn][mm] = mm_with_v[0][nn][mm];
        for (size_t nn = dimRange / 2 + 1; nn < dimRange; ++nn)
          for (size_t mm = dimRange / 2 + 1; mm < dimRange; ++mm)
            ret_pos[nn][mm] = mm_with_v[0][nn][mm];
        // treat middle interval
        // mixed integrals are symmetric, so ret_pos and ret_neg both get half of it
        ret_neg[dimRange / 2 - 1][dimRange / 2] = mm_with_v[0][dimRange / 2 - 1][dimRange / 2] / 2;
        ret_neg[dimRange / 2][dimRange / 2 - 1] = mm_with_v[0][dimRange / 2][dimRange / 2 - 1] / 2;
        ret_pos[dimRange / 2 - 1][dimRange / 2] = mm_with_v[0][dimRange / 2 - 1][dimRange / 2] / 2;
        ret_pos[dimRange / 2][dimRange / 2 - 1] = mm_with_v[0][dimRange / 2][dimRange / 2 - 1] / 2;
        // integral corresponding to constant basis function
        ret_neg[dimRange / 2 - 1][dimRange / 2 - 1] = -std::pow(triangulation_[dimRange / 4], 2) / 2;
        ret_pos[dimRange / 2 - 1][dimRange / 2 - 1] = std::pow(triangulation_[dimRange / 4], 2) / 2;
        // integral corresponding to v basis function
        ret_neg[dimRange / 2][dimRange / 2] = -std::pow(triangulation_[dimRange / 4], 4) / 4;
        ret_pos[dimRange / 2][dimRange / 2] = std::pow(triangulation_[dimRange / 4], 4) / 4;
      } else {
        // if there is an even number of intervals, the matrix is just split up in upper and lower part
        for (size_t nn = 0; nn < dimRange / 2; ++nn)
          for (size_t mm = 0; mm < dimRange / 2; ++mm)
            ret_neg[nn][mm] = mm_with_v[0][nn][mm];
        for (size_t nn = dimRange / 2; nn < dimRange; ++nn)
          for (size_t mm = dimRange / 2; mm < dimRange; ++mm)
            ret_pos[nn][mm] = mm_with_v[0][nn][mm];
      }
    } // nn
    std::cout << "pos: " << XT::Common::to_string(ret_pos) << std::endl;
    std::cout << "neg: " << XT::Common::to_string(ret_neg) << std::endl;
    return ret;
  }

  virtual MatrixType reflection_matrix(const DomainType& n) const
  {
    MatrixType ret(dimRange, dimRange, 0);
    for (size_t ii = 0; ii < dimDomain; ++ii)
      if (XT::Common::FloatCmp::ne(n[ii], 0.))
        if (XT::Common::FloatCmp::ne(std::abs(n[ii]), 1.))
          DUNE_THROW(NotImplemented, "Implemented only for +-e_i where e_i is the i-th canonical basis vector!");
    const auto mass_mat = mass_matrix();
    for (size_t ii = 0; ii < dimRange; ++ii) {
      for (size_t jj = 0; jj < dimRange / 2; ++jj) {
        ret[ii][2 * jj] = mass_mat[ii][dimRange - 1 - (2 * jj + 1)];
        ret[ii][2 * jj + 1] = -mass_mat[ii][dimRange - 1 - 2 * jj];
      }
    }
    ret.rightmultiply(mass_matrix_inverse());
    return ret;
  }

  template <class DiscreteFunctionType>
  VisualizerType<DiscreteFunctionType> visualizer() const
  {
    return [](const DiscreteFunctionType& u_n, const std::string& filename_prefix, const size_t ii) {
      sum_divisible_by_visualizer<DiscreteFunctionType, dimRange>(u_n, filename_prefix, ii, 2);
    };
  }

  RangeFieldType calculate_psi_from_moments(const RangeType& val) const
  {
    RangeFieldType psi(0);
    for (size_t rr = 0; rr < dimRange; rr += 2)
      psi += val[rr];
    return psi;
  }

  static StringifierType stringifier()
  {
    return [](const RangeType& val) {
      RangeFieldType psi(0);
      for (size_t ii = 0; ii < dimRange; ii += 2)
        psi += val[ii];
      return XT::Common::to_string(psi, 15);
    };
  } // ... stringifier()

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
    RangeFieldType u_sum(0.);
    RangeFieldType u_bar_sum(0.);
    for (size_t ii = 0; ii < u.size(); ii += 2) {
      u_sum += u[ii];
      u_bar_sum += u_bar[ii];
    }
    return 2 * std::max(u_sum, u_bar_sum);
    //    return 2 * std::max(u[0], u_bar[0]);
  }

  // get indices of all faces that contain point
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
  using typename BaseType::StringifierType;
  template <class DiscreteFunctionType>
  using VisualizerType = typename BaseType::template VisualizerType<DiscreteFunctionType>;

  using BaseType::barycentre_rule;

  PiecewiseMonomials(const TriangulationType& triangulation, const QuadratureType& quadrature)
    : triangulation_(triangulation)
    , quadrature_(quadrature)
  {
    assert(4 * triangulation_.faces().size() == dimRange);
    calculate_plane_equations();
  }

  PiecewiseMonomials(const size_t refinements = 0,
                     const size_t quadrature_refinements = 4,
                     const QuadratureRule<RangeFieldType, 2>& reference_quadrature_rule = barycentre_rule(),
                     std::vector<Dune::XT::Common::FieldVector<DomainFieldType, dimDomain>> initial_points =
                         {{1., 0., 0.}, {-1., 0., 0.}, {0., 1., 0.}, {0., -1., 0.}, {0., 0., 1.}, {0., 0., -1.}})
    : triangulation_(initial_points, refinements, reference_quadrature_rule)
    , quadrature_(triangulation_.quadrature_rule(quadrature_refinements))
  {
    assert(4 * triangulation_.faces().size() == dimRange);
    calculate_plane_equations();
  }

  PiecewiseMonomials(const size_t refinements,
                     const QuadratureType& quadrature,
                     std::vector<Dune::XT::Common::FieldVector<DomainFieldType, dimDomain>> initial_points =
                         {{1., 0., 0.}, {-1., 0., 0.}, {0., 1., 0.}, {0., -1., 0.}, {0., 0., 1.}, {0., 0., -1.}})
    : triangulation_(initial_points, refinements)
    , quadrature_(quadrature)
  {
    assert(4 * triangulation_.faces().size() == dimRange);
    calculate_plane_equations();
  }

  virtual RangeType evaluate(const DomainType& v) const override final
  {
    size_t dummy;
    return evaluate(v, true, dummy);
  } // ... evaluate(...)

  virtual RangeType evaluate(const DomainType& v, bool split_boundary, size_t& num_faces) const
  {
    RangeType ret(0);
    const auto face_indices = triangulation_.get_face_indices(v);
    num_faces = face_indices.size();
    assert(num_faces);
    for (const auto& face_index : face_indices) {
      ret[4 * face_index] = 1.;
      for (size_t ii = 1; ii < 4; ++ii) {
        assert(4 * face_index + ii < ret.size());
        ret[4 * face_index + ii] = v[ii - 1];
      }
    } // face_indices
    if (split_boundary)
      ret /= RangeFieldType(num_faces);
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

  virtual FieldVector<FieldVector<MatrixType, 2>, dimFlux> kinetic_flux_matrices() const
  {
    FieldVector<FieldVector<MatrixType, 2>, dimFlux> B_kinetic(
        FieldVector<MatrixType, 2>(MatrixType(dimRange, dimRange, 0.)));
    QuadratureType neg_quadrature;
    QuadratureType pos_quadrature;
    neg_quadrature.reserve(quadrature_.size());
    pos_quadrature.reserve(quadrature_.size());
    for (size_t dd = 0; dd < dimFlux; ++dd) {
      neg_quadrature.clear();
      pos_quadrature.clear();
      for (const auto& quad_point : quadrature_) {
        const auto& v = quad_point.position();
        const auto& weight = quad_point.weight();
        if (XT::Common::FloatCmp::eq(v[dd], 0.)) {
          neg_quadrature.emplace_back(v, weight / 2.);
          pos_quadrature.emplace_back(v, weight / 2.);
        } else if (v[dd] > 0.)
          pos_quadrature.emplace_back(v, weight);
        else
          neg_quadrature.emplace_back(v, weight);
      }
      parallel_quadrature(neg_quadrature, B_kinetic[dd][0], dd);
      parallel_quadrature(pos_quadrature, B_kinetic[dd][1], dd);
    }
    return B_kinetic;
  } // ... kinetic_flux_matrices()

  virtual MatrixType reflection_matrix(const DomainType& n) const
  {
    MatrixType ret(dimRange, dimRange, 0);
    size_t direction;
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      if (XT::Common::FloatCmp::ne(n[ii], 0.)) {
        direction = ii;
        if (XT::Common::FloatCmp::ne(std::abs(n[ii]), 1.))
          DUNE_THROW(NotImplemented, "Implemented only for +-e_i where e_i is the i-th canonical basis vector!");
      }
    }
    parallel_quadrature(quadrature_, ret, direction, true);
    ret.rightmultiply(mass_matrix_inverse());
    return ret;
  }

  template <class DiscreteFunctionType>
  VisualizerType<DiscreteFunctionType> visualizer() const
  {
    return [](const DiscreteFunctionType& u_n, const std::string& filename_prefix, const size_t ii) {
      sum_divisible_by_visualizer<DiscreteFunctionType, dimRange>(u_n, filename_prefix, ii, 4);
    };
  }

  RangeFieldType calculate_psi_from_moments(const RangeType& val) const
  {
    RangeFieldType psi(0);
    for (size_t rr = 0; rr < dimRange; rr += 4)
      psi += val[rr];
    return psi;
  }

  static StringifierType stringifier()
  {
    return [](const RangeType& val) {
      RangeFieldType psi(0);
      for (size_t ii = 0; ii < dimRange; ii += 4)
        psi += val[ii];
      return XT::Common::to_string(psi, 15);
    };
  } // ... stringifier()

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
    RangeFieldType u_sum(0.);
    RangeFieldType u_bar_sum(0.);
    for (size_t ii = 0; ii < u.size(); ii += 4) {
      u_sum += u[ii];
      u_bar_sum += u_bar[ii];
    }
    return 2 * std::max(u_sum, u_bar_sum);
  }

  std::vector<size_t> get_face_indices(const DomainType& v) const
  {
    return triangulation_.get_face_indices(v);
  }

  const TriangulationType& triangulation() const
  {
    return triangulation_;
  }

  const QuadratureType& quadrature() const
  {
    return quadrature_;
  }

  // calculates <b(v) dirac(v-dirac_position)>
  RangeType integrate_dirac_at(const DomainType& dirac_position) const
  {
    size_t dummy;
    return evaluate(dirac_position, true, dummy);
  }

  bool obviously_realizable(const RangeType& u) const
  {
    static const RangeFieldType tol = 1e-6;
    for (size_t jj = 0; jj < dimRange / 4; ++jj) {
      const size_t offset = 4 * jj;
      const auto& u0 = u[offset];
      if (u0 < tol)
        return false;
      FieldVector<RangeFieldType, 3> u_norm{u[offset + 1], u[offset + 2], u[offset + 3]};
      u_norm /= u0;
      if ((n_[jj] * u_norm < d_[jj] + tol) || (u_norm.two_norm() > 1 - tol))
        return false;
    } // jj
    return true;
  }

protected:
  using BaseType::parallel_quadrature;

  virtual void calculate_in_thread(const QuadratureType& quadrature,
                                   MatrixType& local_matrix,
                                   const size_t v_index,
                                   const std::vector<size_t>& indices,
                                   const bool reflecting = false) const override
  {
    for (const auto& jj : indices) {
      const auto& quad_point = quadrature[jj];
      const auto& v = quad_point.position();
      auto v_reflected = v;
      if (reflecting)
        v_reflected[v_index] *= -1.;
      size_t num_adjacent_faces;
      const auto basis_evaluated = evaluate(v, false, num_adjacent_faces);
      const auto basis_reflected = evaluate(v_reflected);
      const auto& weight = quad_point.weight();
      const auto factor = (reflecting || v_index == size_t(-1)) ? 1. : v[v_index];
      for (size_t kk = 0; kk < local_matrix.N(); kk += 4)
        for (size_t nn = kk; nn < kk + 4; ++nn)
          for (size_t mm = kk; mm < kk + 4; ++mm)
            local_matrix[nn][mm] += basis_evaluated[nn] * (reflecting ? basis_reflected[mm] : basis_evaluated[mm])
                                    * factor * weight / num_adjacent_faces;
    } // ii
  } // void calculate_in_thread(...)

  void calculate_plane_equations()
  {
    for (const auto& face : triangulation_.faces()) {
      const auto& vertices = face->vertices();
      const auto index = face->index();
      // calculate plane equation defined by three points
      const DomainType v0v1 = vertices[1]->position() - vertices[0]->position();
      const DomainType v0v2 = vertices[2]->position() - vertices[0]->position();
      auto& n = n_[index];
      n = XT::Common::cross_product(v0v1, v0v2);
      if (n * vertices[0]->position() < 0.)
        n *= -1.;
      n /= n.two_norm();
      d_[index] = n * vertices[0]->position();
    }
  }

  using BaseType::integrated_initializer;

  const TriangulationType triangulation_;
  const QuadratureType quadrature_;
  FieldVector<FieldVector<DomainFieldType, 3>, dimRange / 4> n_;
  FieldVector<DomainFieldType, dimRange / 4> d_;
}; // class PiecewiseMonomials<DomainFieldType, 3, ...>


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_PIECEWISEMONOMIALS_HH
