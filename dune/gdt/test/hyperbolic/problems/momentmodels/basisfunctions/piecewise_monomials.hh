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
  static constexpr size_t dimDomain = 1;
  static constexpr size_t dimRange = rangeDim;
  static constexpr size_t dimRangeCols = rangeDimCols;
  static_assert(!(dimRange % 2), "dimRange has to be even!");
  static constexpr size_t num_intervals = dimRange / 2;

private:
  typedef BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols> BaseType;

public:
  using typename BaseType::DomainType;
  using typename BaseType::MatrixType;
  using typename BaseType::QuadraturesType;
  using typename BaseType::RangeType;
  using typename BaseType::StringifierType;
  template <class DiscreteFunctionType>
  using VisualizerType = typename BaseType::template VisualizerType<DiscreteFunctionType>;
  using TriangulationType = typename BaseType::Triangulation1dType;

  static std::string static_id()
  {
    return "pcw";
  }

  PiecewiseMonomials(const QuadraturesType& quadratures)
    : BaseType(quadratures)
    , triangulation_(BaseType::create_1d_triangulation(num_intervals))
  {
  }

  PiecewiseMonomials(const size_t quad_order = 15, const size_t DXTC_DEBUG_ONLY(quad_refinements) = 0)
    : BaseType(BaseType::gauss_lobatto_quadratures(num_intervals, quad_order))
    , triangulation_(BaseType::create_1d_triangulation(num_intervals))
  {
    assert(quad_refinements == 0 && "Refinement of the quadrature intervals not implemented for this basis!");
  }

  virtual RangeType evaluate(const DomainType& v) const override final
  {
    size_t dummy;
    return evaluate(v, false, dummy);
  } // ... evaluate(...)

  // evaluate on interval ii
  virtual RangeType evaluate(const DomainType& v, const size_t ii) const override final
  {
    RangeType ret(0);
    ret[2 * ii] = 1;
    ret[2 * ii + 1] = v[0];
    return ret;
  } // ... evaluate(...)

  virtual RangeType evaluate(const DomainType& v, bool split_boundary, size_t& /*num_faces*/) const
  {
    RangeType ret(0);
    bool boundary = false;
    for (size_t ii = 0; ii < num_intervals; ++ii) {
      if (XT::Common::FloatCmp::eq(v[0], triangulation_[ii]))
        boundary = true;
      if (XT::Common::FloatCmp::ge(v[0], triangulation_[ii])
          && XT::Common::FloatCmp::le(v[0], triangulation_[ii + 1])) {
        ret[2 * ii] = 1;
        ret[2 * ii + 1] = v[0];
      }
    }
    if (XT::Common::FloatCmp::eq(v[0], triangulation_[num_intervals]))
      boundary = true;
    if (split_boundary && boundary)
      ret /= 2.;
    return ret;
  } // ... evaluate(...)

  virtual RangeType integrated() const override final
  {
    RangeType ret(0);
    for (size_t ii = 0; ii < num_intervals; ++ii) {
      ret[2 * ii] = triangulation_[ii + 1] - triangulation_[ii];
      ret[2 * ii + 1] = (std::pow(triangulation_[ii + 1], 2) - std::pow(triangulation_[ii], 2)) / 2.;
    }
    return ret;
  }

  // returns matrix with entries <h_i h_j>
  virtual MatrixType mass_matrix() const override final
  {
    MatrixType M(dimRange, dimRange, 0.);
    for (size_t ii = 0; ii < num_intervals; ++ii) {
      M[2 * ii][2 * ii] = triangulation_[ii + 1] - triangulation_[ii];
      M[2 * ii + 1][2 * ii + 1] = (std::pow(triangulation_[ii + 1], 3) - std::pow(triangulation_[ii], 3)) / 3.;
      M[2 * ii][2 * ii + 1] = (std::pow(triangulation_[ii + 1], 2) - std::pow(triangulation_[ii], 2)) / 2.;
      M[2 * ii + 1][2 * ii] = M[2 * ii][2 * ii + 1];
    }
    return M;
  }

  virtual MatrixType mass_matrix_inverse() const override final
  {
    return tridiagonal_matrix_inverse<RangeFieldType, dimRange>(mass_matrix());
  }

  // returns matrix with entries <v h_i h_j>
  virtual FieldVector<MatrixType, dimDomain> mass_matrix_with_v() const override final
  {
    MatrixType B(dimRange, dimRange, 0.);
    for (size_t ii = 0; ii < num_intervals; ++ii) {
      B[2 * ii][2 * ii] = (std::pow(triangulation_[ii + 1], 2) - std::pow(triangulation_[ii], 2)) / 2.;
      B[2 * ii + 1][2 * ii + 1] = (std::pow(triangulation_[ii + 1], 4) - std::pow(triangulation_[ii], 4)) / 4.;
      B[2 * ii][2 * ii + 1] = (std::pow(triangulation_[ii + 1], 3) - std::pow(triangulation_[ii], 3)) / 3.;
      B[2 * ii + 1][2 * ii] = B[2 * ii][2 * ii + 1];
    }
    return FieldVector<MatrixType, dimDomain>(B);
  }

  // returns matrices with entries <v h_i h_j>_- and <v h_i h_j>_+
  virtual FieldVector<FieldVector<MatrixType, 2>, 1> kinetic_flux_matrices() const override final
  {
    FieldVector<FieldVector<MatrixType, 2>, 1> ret(FieldVector<MatrixType, 2>(MatrixType(dimRange, dimRange, 0.)));
    auto mm_with_v = mass_matrix_with_v();
    auto& ret_neg = ret[0][0];
    auto& ret_pos = ret[0][1];
    for (size_t ii = 0; ii < num_intervals; ++ii) {
      if (num_intervals % 2) {
        // if there is an odd number of intervals, the middle interval is split in 2 parts
        // copy other intervals
        for (size_t nn = 0; nn < num_intervals - 1; ++nn)
          for (size_t mm = 0; mm < num_intervals - 1; ++mm)
            ret_neg[nn][mm] = mm_with_v[0][nn][mm];
        for (size_t nn = num_intervals + 1; nn < dimRange; ++nn)
          for (size_t mm = num_intervals + 1; mm < dimRange; ++mm)
            ret_pos[nn][mm] = mm_with_v[0][nn][mm];
        // treat middle interval
        // mixed integrals are symmetric, so ret_pos and ret_neg both get half of it
        ret_neg[num_intervals - 1][num_intervals] = mm_with_v[0][num_intervals - 1][num_intervals] / 2;
        ret_neg[num_intervals][num_intervals - 1] = mm_with_v[0][num_intervals][num_intervals - 1] / 2;
        ret_pos[num_intervals - 1][num_intervals] = mm_with_v[0][num_intervals - 1][num_intervals] / 2;
        ret_pos[num_intervals][num_intervals - 1] = mm_with_v[0][num_intervals][num_intervals - 1] / 2;
        // integral corresponding to constant basis function
        ret_neg[num_intervals - 1][num_intervals - 1] = -std::pow(triangulation_[num_intervals / 2], 2) / 2;
        ret_pos[num_intervals - 1][num_intervals - 1] = std::pow(triangulation_[num_intervals / 2], 2) / 2;
        // integral corresponding to v basis function
        ret_neg[num_intervals][num_intervals] = -std::pow(triangulation_[num_intervals / 2], 4) / 4;
        ret_pos[num_intervals][num_intervals] = std::pow(triangulation_[num_intervals / 2], 4) / 4;
      } else {
        // if there is an even number of intervals, the matrix is just split up in upper and lower part
        for (size_t nn = 0; nn < num_intervals; ++nn)
          for (size_t mm = 0; mm < num_intervals; ++mm)
            ret_neg[nn][mm] = mm_with_v[0][nn][mm];
        for (size_t nn = num_intervals; nn < dimRange; ++nn)
          for (size_t mm = num_intervals; mm < dimRange; ++mm)
            ret_pos[nn][mm] = mm_with_v[0][nn][mm];
      }
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
    for (size_t ii = 0; ii < dimRange; ++ii) {
      for (size_t jj = 0; jj < num_intervals; ++jj) {
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

  static StringifierType stringifier()
  {
    return [](const RangeType& val) {
      RangeFieldType psi(0);
      for (size_t ii = 0; ii < dimRange; ii += 2)
        psi += val[ii];
      return XT::Common::to_string(psi, 15);
    };
  } // ... stringifier()

  const TriangulationType& triangulation() const
  {
    return triangulation_;
  }

  virtual RangeType alpha_iso() const override final
  {
    RangeType ret(0.);
    for (size_t ii = 0; ii < dimRange; ii += 2)
      ret[ii] = 1.;
    return ret;
  }

  virtual RangeFieldType density(const RangeType& u) const override final
  {
    RangeFieldType ret(0.);
    for (size_t ii = 0; ii < dimRange; ii += 2) {
      ret += u[ii];
    }
    return ret;
  }

  RangeFieldType density(const XT::Common::BlockedFieldVector<RangeFieldType, num_intervals, 2>& u) const
  {
    RangeFieldType ret(0.);
    for (size_t jj = 0; jj < num_intervals; ++jj)
      ret += u.block(jj)[0];
    return ret;
  }

  virtual std::string short_id() const override final
  {
    return "1dpm";
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

template <class DomainFieldType, class RangeFieldType, size_t refinements, size_t dimFlux>
class PiecewiseMonomials<DomainFieldType, 3, RangeFieldType, refinements, 1, dimFlux, 1>
    : public BasisfunctionsInterface<DomainFieldType,
                                     3,
                                     RangeFieldType,
                                     OctaederStatistics<refinements>::num_faces() * 4,
                                     1,
                                     dimFlux>
{
public:
  static const size_t dimDomain = 3;
  static const size_t dimRange = OctaederStatistics<refinements>::num_faces() * 4;
  static const size_t dimRangeCols = 1;

private:
  using BaseType = BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols, dimFlux>;

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

  PiecewiseMonomials(const QuadraturesType& quadratures)
    : BaseType(quadratures)
  {
    triangulation_ = TriangulationType(refinements);
    assert(4 * triangulation_.faces().size() == dimRange);
  }

  PiecewiseMonomials(const size_t quad_refinements, const QuadratureRule<RangeFieldType, 2>& reference_quadrature_rule)
  {
    triangulation_ = TriangulationType(refinements, reference_quadrature_rule);
    quadratures_ = triangulation_.quadrature_rules(quad_refinements);
    assert(4 * triangulation_.faces().size() == dimRange);
  }

  // This constructor is here for compatibility with the one-dimensional basis to simplify testing
  PiecewiseMonomials(const size_t fekete_rule_num = 3,
                     const size_t quad_refinements =
#if HAVE_FEKETE
                         0
#else
                         7
#endif
                     )
  {
#if HAVE_FEKETE
    const QuadratureRule<RangeFieldType, 2> reference_quadrature_rule =
        FeketeQuadrature<DomainFieldType>::get(fekete_rule_num);
#else
    DUNE_UNUSED_PARAMETER(fekete_rule_num);
    const QuadratureRule<RangeFieldType, 2> reference_quadrature_rule = barycentre_rule();
#endif
    triangulation_ = TriangulationType(refinements, reference_quadrature_rule);
    quadratures_ = triangulation_.quadrature_rules(quad_refinements);
    assert(4 * triangulation_.faces().size() == dimRange);
  }

  virtual RangeType evaluate(const DomainType& v) const override final
  {
    size_t dummy;
    return evaluate(v, true, dummy);
  } // ... evaluate(...)

  // evaluate on spherical triangle face_index
  virtual RangeType evaluate(const DomainType& v, const size_t face_index) const override final
  {
    RangeType ret(0);
    ret[4 * face_index] = 1.;
    for (size_t ii = 1; ii < 4; ++ii) {
      assert(4 * face_index + ii < ret.size());
      ret[4 * face_index + ii] = v[ii - 1];
    }
    return ret;
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

  template <class DiscreteFunctionType>
  VisualizerType<DiscreteFunctionType> visualizer() const
  {
    return [](const DiscreteFunctionType& u_n, const std::string& filename_prefix, const size_t ii) {
      sum_divisible_by_visualizer<DiscreteFunctionType, dimRange>(u_n, filename_prefix, ii, 4);
    };
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

  virtual RangeType alpha_iso() const override final
  {
    RangeType ret(0.);
    for (size_t ii = 0; ii < dimRange; ii += 4)
      ret[ii] = 1.;
    return ret;
  }

  virtual RangeFieldType density(const RangeType& u) const override final
  {
    RangeFieldType ret(0.);
    for (size_t ii = 0; ii < dimRange; ii += 4)
      ret += u[ii];
    return ret;
  }

  RangeFieldType density(const XT::Common::BlockedFieldVector<RangeFieldType, dimRange / 4, 4>& u) const
  {
    RangeFieldType ret(0.);
    for (size_t jj = 0; jj < dimRange / 4; ++jj)
      ret += u.block(jj)[0];
    return ret;
  }

  virtual std::string short_id() const override final
  {
    return "3dpm";
  }

  std::vector<size_t> get_face_indices(const DomainType& v) const
  {
    return triangulation_.get_face_indices(v);
  }

  const TriangulationType& triangulation() const
  {
    return triangulation_;
  }

  // calculates <b(v) dirac(v-dirac_position)>
  RangeType integrate_dirac_at(const DomainType& dirac_position) const
  {
    size_t dummy;
    return evaluate(dirac_position, true, dummy);
  }

protected:
  using BaseType::quadratures_;
  using BaseType::triangulation_;
}; // class PiecewiseMonomials<DomainFieldType, 3, ...>


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_PIECEWISEMONOMIALS_HH
