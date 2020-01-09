// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_PARTIALMOMENTS_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_PARTIALMOMENTS_HH

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/null.hpp>
// #include <boost/math/special_functions/lambert_w.hpp>

#if HAVE_QHULL
#  include <dune/xt/common/disable_warnings.hh>
#  include <libqhullcpp/Qhull.h>
#  include <libqhullcpp/QhullFacetList.h>
#  include <dune/xt/common/reenable_warnings.hh>
#endif // HAVE_QHULL

#include <dune/xt/common/fvector.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {


template <class DomainFieldType,
          size_t domainDim,
          class RangeFieldType,
          size_t rangeDim,
          size_t dimFlux,
          EntropyType entropy,
          size_t sizeBlock>
class PartialMomentBasisBase
  : public MomentBasisInterface<DomainFieldType, domainDim, RangeFieldType, rangeDim, 1, dimFlux, entropy>
{
  using BaseType = MomentBasisInterface<DomainFieldType, domainDim, RangeFieldType, rangeDim, 1, dimFlux, entropy>;
  using ThisType = PartialMomentBasisBase;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  static constexpr size_t block_size = sizeBlock;
  static constexpr size_t num_blocks = dimRange / block_size;

  using typename BaseType::DomainType;
  using typename BaseType::DynamicRangeType;
  using typename BaseType::RangeType;
  using typename BaseType::StringifierType;
  using LocalVectorType = XT::Common::FieldVector<RangeFieldType, block_size>;

  template <class... Args>
  PartialMomentBasisBase(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}

  // evaluate on spherical triangle face_index
  DynamicRangeType evaluate(const DomainType& v, const size_t face_index) const override final
  {
    DynamicRangeType ret(dimRange, 0);
    const auto local_eval = evaluate_on_face(v, face_index);
    for (size_t ii = 0; ii < block_size; ++ii)
      ret[block_size * face_index + ii] = local_eval[ii];
    return ret;
  } // ... evaluate(...)

  // evaluate on spherical triangle face_index
  LocalVectorType evaluate_on_face(const DomainType& v, const size_t /*face_index*/) const
  {
    LocalVectorType ret;
    ret[0] = 1.;
    for (size_t ii = 1; ii < block_size; ++ii)
      ret[ii] = v[ii - 1];
    return ret;
  } // ... evaluate(...)

  static StringifierType stringifier()
  {
    return [](const DynamicRangeType& val) {
      RangeFieldType psi(0);
      for (size_t ii = 0; ii < dimRange; ii += block_size)
        psi += val[ii];
      return XT::Common::to_string(psi, 15);
    };
  } // ... stringifier()

  DynamicRangeType alpha_one() const override final
  {
    DynamicRangeType ret(dimRange, 0.);
    for (size_t ii = 0; ii < dimRange; ii += block_size)
      ret[ii] = 1.;
    return ret;
  }

  RangeFieldType density(const RangeType& u) const override final
  {
    RangeFieldType ret(0.);
    for (size_t ii = 0; ii < dimRange; ii += block_size)
      ret += u[ii];
    return ret;
  }

  RangeFieldType density(const DynamicRangeType& u) const override final
  {
    RangeFieldType ret(0.);
    for (size_t ii = 0; ii < dimRange; ii += block_size)
      ret += u[ii];
    return ret;
  }

  RangeFieldType density(const XT::Common::BlockedFieldVector<RangeFieldType, num_blocks, block_size>& u) const
  {
    RangeFieldType ret(0.);
    for (size_t jj = 0; jj < num_blocks; ++jj)
      ret += u.block(jj)[0];
    return ret;
  }

  RangeFieldType min_density(const XT::Common::BlockedFieldVector<RangeFieldType, num_blocks, block_size>& u) const
  {
    RangeFieldType ret(u.block(0)[0]);
    for (size_t jj = 1; jj < num_blocks; ++jj)
      ret = std::min(ret, u.block(jj)[0]);
    return ret;
  }

  using BaseType::u_iso;

  // For the partial moments, we might not be able to solve the optimization problem for some moments where the density
  // on one interval/spherical triangle is very low. The overall density might be much higher than the density on that
  // triangle, so we specialize this function.
  void ensure_min_density(DynamicRangeType& u, const RangeFieldType min_density) const override final
  {
    const auto u_iso_min = u_iso() * min_density;
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      const auto block_density = u[block_size * jj];
      if (block_density < u_iso_min[block_size * jj]) {
        for (size_t ii = 0; ii < block_size; ++ii)
          u[block_size * jj + ii] = u_iso_min[block_size * jj + ii];
      }
    }
  }

  // For the partial moments, we might not be able to solve the optimization problem for some moments where the density
  // on one interval/spherical triangle is very low. The overall density might be much higher than the density on that
  // triangle, so we specialize this function.
  void ensure_min_density(RangeType& u, const RangeFieldType min_density) const override final
  {
    const auto u_iso_min = u_iso() * min_density;
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      const auto block_density = u[block_size * jj];
      if (block_density < u_iso_min[block_size * jj]) {
        for (size_t ii = 0; ii < block_size; ++ii)
          u[block_size * jj + ii] = u_iso_min[block_size * jj + ii];
      }
    }
  }

  using BaseType::has_fixed_sign;

  DynamicRangeType
  get_u(const std::function<RangeFieldType(DomainType, std::array<int, dimDomain>)>& psi) const override final
  {
    DynamicRangeType ret(dimRange, 0);
    for (size_t jj = 0; jj < quadratures_.size(); ++jj) {
      const size_t offset = jj * block_size;
      for (const auto& quad_point : quadratures_[jj]) {
        const auto& v = quad_point.position();
        const auto basis_val = evaluate_on_face(v, jj);
        const auto psi_val = psi(v, has_fixed_sign(jj));
        const auto factor = psi_val * quad_point.weight();
        for (size_t ii = 0; ii < block_size; ++ii)
          ret[offset + ii] += basis_val[ii] * factor;
      }
    }
    return ret;
  }

  std::string short_id() const override final
  {
    return "pm";
  }

  std::string mn_name() const override final
  {
    return "pmm" + XT::Common::to_string(dimRange);
  }

  std::string pn_name() const override final
  {
    return "pmp" + XT::Common::to_string(dimRange);
  }

protected:
  using BaseType::quadratures_;
}; // class PartialMomentBasisBase<...>

template <class DomainFieldType,
          size_t domainDim,
          class RangeFieldType,
          size_t rangeDim,
          size_t dimFlux,
          EntropyType entropy,
          size_t sizeBlock>
constexpr size_t
    PartialMomentBasisBase<DomainFieldType, domainDim, RangeFieldType, rangeDim, dimFlux, entropy, sizeBlock>::
        num_blocks;

template <class DomainFieldType,
          size_t domainDim,
          class RangeFieldType,
          size_t rangeDim,
          size_t dimFlux,
          EntropyType entropy,
          size_t sizeBlock>
constexpr size_t
    PartialMomentBasisBase<DomainFieldType, domainDim, RangeFieldType, rangeDim, dimFlux, entropy, sizeBlock>::
        block_size;


template <class DomainFieldType,
          size_t dimDomain,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols = 1,
          size_t dimFlux = dimDomain,
          size_t order = 1,
          EntropyType entropy = EntropyType::MaxwellBoltzmann>
class PartialMomentBasis
{
  //  static_assert(false, "Not implemented for this combination of dimension and order!");
};

template <class DomainFieldType,
          class RangeFieldType,
          size_t rangeDim,
          size_t rangeDimCols,
          size_t dimFlux,
          EntropyType entropy>
class PartialMomentBasis<DomainFieldType, 1, RangeFieldType, rangeDim, rangeDimCols, dimFlux, 1, entropy>
  : public PartialMomentBasisBase<DomainFieldType, 1, RangeFieldType, rangeDim, dimFlux, entropy, 2>
{
  using BaseType = PartialMomentBasisBase<DomainFieldType, 1, RangeFieldType, rangeDim, dimFlux, entropy, 2>;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  static_assert(!(dimRange % 2), "dimRange has to be even!");
  static constexpr size_t num_intervals = BaseType::num_blocks;

  using typename BaseType::DomainType;
  using typename BaseType::DynamicRangeType;
  using typename BaseType::MatrixType;
  using typename BaseType::QuadraturesType;
  using typename BaseType::RangeType;
  using typename BaseType::SphericalTriangulationType;
  using typename BaseType::StringifierType;
  using LocalVectorType = FieldVector<RangeFieldType, 2>;
  using PartitioningType = typename BaseType::Partitioning1dType;

  static size_t default_quad_order()
  {
    return 15;
  }

  using BaseType::default_quad_refinements;

  PartialMomentBasis(const QuadraturesType& quadratures)
    : BaseType(quadratures)
    , partitioning_(BaseType::create_1d_partitioning(num_intervals))
  {
    BaseType::initialize_base_values();
  }

  PartialMomentBasis(const SphericalTriangulationType& /*triangulation*/, const QuadraturesType& quadratures)
    : PartialMomentBasis(quadratures)
  {}

  PartialMomentBasis(const size_t quad_order = default_quad_order(),
                     const size_t quad_refinements = default_quad_refinements())
    : BaseType(BaseType::gauss_lobatto_quadratures(num_intervals, quad_order, quad_refinements))
    , partitioning_(BaseType::create_1d_partitioning(num_intervals))
  {
    BaseType::initialize_base_values();
  }

  using BaseType::evaluate;

  DynamicRangeType evaluate(const DomainType& v) const override final
  {
    for (size_t ii = 0; ii < num_intervals; ++ii)
      if (XT::Common::FloatCmp::ge(v[0], partitioning_[ii]) && XT::Common::FloatCmp::le(v[0], partitioning_[ii + 1]))
        return evaluate(v, ii);
    DUNE_THROW(Dune::MathError, "v has to be between -1 and 1!");
    return DynamicRangeType();
  } // ... evaluate(...)

  bool adjust_alpha_to_ensure_min_density(RangeType& alpha, const RangeFieldType psi_min) const override final
  {
    bool changed = false;
    const auto alpha_min = std::log(psi_min);
    for (size_t ii = 0; ii < num_intervals; ++ii) {
      auto& alpha_0 = alpha[2 * ii];
      auto& alpha_1 = alpha[2 * ii + 1];
      const auto mu_i = partitioning_[ii];
      const auto mu_ip1 = partitioning_[ii + 1];
      const auto alpha_left = alpha_0 + mu_i * alpha_1;
      const auto alpha_right = alpha_0 + mu_ip1 * alpha_1;
      const auto max_alpha = std::max(alpha_left, alpha_right);
      if (max_alpha < alpha_min) {
        alpha_0 = alpha_min;
        alpha_1 = 0.;
        changed = true;
      }
#if 0
      const bool min_is_left = alpha_left < alpha_right;
      const auto alpha_min_ii = min_is_left ? alpha_left : alpha_right;
      const auto alpha_max_ii = min_is_left ? alpha_right : alpha_left;
      if (alpha_min_ii < alpha_min) {
        if (XT::Common::FloatCmp::le(alpha_max_ii, alpha_min)) {
          alpha_0 = alpha_min;
          alpha_1 = 0.;
          changed = true;
        } else {
          // We know that alpha_1 != 0 because alpha_max_ii > alpha_min and alpha_min_ii < alpha_min
          const auto rho_ii = -(std::exp(alpha_0 + mu_i * alpha_1) - std::exp(alpha_0 + mu_ip1 * alpha_1)) / alpha_1;
          if (std::isnan(rho_ii) || std::isinf(rho_ii))
            DUNE_THROW(Dune::MathError, "Inf or nan in rho!");
          const auto h = (mu_ip1 - mu_i);
          if (XT::Common::FloatCmp::le(rho_ii, psi_min * h)) {
            alpha_0 = alpha_min;
            alpha_1 = 0.;
            changed = true;
            continue;
          }
          // get positive slope
#  if 0
        // Set minimum to alpha_min, leave max alpha unchanged
        alpha_1 = (alpha_max_ii - alpha_min_ii) / h;
#  else
          // Set minimum to alpha_min, leave rho unchanged
          alpha_1 = -1. / h * boost::math::lambert_wm1(-h * psi_min / rho_ii * std::exp(-h * psi_min / rho_ii))
                    - psi_min / rho_ii;
#  endif
          if (!min_is_left)
            alpha_1 *= -1.; // slope has to be negative in this case
          alpha_0 = alpha_min - (min_is_left ? mu_i : mu_ip1) * alpha_1;
          changed = true;
        }
      }
#endif
    } // ii
    return changed;
  }

  // returns matrix with entries <h_i h_j>
  MatrixType mass_matrix() const override final
  {
    MatrixType M(dimRange, dimRange, 0.);
    for (size_t ii = 0; ii < num_intervals; ++ii) {
      M[2 * ii][2 * ii] = partitioning_[ii + 1] - partitioning_[ii];
      M[2 * ii + 1][2 * ii + 1] = (std::pow(partitioning_[ii + 1], 3) - std::pow(partitioning_[ii], 3)) / 3.;
      M[2 * ii][2 * ii + 1] = (std::pow(partitioning_[ii + 1], 2) - std::pow(partitioning_[ii], 2)) / 2.;
      M[2 * ii + 1][2 * ii] = M[2 * ii][2 * ii + 1];
    }
    return M;
  }

  MatrixType mass_matrix_inverse() const override final
  {
    return tridiagonal_matrix_inverse<RangeFieldType, dimRange>(mass_matrix());
  }

  // returns matrix with entries <v h_i h_j>
  FieldVector<MatrixType, dimDomain> flux_matrix() const override final
  {
    MatrixType B(dimRange, dimRange, 0.);
    for (size_t ii = 0; ii < num_intervals; ++ii) {
      B[2 * ii][2 * ii] = (std::pow(partitioning_[ii + 1], 2) - std::pow(partitioning_[ii], 2)) / 2.;
      B[2 * ii + 1][2 * ii + 1] = (std::pow(partitioning_[ii + 1], 4) - std::pow(partitioning_[ii], 4)) / 4.;
      B[2 * ii][2 * ii + 1] = (std::pow(partitioning_[ii + 1], 3) - std::pow(partitioning_[ii], 3)) / 3.;
      B[2 * ii + 1][2 * ii] = B[2 * ii][2 * ii + 1];
    }
    return FieldVector<MatrixType, dimDomain>(B);
  }

  // returns V M^-1 where the matrix V has entries <v h_i h_j>_- and <v h_i h_j>_+
  FieldVector<FieldVector<MatrixType, 2>, 1> kinetic_flux_matrices() const override final
  {
    FieldVector<FieldVector<MatrixType, 2>, 1> ret(FieldVector<MatrixType, 2>(MatrixType(dimRange, dimRange, 0.)));
    auto mm_with_v = flux_matrix();
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
        ret_neg[num_intervals - 1][num_intervals - 1] = -std::pow(partitioning_[num_intervals / 2], 2) / 2;
        ret_pos[num_intervals - 1][num_intervals - 1] = std::pow(partitioning_[num_intervals / 2], 2) / 2;
        // integral corresponding to v basis function
        ret_neg[num_intervals][num_intervals] = -std::pow(partitioning_[num_intervals / 2], 4) / 4;
        ret_pos[num_intervals][num_intervals] = std::pow(partitioning_[num_intervals / 2], 4) / 4;
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
    // apply M^{-1} from the right
    const auto M = std::make_unique<XT::Common::FieldMatrix<RangeFieldType, dimRange, dimRange>>(mass_matrix());
    MatrixType tmp_mat = ret_neg;
    for (size_t rr = 0; rr < dimRange; ++rr)
      M->solve(ret_neg[rr], tmp_mat[rr]);
    tmp_mat = ret_pos;
    for (size_t rr = 0; rr < dimRange; ++rr)
      M->solve(ret_pos[rr], tmp_mat[rr]);
    return ret;
  }

  MatrixType reflection_matrix(const DomainType& n) const override final
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

  const PartitioningType& partitioning() const
  {
    return partitioning_;
  }

  // get indices of all faces that contain point
  std::vector<size_t> get_face_indices(const DomainType& v) const
  {
    std::vector<size_t> face_indices;
    for (size_t jj = 0; jj < partitioning_.size() - 1; ++jj)
      if (XT::Common::FloatCmp::ge(v[0], partitioning_[jj]) && XT::Common::FloatCmp::le(v[0], partitioning_[jj + 1]))
        face_indices.push_back(jj);
    assert(face_indices.size());
    return face_indices;
  }

  DynamicRangeType integrated_initializer(const QuadraturesType& /*quadratures*/) const override final
  {
    DynamicRangeType ret(dimRange, 0);
    for (size_t ii = 0; ii < num_intervals; ++ii) {
      ret[2 * ii] = partitioning_[ii + 1] - partitioning_[ii];
      ret[2 * ii + 1] = (std::pow(partitioning_[ii + 1], 2) - std::pow(partitioning_[ii], 2)) / 2.;
    }
    return ret;
  }

  std::array<int, dimDomain> has_fixed_sign(const size_t index) const override final
  {
    return BaseType::interval_has_fixed_sign(index, num_intervals);
  }

private:
  const PartitioningType partitioning_;
  using BaseType::quadratures_;
}; // class PartialMomentBasis<DomainFieldType, 1, ...>

template <class DomainFieldType,
          class RangeFieldType,
          size_t rangeDim,
          size_t rangeDimCols,
          size_t dimFlux,
          EntropyType entropy>
constexpr size_t
    PartialMomentBasis<DomainFieldType, 1, RangeFieldType, rangeDim, rangeDimCols, dimFlux, 1, entropy>::num_intervals;

template <class DomainFieldType, class RangeFieldType, size_t refinements, size_t dimFlux, EntropyType entropy>
class PartialMomentBasis<DomainFieldType, 3, RangeFieldType, refinements, 1, dimFlux, 1, entropy>
  : public PartialMomentBasisBase<DomainFieldType,
                                  3,
                                  RangeFieldType,
                                  OctaederStatistics<refinements>::num_faces() * 4,
                                  dimFlux,
                                  entropy,
                                  4>
{
  using BaseType = PartialMomentBasisBase<DomainFieldType,
                                          3,
                                          RangeFieldType,
                                          OctaederStatistics<refinements>::num_faces() * 4,
                                          dimFlux,
                                          entropy,
                                          4>;
  using ThisType = PartialMomentBasis;

public:
  using BaseType::block_size;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::num_blocks;
  static constexpr size_t num_refinements = refinements;

  using TriangulationType = typename BaseType::SphericalTriangulationType;
  using typename BaseType::DomainType;
  using typename BaseType::DynamicRangeType;
  using typename BaseType::MatrixType;
  using typename BaseType::QuadraturesType;
  using typename BaseType::QuadratureType;
  using typename BaseType::RangeType;
  using typename BaseType::StringifierType;
  using BlockRangeType = XT::Common::FieldVector<RangeFieldType, block_size>;
  using BlockPlaneCoefficientsType = typename std::vector<std::pair<BlockRangeType, RangeFieldType>>;
  using PlaneCoefficientsType = XT::Common::FieldVector<BlockPlaneCoefficientsType, num_blocks>;
  using BlockMatrixType = XT::Common::BlockedFieldMatrix<RangeFieldType, num_blocks, block_size, block_size>;
  using LocalMatrixType = typename BlockMatrixType::BlockType;
  using LocalVectorType = XT::Common::FieldVector<RangeFieldType, block_size>;

  using BaseType::barycentre_rule;

  static size_t default_quad_order()
  {
    return refinements == 0 ? 15 : 9;
  }

  using BaseType::default_quad_refinements;

  PartialMomentBasis(const QuadraturesType& quadratures)
    : BaseType(refinements, quadratures)
  {
    assert(4 * triangulation_.faces().size() == dimRange);
    BaseType::initialize_base_values();
  }

  PartialMomentBasis(const size_t quad_refinements, const QuadratureRule<RangeFieldType, 2>& reference_quadrature_rule)
    : BaseType(refinements)
  {
    quadratures_ = triangulation_.quadrature_rules(quad_refinements, reference_quadrature_rule);
    assert(4 * triangulation_.faces().size() == dimRange);
    BaseType::initialize_base_values();
  }

  PartialMomentBasis(const TriangulationType& triangulation, const QuadraturesType& quadratures)
    : BaseType(triangulation, quadratures)
  {
    assert(4 * triangulation_.faces().size() == dimRange);
    BaseType::initialize_base_values();
  }

  PartialMomentBasis(const size_t quad_order = default_quad_order(),
                     const size_t quad_refinements = default_quad_refinements())
    : BaseType(refinements)
  {
    const QuadratureRule<RangeFieldType, 2> reference_quadrature_rule =
        XT::Data::FeketeQuadrature<DomainFieldType>::get(quad_order);
    quadratures_ = triangulation_.quadrature_rules(quad_refinements, reference_quadrature_rule);
    assert(4 * triangulation_.faces().size() == dimRange);
    BaseType::initialize_base_values();
  }

  using BaseType::evaluate;

  DynamicRangeType evaluate(const DomainType& v) const override final
  {
    DynamicRangeType ret(dimRange, 0);
    const auto face_indices = triangulation_.get_face_indices(v);
    if (face_indices.size() != 1)
      DUNE_THROW(Dune::MathError,
                 "Either v is not on the unit sphere or on a boundary between triangles where pointwise evaluation is "
                 "not defined for this basis!");
    return evaluate(v, face_indices[0]);
  } // ... evaluate(...)

  RangeFieldType unit_ball_volume() const override final
  {
    return BaseType::unit_ball_volume_quad();
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
  DynamicRangeType integrate_dirac_at(const DomainType& dirac_position) const
  {
    size_t num_faces;
    auto ret = evaluate(dirac_position, true, num_faces);
    if (dirac_position == DomainType{1, 0, 0})
      DXT_ASSERT(num_faces == 4);
    return ret;
  }

  const PlaneCoefficientsType& plane_coefficients() const
  {
    return plane_coefficients_;
  }

  // calculate half space representation of realizable set
  void calculate_plane_coefficients() const
  {
    XT::Common::FieldVector<std::vector<XT::Common::FieldVector<RangeFieldType, block_size>>, num_blocks> points;
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      points[jj].resize(quadratures_[jj].size() + 1);
      for (size_t ll = 0; ll < quadratures_[jj].size(); ++ll) {
        const auto val = evaluate(quadratures_[jj][ll].position(), jj);
        for (size_t ii = 0; ii < block_size; ++ii)
          points[jj][ll][ii] = val[block_size * jj + ii];
      } // ll
      points[jj][quadratures_[jj].size()] = XT::Common::FieldVector<RangeFieldType, block_size>(0.);
    }
    std::vector<std::thread> threads(num_blocks);
    // Calculate plane coefficients for each block in a separate thread
    for (size_t jj = 0; jj < num_blocks; ++jj)
      threads[jj] = std::thread(&ThisType::calculate_plane_coefficients_block, this, std::ref(points[jj]), jj);
    // Join the threads with the main thread
    for (size_t jj = 0; jj < num_blocks; ++jj)
      threads[jj].join();
  }

  std::unique_ptr<BlockMatrixType> block_mass_matrix() const
  {
    auto block_matrix = std::make_unique<BlockMatrixType>();
    parallel_quadrature_blocked(quadratures_, *block_matrix, size_t(-1));
    return block_matrix;
  } // ... mass_matrix()

  MatrixType mass_matrix() const override final
  {
    return block_mass_matrix()->convert_to_dynamic_matrix();
  } // ... mass_matrix()

  FieldVector<MatrixType, dimFlux> flux_matrix() const override final
  {
    FieldVector<MatrixType, dimFlux> B(MatrixType(dimRange, dimRange, 0));
    BlockMatrixType block_matrix;
    for (size_t dd = 0; dd < dimFlux; ++dd) {
      parallel_quadrature_blocked(quadratures_, block_matrix, dd);
      B[dd] = block_matrix.convert_to_dynamic_matrix();
    }
    return B;
  }

  // returns V M^-1 where V has entries <v h_i h_j>_- and <v h_i h_j>_+
  FieldVector<FieldVector<MatrixType, 2>, dimFlux> kinetic_flux_matrices() const override final
  {
    FieldVector<FieldVector<MatrixType, 2>, dimFlux> B_kinetic(
        FieldVector<MatrixType, 2>(MatrixType(dimRange, dimRange, 0.)));
    BlockMatrixType block_matrix;
    const auto mass_matrix = block_mass_matrix();
    for (size_t dd = 0; dd < dimFlux; ++dd) {
      QuadraturesType neg_quadratures(quadratures_.size());
      QuadraturesType pos_quadratures(quadratures_.size());
      BaseType::get_pos_and_neg_quadratures(neg_quadratures, pos_quadratures, dd);
      parallel_quadrature_blocked(neg_quadratures, block_matrix, dd);
      apply_invM_from_right(*mass_matrix, block_matrix, B_kinetic[dd][0]);
      parallel_quadrature_blocked(pos_quadratures, block_matrix, dd);
      apply_invM_from_right(*mass_matrix, block_matrix, B_kinetic[dd][1]);
    } // dd
    return B_kinetic;
  } // ... kinetic_flux_matrices()

  std::array<int, dimDomain> has_fixed_sign(const size_t index) const override final
  {
    return BaseType::triangle_has_fixed_sign(index);
  }

private:
  void calculate_plane_coefficients_block(std::vector<XT::Common::FieldVector<RangeFieldType, block_size>>& points,
                                          const size_t jj) const
  {
#if HAVE_QHULL
    orgQhull::Qhull qhull;
    // ignore output
    boost::iostreams::stream<boost::iostreams::null_sink> null_ostream((boost::iostreams::null_sink()));
    qhull.setOutputStream(&null_ostream);
    qhull.setErrorStream(&null_ostream);
    // calculate convex hull
    assert(points.size() < std::numeric_limits<int>::max());
    qhull.runQhull(
        "Realizable set", static_cast<int>(block_size), static_cast<int>(points.size()), &(points[0][0]), "Qt T1");
    const auto facet_end = qhull.endFacet();
    BlockPlaneCoefficientsType block_plane_coefficients(qhull.facetList().count());
    //    std::cout << "num_vertices: " << qhull.vertexList().count() << std::endl;
    size_t ll = 0;
    for (auto facet = qhull.beginFacet(); facet != facet_end; facet = facet.next(), ++ll) {
      for (size_t ii = 0; ii < block_size; ++ii)
        block_plane_coefficients[ll].first[ii] = *(facet.hyperplane().coordinates() + ii);
      block_plane_coefficients[ll].second = -facet.hyperplane().offset();
    } // ii
    // discard duplicate facets (qhull triangulates output, so there may be several facets on the same hyperplane)
    using CoeffType = typename BlockPlaneCoefficientsType::value_type;
    std::sort(block_plane_coefficients.begin(),
              block_plane_coefficients.end(),
              [](const CoeffType& first, const CoeffType& second) {
                // Check component-wise if first.a[ii] < second.a[ii]. If they are equal, check next component.
                for (size_t ii = 0; ii < block_size; ++ii) {
                  if (XT::Common::FloatCmp::lt(first.first[ii], second.first[ii]))
                    return true;
                  else if (XT::Common::FloatCmp::gt(first.first[ii], second.first[ii]))
                    return false;
                }
                // first.a and second.a are equal, check first.b and second.b
                if (XT::Common::FloatCmp::lt(first.second, second.second))
                  return true;
                return false;
              });
    static const auto pair_float_cmp = [](const CoeffType& first, const CoeffType& second) {
      return XT::Common::FloatCmp::eq(first.first, second.first)
             && XT::Common::FloatCmp::eq(first.second, second.second);
    };
    block_plane_coefficients.erase(
        std::unique(block_plane_coefficients.begin(), block_plane_coefficients.end(), pair_float_cmp),
        block_plane_coefficients.end());
    // discard facets belonging to the condition u0 < 1, i.e. with a = (1, 0, 0, ...) and b = 1
    CoeffType coeff_to_erase{BlockRangeType(0), 1.};
    coeff_to_erase.first[0] = 1.;
    auto coeff_to_erase_it =
        std::find_if(block_plane_coefficients.begin(),
                     block_plane_coefficients.end(),
                     [&coeff_to_erase](const CoeffType& coeff) { return pair_float_cmp(coeff, coeff_to_erase); });
    if (coeff_to_erase_it == block_plane_coefficients.end())
      DUNE_THROW(Dune::MathError, "There should be such a coefficient!");
    block_plane_coefficients.erase(coeff_to_erase_it);
    plane_coefficients_[jj] = block_plane_coefficients;
#else // HAVE_QHULL
    DUNE_UNUSED_PARAMETER(points);
    DUNE_UNUSED_PARAMETER(jj);
    DUNE_THROW(Dune::NotImplemented, "You are missing Qhull!");
#endif // HAVE_QHULL
  }

  void
  parallel_quadrature_blocked(const QuadraturesType& quadratures, BlockMatrixType& matrix, const size_t v_index) const
  {
    const size_t num_threads = std::min(XT::Common::threadManager().max_threads(), num_blocks);
    std::vector<std::thread> threads(num_threads);
    // Launch a group of threads
    size_t blocks_done = 0;
    while (blocks_done < num_blocks) {
      size_t inner_loop_count = std::min(num_threads, num_blocks - blocks_done);
      for (size_t jj = 0; jj < inner_loop_count; ++jj)
        threads[jj] = std::thread(&ThisType::calculate_block_in_thread,
                                  this,
                                  std::cref(quadratures[blocks_done + jj]),
                                  std::ref(matrix.block(blocks_done + jj)),
                                  v_index,
                                  blocks_done + jj);
      // Join the threads with the main thread
      for (size_t jj = 0; jj < inner_loop_count; ++jj)
        threads[jj].join();
      blocks_done += inner_loop_count;
    }
  } // void parallel_quadrature(...)

  void calculate_block_in_thread(const QuadratureType& quadrature,
                                 LocalMatrixType& local_matrix,
                                 const size_t v_index,
                                 const size_t jj) const
  {
    local_matrix *= 0.;
    for (const auto& quad_point : quadrature) {
      const auto& v = quad_point.position();
      const auto basis_evaluated = evaluate(v, jj);
      const auto& weight = quad_point.weight();
      const auto factor = v_index == size_t(-1) ? 1. : v[v_index];
      for (size_t nn = 0; nn < local_matrix.N(); ++nn)
        for (size_t mm = 0; mm < local_matrix.M(); ++mm)
          local_matrix[nn][mm] +=
              basis_evaluated[jj * block_size + nn] * basis_evaluated[jj * block_size + mm] * factor * weight;
    }
  } // void calculate_in_thread(...)

  void apply_invM_from_right(const BlockMatrixType& M, BlockMatrixType& V, MatrixType& ret) const
  {
    LocalVectorType local_vec;
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      for (size_t rr = 0; rr < block_size; ++rr) {
        M.block(jj).solve(local_vec, V.block(jj)[rr]);
        V.block(jj)[rr] = local_vec;
      } // rr
    } // jj
    ret = V.convert_to_dynamic_matrix();
  }

  using BaseType::quadratures_;
  using BaseType::triangulation_;
  mutable PlaneCoefficientsType plane_coefficients_;
}; // class PartialMomentBasis<DomainFieldType, 3, ...>

template <class DomainFieldType, class RangeFieldType, size_t refinements, size_t dimFlux, EntropyType entropy>
constexpr size_t
    PartialMomentBasis<DomainFieldType, 3, RangeFieldType, refinements, 1, dimFlux, 1, entropy>::num_refinements;


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_PARTIALMOMENTS_HH
