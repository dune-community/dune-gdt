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

#include <string>
#include <vector>

#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/common/numeric.hh>

#include <dune/gdt/test/momentmodels/triangulation.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {


template <class DomainFieldType,
          size_t domainDim,
          class RangeFieldType,
          size_t rangeDim,
          size_t fluxDim,
          EntropyType entropy>
class HatFunctionMomentBasisBase
  : public MomentBasisInterface<DomainFieldType, domainDim, RangeFieldType, rangeDim, 1, fluxDim, entropy>
{
  using BaseType = MomentBasisInterface<DomainFieldType, domainDim, RangeFieldType, rangeDim, 1, fluxDim, entropy>;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;

  using typename BaseType::DynamicRangeType;
  using typename BaseType::RangeType;
  using typename BaseType::StringifierType;

  template <class... Args>
  HatFunctionMomentBasisBase(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}

  static StringifierType stringifier()
  {
    return [](const RangeType& val) {
      RangeFieldType psi(0);
      for (const auto& entry : val)
        psi += entry;
      return XT::Common::to_string(psi, 15);
    };
  } // ... stringifier()

  DynamicRangeType alpha_one() const override final
  {
    return DynamicRangeType(dimRange, 1.);
  }

  RangeFieldType density(const DynamicRangeType& u) const override final
  {
    return XT::Common::reduce(&(u[0]), &(u[0]) + dimRange, RangeFieldType(0));
  }

  RangeFieldType density(const RangeType& u) const override final
  {
    return XT::Common::reduce(&(u[0]), &(u[0]) + dimRange, RangeFieldType(0));
  }

  template <class Vec>
  std::enable_if_t<XT::Common::is_vector<Vec>::value, void> alpha_one(Vec& ret) const
  {
    for (size_t ii = 0; ii < ret.size(); ++ii)
      XT::Common::VectorAbstraction<Vec>::set_entry(ret, ii, 1.);
  }

  template <class Vec>
  std::enable_if_t<XT::Common::is_vector<Vec>::value && !std::is_same<Vec, DynamicRangeType>::value, RangeFieldType>
  density(const Vec& u) const
  {
    RangeFieldType ret(0.);
    for (size_t ii = 0; ii < u.size(); ++ii)
      ret += XT::Common::VectorAbstraction<Vec>::get_entry(u, ii);
    return ret;
  }

  using BaseType::u_iso;

  template <class Vec>
  std::enable_if_t<XT::Common::is_vector<Vec>::value, void> u_iso(Vec& ret) const
  {
    const auto& ret_range = u_iso();
    using V = XT::Common::VectorAbstraction<Vec>;
    for (size_t ii = 0; ii < dimRange; ++ii)
      V::set_entry(ret, ii, ret_range[ii]);
  }

  void ensure_min_density(DynamicRangeType& u, const RangeFieldType psi_min) const override final
  {
    const auto u_iso_min = u_iso() * psi_min;
    for (size_t ii = 0; ii < dimRange; ++ii)
      if (u[ii] < u_iso_min[ii])
        u[ii] = u_iso_min[ii];
  }

  void ensure_min_density(RangeType& u, const RangeFieldType psi_min) const override final
  {
    const auto u_iso_min = u_iso() * psi_min;
    for (size_t ii = 0; ii < dimRange; ++ii)
      if (u[ii] < u_iso_min[ii])
        u[ii] = u_iso_min[ii];
  }

  bool adjust_alpha_to_ensure_min_density(RangeType& alpha, const RangeFieldType psi_min) const override final
  {
    bool changed = false;
    const auto alpha_min = std::log(psi_min);
    for (size_t ii = 0; ii < dimRange; ++ii) {
      if (alpha[ii] < alpha_min) {
        alpha[ii] = alpha_min;
        changed = true;
      }
    }
    return changed;
  }

  std::string short_id() const override final
  {
    return "hf";
  }

  std::string mn_name() const override final
  {
    return "hfm" + XT::Common::to_string(dimRange);
  }

  std::string pn_name() const override final
  {
    return "hfp" + XT::Common::to_string(dimRange);
  }
}; // class HatFunctionMomentBasisBase


template <class DomainFieldType,
          size_t dimDomain,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols = 1,
          size_t dimFlux = dimDomain,
          EntropyType entropy = EntropyType::MaxwellBoltzmann>
class HatFunctionMomentBasis
{
  //  static_assert(false, "Not implemented for this dimension!");
};

template <class DomainFieldType,
          class RangeFieldType,
          size_t rangeDim,
          size_t rangeDimCols,
          size_t fluxDim,
          EntropyType entropy>
class HatFunctionMomentBasis<DomainFieldType, 1, RangeFieldType, rangeDim, rangeDimCols, fluxDim, entropy>
  : public HatFunctionMomentBasisBase<DomainFieldType, 1, RangeFieldType, rangeDim, fluxDim, entropy>
{
public:
  static const size_t dimDomain = 1;
  static const size_t dimRange = rangeDim;
  static const size_t dimRangeCols = rangeDimCols;
  static const size_t num_intervals = dimRange - 1;

private:
  using BaseType = HatFunctionMomentBasisBase<DomainFieldType, 1, RangeFieldType, rangeDim, fluxDim, entropy>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::DynamicRangeType;
  using typename BaseType::MatrixType;
  using typename BaseType::QuadraturesType;
  using typename BaseType::RangeType;
  using typename BaseType::SphericalTriangulationType;
  using typename BaseType::StringifierType;
  using PartitioningType = typename BaseType::Partitioning1dType;

  static size_t default_quad_order()
  {
    return 15;
  }

  using BaseType::default_quad_refinements;

  HatFunctionMomentBasis(const QuadraturesType& quadratures)
    : BaseType(quadratures)
    , partitioning_(BaseType::create_1d_partitioning(num_intervals))
  {
    BaseType::initialize_base_values();
  }

  HatFunctionMomentBasis(const SphericalTriangulationType& /*triangulation*/, const QuadraturesType& quadratures)
    : HatFunctionMomentBasis(quadratures)
  {}

  HatFunctionMomentBasis(const size_t quad_order = default_quad_order(),
                         const size_t quad_refinements = default_quad_refinements())
    : BaseType(BaseType::gauss_lobatto_quadratures(num_intervals, quad_order, quad_refinements))
    , partitioning_(BaseType::create_1d_partitioning(num_intervals))
  {
    BaseType::initialize_base_values();
  }

  DynamicRangeType evaluate(const DomainType& v) const override final
  {
    DynamicRangeType ret(dimRange, 0);
    const auto& mu = partitioning_;
    for (size_t ii = 0; ii < dimRange; ++ii) {
      if (ii < num_intervals && XT::Common::FloatCmp::ge(v[0], mu[ii]) && XT::Common::FloatCmp::le(v[0], mu[ii + 1]))
        ret[ii] = (v - mu[ii + 1]) / (mu[ii] - mu[ii + 1]);
      if (ii > 0 && XT::Common::FloatCmp::ge(v[0], mu[ii - 1]) && XT::Common::FloatCmp::le(v[0], mu[ii]))
        ret[ii] = (v - mu[ii - 1]) / (mu[ii] - mu[ii - 1]);
    }
    return ret;
  } // ... evaluate(...)

  DynamicRangeType evaluate(const DomainType& v, const size_t interval_index) const override final
  {
    DynamicRangeType ret(dimRange, 0);
    const auto& mu = partitioning_;
    ret[interval_index] = (v - mu[interval_index + 1]) / (mu[interval_index] - mu[interval_index + 1]);
    ret[interval_index + 1] = (v - mu[interval_index]) / (mu[interval_index + 1] - mu[interval_index]);
    return ret;
  } // ... evaluate(...)

  XT::Common::FieldVector<RangeFieldType, 2> evaluate_on_interval(const DomainType& v,
                                                                  const size_t interval_index) const
  {
    XT::Common::FieldVector<RangeFieldType, 2> ret;
    const auto& mu = partitioning_;
    ret[0] = (v - mu[interval_index + 1]) / (mu[interval_index] - mu[interval_index + 1]);
    ret[1] = (v - mu[interval_index]) / (mu[interval_index + 1] - mu[interval_index]);
    return ret;
  } // ... evaluate(...)

  // returns matrix with entries <h_i h_j>
  MatrixType mass_matrix() const override final
  {
    MatrixType ret(dimRange, dimRange, 0);
    const auto& mu = partitioning_;
    ret[0][0] = (mu[1] - mu[0]) / 3.;
    for (size_t rr = 0; rr < dimRange; ++rr) {
      if (rr > 0 && rr < num_intervals)
        ret[rr][rr] = (mu[rr + 1] - mu[rr - 1]) / 3.;
      if (rr > 0)
        ret[rr][rr - 1] = (mu[rr] - mu[rr - 1]) / 6.;
      if (rr < num_intervals)
        ret[rr][rr + 1] = (mu[rr + 1] - mu[rr]) / 6.;
    }
    ret[num_intervals][num_intervals] = (mu[num_intervals] - mu[num_intervals - 1]) / 3.;
    return ret;
  }

  MatrixType mass_matrix_inverse() const override final
  {
    return tridiagonal_matrix_inverse<RangeFieldType, dimRange>(mass_matrix());
  }

  // returns matrix with entries <v h_i h_j>
  FieldVector<MatrixType, 1> flux_matrix() const override final
  {
    MatrixType ret(dimRange, dimRange, 0.);
    const auto& mu = partitioning_;
    ret[0][0] = (std::pow(mu[1], 2) + 2 * mu[1] * mu[0] - 3 * std::pow(mu[0], 2)) / 12.;
    for (size_t rr = 0; rr < dimRange; ++rr) {
      if (rr > 0 && rr < num_intervals)
        ret[rr][rr] =
            (std::pow(mu[rr + 1], 2) + 2 * mu[rr + 1] * mu[rr] - 2 * mu[rr] * mu[rr - 1] - std::pow(mu[rr - 1], 2))
            / 12.;
      if (rr > 0)
        ret[rr][rr - 1] = (std::pow(mu[rr], 2) - std::pow(mu[rr - 1], 2)) / 12.;
      if (rr < num_intervals)
        ret[rr][rr + 1] = (std::pow(mu[rr + 1], 2) - std::pow(mu[rr], 2)) / 12.;
    }
    ret[num_intervals][num_intervals] =
        (3 * std::pow(mu[num_intervals], 2) - 2 * mu[num_intervals] * mu[num_intervals - 1]
         - std::pow(mu[num_intervals - 1], 2))
        / 12.;
    return ret;
  }

  // returns V M^-1 where the matrix V has entries <v h_i h_j>_- and <v h_i h_j>_+
  FieldVector<FieldVector<MatrixType, 2>, 1> kinetic_flux_matrices() const override final
  {
    FieldVector<FieldVector<MatrixType, 2>, 1> ret(FieldVector<MatrixType, 2>(MatrixType(dimRange, dimRange, 0.)));
    auto mm_with_v = flux_matrix();
    auto& ret_neg = ret[0][0];
    auto& ret_pos = ret[0][1];
    size_t N = dimRange;
    const auto& mu = partitioning_;
    for (size_t nn = 0; nn < N; ++nn) {
      for (size_t mm = (nn > 0 ? nn - 1 : 0); mm <= (nn < N - 1 ? nn + 1 : nn); ++mm) {
        if (N % 2) {
          if (nn < N / 2 || mm < N / 2)
            ret_neg[nn][mm] = mm_with_v[0][nn][mm];
          else if (nn > N / 2 || mm > N / 2)
            ret_pos[nn][mm] = mm_with_v[0][nn][mm];
          else { // nn == mm == N/2
            ret_neg[nn][mm] = -std::pow(mu[mm - 1], 2) / 12.;
            ret_pos[nn][mm] = std::pow(mu[mm + 1], 2) / 12.;
          }
        } else {
          if (nn < N / 2 - 1 || mm < N / 2 - 1)
            ret_neg[nn][mm] = mm_with_v[0][nn][mm];
          else if (nn > N / 2 || mm > N / 2)
            ret_pos[nn][mm] = mm_with_v[0][nn][mm];
          else if (nn == N / 2 && mm == nn) {
            ret_neg[nn][mm] = -std::pow(mu[mm], 2) / 48.;
            ret_pos[nn][mm] = (5 * std::pow(mu[mm], 2) + 8 * mu[mm] * mu[mm + 1] + 4 * std::pow(mu[mm + 1], 2)) / 48.;
          } else if (nn == N / 2 - 1 && mm == nn) {
            ret_neg[nn][mm] = (-5 * std::pow(mu[mm], 2) - 8 * mu[mm] * mu[mm - 1] - 4 * std::pow(mu[mm - 1], 2)) / 48.;
            ret_pos[nn][mm] = std::pow(mu[mm], 2) / 48.;
          } else { // (((mm == N / 2 && nn == N / 2 - 1) || (mm == N / 2 - 1 && nn == N / 2))) {
            ret_neg[nn][mm] = -std::pow(mu[mm], 2) / 16.;
            ret_pos[nn][mm] = std::pow(mu[mm], 2) / 16.;
          }
        } // else (N % 2)
      } // mm
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
    for (size_t ii = 0; ii < dimRange; ++ii)
      for (size_t jj = 0; jj < dimRange; ++jj)
        ret[ii][jj] = mass_mat[ii][num_intervals - jj];
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

  // get indices of all faces that contain point v
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
    const auto& mu = partitioning_;
    ret[0] = mu[1] - mu[0];
    for (size_t ii = 1; ii < num_intervals; ++ii)
      ret[ii] = mu[ii + 1] - mu[ii - 1];
    ret[num_intervals] = mu[num_intervals] - mu[num_intervals - 1];
    ret *= 0.5;
    return ret;
  }

  const PartitioningType& partitioning() const
  {
    return partitioning_;
  }

  DynamicRangeType
  get_u(const std::function<RangeFieldType(DomainType, std::array<int, dimDomain>)>& psi) const override final
  {
    DynamicRangeType ret(dimRange, 0.);
    const auto merged_quads = XT::Data::merged_quadrature(quadratures_);
    for (auto it = merged_quads.begin(); it != merged_quads.end(); ++it) {
      const auto& quad_point = *it;
      const size_t interval_index = it.first_index();
      const auto& v = quad_point.position();
      const auto val = evaluate_on_interval(v, interval_index);
      const auto factor = psi(v, has_fixed_sign(interval_index)) * quad_point.weight();
      for (size_t ii = 0; ii < 2; ++ii)
        ret[interval_index + ii] += val[ii] * factor;
    }
    return ret;
  }

  std::array<int, dimDomain> has_fixed_sign(const size_t index) const override final
  {
    return BaseType::interval_has_fixed_sign(index, num_intervals);
  }

private:
  using BaseType::quadratures_;
  const PartitioningType partitioning_;
}; // class HatFunctionMomentBasis<DomainFieldType, 1, ...>

template <class DomainFieldType,
          class RangeFieldType,
          size_t rangeDim,
          size_t rangeDimCols,
          size_t fluxDim,
          EntropyType entropy>
constexpr size_t
    HatFunctionMomentBasis<DomainFieldType, 1, RangeFieldType, rangeDim, rangeDimCols, fluxDim, entropy>::dimRange;

template <class DomainFieldType, class RangeFieldType, size_t refinements, size_t fluxDim, EntropyType entropy>
class HatFunctionMomentBasis<DomainFieldType, 3, RangeFieldType, refinements, 1, fluxDim, entropy>
  : public HatFunctionMomentBasisBase<DomainFieldType,
                                      3,
                                      RangeFieldType,
                                      OctaederStatistics<refinements>::num_vertices(),
                                      fluxDim,
                                      entropy>
{
public:
  static constexpr size_t dimDomain = 3;
  static constexpr size_t dimRange = OctaederStatistics<refinements>::num_vertices();
  static constexpr size_t dimRangeCols = 1;
  static constexpr size_t dimFlux = fluxDim;
  static constexpr size_t num_refinements = refinements;

private:
  using BaseType = HatFunctionMomentBasisBase<DomainFieldType,
                                              3,
                                              RangeFieldType,
                                              OctaederStatistics<refinements>::num_vertices(),
                                              fluxDim,
                                              entropy>;
  using ThisType = HatFunctionMomentBasis;

public:
  using TriangulationType = typename BaseType::SphericalTriangulationType;
  using typename BaseType::DomainType;
  using typename BaseType::DynamicRangeType;
  using typename BaseType::MatrixType;
  using typename BaseType::MergedQuadratureIterator;
  using typename BaseType::QuadraturesType;
  using typename BaseType::RangeType;
  using typename BaseType::StringifierType;
  using LocalMatrixType = XT::Common::FieldMatrix<RangeFieldType, 3, 3>;

  using BaseType::barycentre_rule;

  static size_t default_quad_order()
  {
    return refinements == 0 ? 15 : 9;
  }

  using BaseType::default_quad_refinements;

  HatFunctionMomentBasis(const QuadraturesType& quadratures)
    : BaseType(refinements, quadratures)
  {
    assert(triangulation_.vertices().size() == dimRange);
    BaseType::initialize_base_values();
  }

  HatFunctionMomentBasis(const TriangulationType& triangulation, const QuadraturesType& quadratures)
    : BaseType(triangulation, quadratures)
  {
    assert(triangulation_.vertices().size() == dimRange);
    BaseType::initialize_base_values();
  }

  HatFunctionMomentBasis(const size_t quad_refinements,
                         const QuadratureRule<RangeFieldType, 2>& reference_quadrature_rule)
    : BaseType(refinements)
  {
    quadratures_ = triangulation_.quadrature_rules(quad_refinements, reference_quadrature_rule);
    assert(triangulation_.vertices().size() == dimRange);
    BaseType::initialize_base_values();
  }

  HatFunctionMomentBasis(const size_t quad_order = default_quad_order(),
                         const size_t quad_refinements = default_quad_refinements())
    : BaseType(refinements)
  {
    const QuadratureRule<RangeFieldType, 2> reference_quadrature_rule =
        XT::Data::FeketeQuadrature<DomainFieldType>::get(quad_order);
    quadratures_ = triangulation_.quadrature_rules(quad_refinements, reference_quadrature_rule);
    assert(triangulation_.vertices().size() == dimRange);
    BaseType::initialize_base_values();
  }

  DynamicRangeType evaluate(const DomainType& v) const override
  {
    DynamicRangeType ret(dimRange, 0);
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

  DynamicRangeType evaluate(const DomainType& v, const size_t face_index) const override final
  {
    DynamicRangeType ret(dimRange, 0);
    auto barycentric_coords = evaluate_on_face(v, face_index);
    const auto& vertices = triangulation_.faces()[face_index]->vertices();
    for (size_t ii = 0; ii < 3; ++ii)
      ret[vertices[ii]->index()] = barycentric_coords[ii];
    return ret;
  } // ... evaluate(...)

  DomainType evaluate_on_face(const DomainType& v, const size_t face_index) const
  {
    DomainType ret(0);
    const auto& face = triangulation_.faces()[face_index];
    const auto& vertices = face->vertices();
    bool success = calculate_barycentric_coordinates(v, vertices, ret);
    assert(success);
#ifdef NDEBUG
    static_cast<void>(success);
#endif
    return ret;
  } // ... evaluate(...)

  const TriangulationType& triangulation() const
  {
    return triangulation_;
  }

  RangeFieldType unit_ball_volume() const override final
  {
    return BaseType::unit_ball_volume_quad();
  }

  // get indices of all faces that contain point v
  std::vector<size_t> get_face_indices(const DomainType& v) const
  {
    return triangulation_.get_face_indices(v);
  }

  // calculates <b(v) dirac(v-dirac_position)>
  DynamicRangeType integrate_dirac_at(const DomainType& dirac_position) const
  {
    return evaluate(dirac_position);
  }

  using BaseType::u_iso;

  DynamicRangeType
  get_u(const std::function<RangeFieldType(DomainType, std::array<int, dimDomain>)>& psi) const override final
  {
    DynamicRangeType ret(dimRange, 0.);
    const auto merged_quads = XT::Data::merged_quadrature(quadratures_);
    for (auto it = merged_quads.begin(); it != merged_quads.end(); ++it) {
      const auto face_index = it.first_index();
      const auto& vertices = triangulation_.faces()[face_index]->vertices();
      const auto& quad_point = *it;
      const auto& v = quad_point.position();
      const auto val = evaluate_on_face(v, face_index);
      const auto factor = psi(v, has_fixed_sign(it)) * quad_point.weight();
      for (size_t ii = 0; ii < 3; ++ii)
        ret[vertices[ii]->index()] += val[ii] * factor;
    }
    return ret;
  }

  std::array<int, dimDomain> has_fixed_sign(const size_t index) const override final
  {
    return BaseType::triangle_has_fixed_sign(index);
  }

protected:
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
    if (XT::Common::FloatCmp::lt(ret[0], 0., 1e-14, 1e-14) || XT::Common::FloatCmp::lt(ret[1], 0., 1e-14, 1e-14)
        || XT::Common::FloatCmp::lt(ret[2], 0., 1e-14, 1e-14))
      return false;
    // we already checked the values are close to 0, now remove negative values stemming from numerical inaccuracies
    for (size_t ii = 0; ii < 3; ++ii)
      if (ret[ii] < 0.)
        ret[ii] = 0.;
    return true;
  } // bool calculate_barycentric_coordinates(...)

  static std::vector<size_t> create_face_decomposition(const size_t num_faces, const size_t num_threads)
  {
    std::vector<size_t> ret(num_threads + 1);
    for (size_t ii = 0; ii < num_threads; ++ii)
      ret[ii] = num_faces / num_threads * ii;
    ret[num_threads] = num_faces;
    return ret;
  }

  void parallel_quadrature(const QuadraturesType& quadratures,
                           MatrixType& matrix,
                           const size_t v_index,
                           const bool reflecting = false) const override final
  {
    const auto& faces = triangulation_.faces();
    size_t num_threads = std::min(XT::Common::threadManager().max_threads(), faces.size());
    const auto decomposition = create_face_decomposition(faces.size(), num_threads);
    std::vector<std::thread> threads(num_threads);
    // Launch a group of threads
    std::vector<LocalMatrixType> local_matrices(faces.size(), LocalMatrixType(0.));
    for (size_t ii = 0; ii < num_threads; ++ii)
      threads[ii] = std::thread(&ThisType::calculate_in_thread_hat,
                                this,
                                std::ref(local_matrices),
                                std::cref(quadratures),
                                v_index,
                                std::cref(decomposition),
                                ii,
                                reflecting);
    // Join the threads with the main thread
    for (size_t ii = 0; ii < num_threads; ++ii)
      threads[ii].join();
    // add local matrices
    matrix *= 0.;
    for (size_t ii = 0; ii < num_threads; ++ii) {
      for (size_t face_index = decomposition[ii]; face_index < decomposition[ii + 1]; ++face_index) {
        const auto& face = faces[face_index];
        const auto& vertices = face->vertices();
        for (size_t nn = 0; nn < 3; ++nn)
          for (size_t mm = 0; mm < 3; ++mm)
            matrix[vertices[nn]->index()][vertices[mm]->index()] += local_matrices[face_index][nn][mm];
      } // faces
    } // threads
  } // void parallel_quadrature(...)

  void calculate_in_thread_hat(std::vector<LocalMatrixType>& local_matrices,
                               const QuadraturesType& quadratures,
                               const size_t v_index,
                               const std::vector<size_t>& decomposition,
                               const size_t ii,
                               const bool reflecting) const
  {
    const auto& reflected_indices = triangulation_.reflected_face_indices();
    for (size_t face_index = decomposition[ii]; face_index < decomposition[ii + 1]; ++face_index) {
      for (const auto& quad_point : quadratures[face_index]) {
        const auto& v = quad_point.position();
        const auto basis_evaluated = evaluate_on_face(v, face_index);
        auto basis_reflected = basis_evaluated;
        if (reflecting) {
          auto v_reflected = v;
          v_reflected[v_index] *= -1.;
          const size_t reflected_index = reflected_indices.size() ? reflected_indices[face_index][v_index] : 0;
          basis_reflected = evaluate_on_face(v_reflected, reflected_index);
        }
        const auto& weight = quad_point.weight();
        const auto factor = (reflecting || v_index == size_t(-1)) ? 1. : v[v_index];
        for (size_t nn = 0; nn < 3; ++nn)
          for (size_t mm = 0; mm < 3; ++mm)
            local_matrices[face_index][nn][mm] +=
                basis_evaluated[nn] * (reflecting ? basis_reflected[mm] : basis_evaluated[mm]) * factor * weight;
      } // quad_points
    } // faces
  } // void calculate_in_thread(...)

  void integrated_initializer_thread(DynamicRangeType& local_range,
                                     const std::vector<MergedQuadratureIterator>& decomposition,
                                     const size_t ii) const override final
  {
    const auto& faces = triangulation_.faces();
    for (auto it = decomposition[ii]; it != decomposition[ii + 1]; ++it) {
      const auto face_index = it.first_index();
      const auto& vertices = faces[face_index]->vertices();
      const auto& quad_point = *it;
      DomainType basis_evaluated = evaluate_on_face(quad_point.position(), face_index);
      basis_evaluated *= quad_point.weight();
      for (size_t nn = 0; nn < 3; ++nn)
        local_range[vertices[nn]->index()] += basis_evaluated[nn];
    } // jj
  } // void calculate_in_thread(...)

  using BaseType::quadratures_;
  using BaseType::triangulation_;
}; // class HatFunctionMomentBasis<DomainFieldType, 3, ...>

template <class DomainFieldType, class RangeFieldType, size_t rangeDim, size_t fluxDim, EntropyType entropy>
constexpr size_t HatFunctionMomentBasis<DomainFieldType, 3, RangeFieldType, rangeDim, 1, fluxDim, entropy>::dimRange;

template <class DomainFieldType, class RangeFieldType, size_t rangeDim, size_t fluxDim, EntropyType entropy>
constexpr size_t
    HatFunctionMomentBasis<DomainFieldType, 3, RangeFieldType, rangeDim, 1, fluxDim, entropy>::num_refinements;


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_HATFUNCTIONS_HH
