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

#if HAVE_QHULL
#  include <dune/xt/common/disable_warnings.hh>
#  include <libqhullcpp/Qhull.h>
#  include <libqhullcpp/QhullFacetList.h>
#  include <dune/xt/common/reenable_warnings.hh>
#endif // HAVE_QHULL

#include <dune/xt/common/fvector.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


template <class DomainFieldType,
          size_t dimDomain,
          class RangeFieldType,
          size_t dimRange,
          size_t dimRangeCols = 1,
          size_t dimFlux = dimDomain,
          size_t order = 1>
class PartialMomentBasis
{
  //  static_assert(false, "Not implemented for this combination of dimension and order!");
};

template <class DomainFieldType, class RangeFieldType, size_t rangeDim, size_t rangeDimCols, size_t dimFlux>
class PartialMomentBasis<DomainFieldType, 1, RangeFieldType, rangeDim, rangeDimCols, dimFlux, 1>
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

  PartialMomentBasis(const QuadraturesType& quadratures)
    : BaseType(quadratures)
    , triangulation_(BaseType::create_1d_triangulation(num_intervals))
  {
    BaseType::initialize_base_values();
  }

  PartialMomentBasis(const size_t quad_order = 15, const size_t DXTC_DEBUG_ONLY(quad_refinements) = 0)
    : BaseType(BaseType::gauss_lobatto_quadratures(num_intervals, quad_order))
    , triangulation_(BaseType::create_1d_triangulation(num_intervals))
  {
    assert(quad_refinements == 0 && "Refinement of the quadrature intervals not implemented for this basis!");
    BaseType::initialize_base_values();
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

  virtual RangeType integrated(const bool /*use_fine_quadratures*/ = false) const override final
  {
    RangeType ret(0);
    for (size_t ii = 0; ii < num_intervals; ++ii) {
      ret[2 * ii] = triangulation_[ii + 1] - triangulation_[ii];
      ret[2 * ii + 1] = (std::pow(triangulation_[ii + 1], 2) - std::pow(triangulation_[ii], 2)) / 2.;
    }
    return ret;
  }

  virtual RangeType get_moment_vector(const std::function<RangeFieldType(DomainType, bool)>& psi) const override final
  {
    RangeType ret(0.);
    const auto merged_quadratures = quadratures_.merged();
    for (auto it = merged_quadratures.begin(); it != merged_quadratures.end(); ++it) {
      const auto& quad_point = *it;
      const auto& v = quad_point.position();
      ret += evaluate(v, it.first_index()) * psi(v, is_negative(it)) * quad_point.weight();
    }
    return ret;
  }

  virtual bool is_negative(
      const typename MergedQuadrature<RangeFieldType, dimDomain>::MergedQuadratureIterator& it) const override final
  {
    return it.first_index() < num_intervals / 2;
  }

  // returns matrix with entries <h_i h_j>
  virtual MatrixType mass_matrix(const bool /*use_fine_quadratures*/ = false) const override final
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

  virtual MatrixType mass_matrix_inverse(const bool /*use_fine_quadratures*/ = false) const override final
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

  // returns V M^-1 where the matrix V has entries <v h_i h_j>_- and <v h_i h_j>_+
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

  RangeFieldType min_density(const XT::Common::BlockedFieldVector<RangeFieldType, num_intervals, 2>& u) const
  {
    RangeFieldType ret(u.block(0)[0]);
    for (size_t jj = 1; jj < num_intervals; ++jj)
      ret = std::min(ret, u.block(jj)[0]);
    return ret;
  }

  using BaseType::u_iso;

  // For the partial moments, we might not be able to solve the optimization problem for some moments where the density
  // on one interval/spherical triangle is very low. The overall density might be much higher than the density on that
  // triangle, so we specialize this function.
  virtual void ensure_min_density(RangeType& u, const RangeFieldType min_density) const override final
  {
    const auto u_iso_min = u_iso() * min_density;
    for (size_t jj = 0; jj < num_intervals; ++jj) {
      if (u[2 * jj] < u_iso_min[2 * jj]) {
        u[2 * jj] = u_iso_min[2 * jj];
        u[2 * jj + 1] = u_iso_min[2 * jj + 1];
      }
    }
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
  using BaseType::quadratures_;
}; // class PartialMomentBasis<DomainFieldType, 1, ...>

template <class DomainFieldType, class RangeFieldType, size_t refinements, size_t dimFlux>
class PartialMomentBasis<DomainFieldType, 3, RangeFieldType, refinements, 1, dimFlux, 1>
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
  static constexpr size_t block_size = 4;
  static constexpr size_t num_blocks = dimRange / block_size;

private:
  using BaseType = BasisfunctionsInterface<DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols, dimFlux>;
  using ThisType = PartialMomentBasis;

public:
  typedef SphericalTriangulation<DomainFieldType> TriangulationType;
  using typename BaseType::DomainType;
  using typename BaseType::MatrixType;
  using typename BaseType::QuadraturesType;
  using typename BaseType::QuadratureType;
  using typename BaseType::RangeType;
  using typename BaseType::StringifierType;
  template <class DiscreteFunctionType>
  using VisualizerType = typename BaseType::template VisualizerType<DiscreteFunctionType>;
  using BlockRangeType = XT::Common::FieldVector<RangeFieldType, block_size>;
  using BlockPlaneCoefficientsType = typename std::vector<std::pair<BlockRangeType, RangeFieldType>>;
  using PlaneCoefficientsType = XT::Common::FieldVector<BlockPlaneCoefficientsType, num_blocks>;
  using BlockMatrixType = XT::Common::BlockedFieldMatrix<RangeFieldType, num_blocks, block_size, block_size>;
  using LocalMatrixType = typename BlockMatrixType::BlockType;
  using LocalVectorType = XT::Common::FieldVector<RangeFieldType, block_size>;

  using BaseType::barycentre_rule;

  PartialMomentBasis(const QuadraturesType& quadratures, const QuadraturesType& fine_quadratures)
    : BaseType(refinements, quadratures, fine_quadratures)
  {
    assert(4 * triangulation_.faces().size() == dimRange);
    BaseType::initialize_base_values();
  }

  PartialMomentBasis(const size_t quad_refinements, const QuadratureRule<RangeFieldType, 2>& reference_quadrature_rule)
    : BaseType(refinements)
  {
    quadratures_ = triangulation_.quadrature_rules(quad_refinements, reference_quadrature_rule);
    fine_quadratures_ = triangulation_.quadrature_rules(quad_refinements + 3, reference_quadrature_rule);
    assert(4 * triangulation_.faces().size() == dimRange);
    BaseType::initialize_base_values();
  }

  // This constructor is here for compatibility with the one-dimensional basis to simplify testing
  PartialMomentBasis(const size_t fekete_rule_num = 3,
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
    fine_quadratures_ = triangulation_.quadrature_rules(quad_refinements +
#if HAVE_FEKETE
                                                            3
#else
                                                            0
#endif
                                                        ,
                                                        reference_quadrature_rule);
    assert(4 * triangulation_.faces().size() == dimRange);
    BaseType::initialize_base_values();
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

  virtual RangeFieldType unit_ball_volume(const bool use_fine_quadratures = false) const override final
  {
    return BaseType::unit_ball_volume_quad(use_fine_quadratures);
  }

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

  RangeFieldType min_density(const XT::Common::BlockedFieldVector<RangeFieldType, dimRange / 4, 4>& u) const
  {
    RangeFieldType ret(u.block(0)[0]);
    for (size_t jj = 1; jj < dimRange / 4; ++jj)
      ret = std::min(ret, u.block(jj)[0]);
    return ret;
  }

  using BaseType::u_iso;

  // For the partial moments, we might not be able to solve the optimization problem for some moments where the density
  // on one interval/spherical triangle is very low. The overall density might be much higher than the density on that
  // triangle, so we specialize this function.
  virtual void ensure_min_density(RangeType& u, const RangeFieldType min_density) const override final
  {
    const auto u_iso_min = u_iso() * min_density;
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      const auto block_density = u[4 * jj];
      if (block_density < u_iso_min[4 * jj]) {
        for (size_t ii = 0; ii < block_size; ++ii)
          u[4 * jj + ii] = u_iso_min[4 * jj + ii];
      }
    }
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

  BlockMatrixType block_mass_matrix(const bool use_fine_quadratures = false) const
  {
    BlockMatrixType block_matrix;
    parallel_quadrature_blocked(use_fine_quadratures ? fine_quadratures_ : quadratures_, block_matrix, size_t(-1));
    return block_matrix;
  } // ... mass_matrix()

  virtual MatrixType mass_matrix(const bool use_fine_quadratures = false) const override final
  {
    return block_mass_matrix(use_fine_quadratures).convert_to_dynamic_matrix();
  } // ... mass_matrix()

  virtual FieldVector<MatrixType, dimFlux> mass_matrix_with_v() const override final
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
  virtual FieldVector<FieldVector<MatrixType, 2>, dimFlux> kinetic_flux_matrices() const override final
  {
    FieldVector<FieldVector<MatrixType, 2>, dimFlux> B_kinetic(
        FieldVector<MatrixType, 2>(MatrixType(dimRange, dimRange, 0.)));
    BlockMatrixType block_matrix;
    const BlockMatrixType mass_matrix = block_mass_matrix();
    for (size_t dd = 0; dd < dimFlux; ++dd) {
      QuadraturesType neg_quadratures(quadratures_.size());
      QuadraturesType pos_quadratures(quadratures_.size());
      BaseType::get_pos_and_neg_quadratures(neg_quadratures, pos_quadratures, dd);
      parallel_quadrature_blocked(neg_quadratures, block_matrix, dd);
      apply_invM_from_right(mass_matrix, block_matrix, B_kinetic[dd][0]);
      parallel_quadrature_blocked(pos_quadratures, block_matrix, dd);
      apply_invM_from_right(mass_matrix, block_matrix, B_kinetic[dd][1]);
    } // dd
    return B_kinetic;
  } // ... kinetic_flux_matrices()


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
    DUNE_THROW(Dune::NotImplemented, "You are missing Qhull!");
#endif // HAVE_QHULL
  }

  void
  parallel_quadrature_blocked(const QuadraturesType& quadratures, BlockMatrixType& matrix, const size_t v_index) const
  {
    size_t num_threads = num_blocks;
    std::vector<std::thread> threads(num_threads);
    // Launch a group of threads
    for (size_t jj = 0; jj < num_threads; ++jj)
      threads[jj] = std::thread(&ThisType::calculate_block_in_thread,
                                this,
                                std::cref(quadratures[jj]),
                                std::ref(matrix.block(jj)),
                                v_index,
                                jj);
    // Join the threads with the main thread
    for (size_t jj = 0; jj < num_threads; ++jj)
      threads[jj].join();
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

  using BaseType::fine_quadratures_;
  using BaseType::quadratures_;
  using BaseType::triangulation_;
  mutable PlaneCoefficientsType plane_coefficients_;
}; // class PartialMomentBasis<DomainFieldType, 3, ...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_PARTIALMOMENTS_HH
