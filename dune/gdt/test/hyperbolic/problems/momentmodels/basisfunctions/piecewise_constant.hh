// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_PIECEWISECONSTANT_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_PIECEWISECONSTANT_HH

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
class PiecewiseConstant
{
  //  static_assert(false, "Not implemented for this combination of dimension and order!");
};

template <class DomainFieldType, class RangeFieldType, size_t rangeDim, size_t rangeDimCols, size_t dimFlux>
class PiecewiseConstant<DomainFieldType, 1, RangeFieldType, rangeDim, rangeDimCols, dimFlux, 1>
    : public BasisfunctionsInterface<DomainFieldType, 1, RangeFieldType, rangeDim, rangeDimCols, dimFlux>
{
public:
  static const size_t dimDomain = 1;
  static const size_t dimRange = rangeDim;
  static const size_t dimRangeCols = rangeDimCols;

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

  PiecewiseConstant(const TriangulationType& triangulation = BaseType::create_1d_triangulation(dimRange),
                    const QuadraturesType& /*quadratures*/ = QuadraturesType())
    : triangulation_(triangulation)
  {
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
    for (size_t ii = 0; ii < dimRange; ++ii) {
      if (XT::Common::FloatCmp::eq(v[0], triangulation_[ii]))
        boundary = true;
      if (XT::Common::FloatCmp::ge(v[0], triangulation_[ii])
          && XT::Common::FloatCmp::le(v[0], triangulation_[ii + 1])) {
        ret[ii] = 1;
      }
    }
    if (XT::Common::FloatCmp::eq(v[0], triangulation_[dimRange]))
      boundary = true;
    if (split_boundary && boundary)
      ret /= 2.;
    return ret;
  } // ... evaluate(...)

  virtual RangeType integrated() const override final
  {
    RangeType ret(0);
    for (size_t ii = 0; ii < dimRange; ++ii)
      ret[ii] = triangulation_[ii + 1] - triangulation_[ii];
    return ret;
  }

  // returns matrix with entries <h_i h_j>
  virtual MatrixType mass_matrix() const override
  {
    MatrixType M(dimRange, dimRange, 0.);
    for (size_t ii = 0; ii < dimRange; ++ii)
      M[ii][ii] = triangulation_[ii + 1] - triangulation_[ii];
    return M;
  }

  virtual MatrixType mass_matrix_inverse() const override
  {
    auto ret = mass_matrix();
    for (size_t ii = 0; ii < dimRange; ++ii)
      ret[ii][ii] = 1. / ret[ii][ii];
    return ret;
  }

  // returns matrix with entries <v h_i h_j>
  virtual FieldVector<MatrixType, dimDomain> mass_matrix_with_v() const override
  {
    MatrixType B(dimRange, dimRange, 0.);
    for (size_t ii = 0; ii < dimRange; ++ii)
      B[ii][ii] = (std::pow(triangulation_[ii + 1], 2) - std::pow(triangulation_[ii], 2)) / 2.;
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
      if (dimRange % 2) {
        // if there is an odd number of intervals, the middle interval is split in 2 parts
        // copy other intervals
        for (size_t nn = 0; nn < dimRange / 2; ++nn)
          for (size_t mm = 0; mm < dimRange / 2; ++mm)
            ret_neg[nn][mm] = mm_with_v[0][nn][mm];
        for (size_t nn = dimRange / 2 + 1; nn < dimRange; ++nn)
          for (size_t mm = dimRange / 2 + 1; mm < dimRange; ++mm)
            ret_pos[nn][mm] = mm_with_v[0][nn][mm];
        ret_neg[dimRange / 2][dimRange / 2] = -std::pow(triangulation_[dimRange / 2], 2) / 2;
        ret_pos[dimRange / 2][dimRange / 2] = std::pow(triangulation_[dimRange / 2], 2) / 2;
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
      for (size_t jj = 0; jj < dimRange; ++jj) {
        ret[ii][jj] = mass_mat[ii][dimRange - 1 - jj];
      }
    }
    ret.rightmultiply(mass_matrix_inverse());
    return ret;
  }

  template <class DiscreteFunctionType>
  VisualizerType<DiscreteFunctionType> visualizer() const
  {
    return [](const DiscreteFunctionType& u_n, const std::string& filename_prefix, const size_t ii) {
      sum_visualizer<DiscreteFunctionType, dimRange>(u_n, filename_prefix, ii);
    };
  }

  static StringifierType stringifier()
  {
    return [](const RangeType& val) {
      RangeFieldType psi(0);
      for (size_t ii = 0; ii < dimRange; ++ii)
        psi += val[ii];
      return XT::Common::to_string(psi, 15);
    };
  } // ... stringifier()

  const TriangulationType& triangulation() const
  {
    return triangulation_;
  }

  RangeType alpha_iso()
  {
    return RangeType(1.);
  }

  RangeFieldType density(const RangeType& u) const
  {
    return std::accumulate(u.begin(), u.end(), RangeFieldType(0.));
  }

  virtual std::string short_id() const override final
  {
    return "pconst";
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
}; // class PiecewiseConstant<DomainFieldType, 1, ...>


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_PIECEWISECONSTANT_HH
