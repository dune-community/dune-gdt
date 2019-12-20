// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#ifndef DUNE_GDT_MOMENTMODELS_BASISFUNCTIONS_INTERFACE_HH
#define DUNE_GDT_MOMENTMODELS_BASISFUNCTIONS_INTERFACE_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/xt/common/math.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/tuple.hh>

#include <dune/xt/data/quadratures.hh>
#include <dune/xt/data/spherical_quadratures.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/test/momentmodels/triangulation.hh>

namespace Dune {
namespace GDT {


// see https://en.wikipedia.org/wiki/Tridiagonal_matrix#Inversion
template <class FieldType, int rows>
Dune::DynamicMatrix<FieldType> tridiagonal_matrix_inverse(const DynamicMatrix<FieldType>& matrix)
{
  typedef Dune::DynamicMatrix<FieldType> MatrixType;
  size_t cols = rows;
#ifndef NDEBUG
  for (size_t rr = 0; rr < rows; ++rr)
    for (size_t cc = 0; cc < cols; ++cc)
      if ((cc > rr + 1 || cc + 1 < rr) && XT::Common::FloatCmp::ne(matrix[rr][cc], 0.))
        DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong, "Matrix has to be tridiagonal!");
#endif // NDEBUG
  MatrixType ret(rows, rows, 0);
  Dune::FieldVector<FieldType, rows + 1> a(0), b(0), c(0), theta(0);
  Dune::FieldVector<FieldType, rows + 2> phi(0);
  for (size_t ii = 1; ii < rows + 1; ++ii) {
    a[ii] = matrix[ii - 1][ii - 1];
    if (ii < rows) {
      b[ii] = matrix[ii - 1][ii];
      c[ii] = matrix[ii][ii - 1];
    }
  }
  theta[0] = 1;
  theta[1] = a[1];
  for (size_t ii = 2; ii < rows + 1; ++ii)
    theta[ii] = a[ii] * theta[ii - 1] - b[ii - 1] * c[ii - 1] * theta[ii - 2];
  phi[rows + 1] = 1;
  phi[rows] = a[rows];
  for (size_t ii = rows - 1; ii > 0; --ii)
    phi[ii] = a[ii] * phi[ii + 1] - b[ii] * c[ii] * phi[ii + 2];
  for (size_t ii = 1; ii < rows + 1; ++ii) {
    for (size_t jj = 1; jj < cols + 1; ++jj) {
      if (ii == jj)
        ret[ii - 1][jj - 1] = theta[ii - 1] * phi[jj + 1] / theta[rows];
      else if (ii < jj) {
        ret[ii - 1][jj - 1] = std::pow(-1, ii + jj) * theta[ii - 1] * phi[jj + 1] / theta[rows];
        for (size_t kk = ii; kk < jj; ++kk)
          ret[ii - 1][jj - 1] *= b[kk];
      } else if (ii > jj) {
        ret[ii - 1][jj - 1] = std::pow(-1, ii + jj) * theta[jj - 1] * phi[ii + 1] / theta[rows];
        for (size_t kk = jj; kk < ii; ++kk)
          ret[ii - 1][jj - 1] *= c[kk];
      }
    } // jj
  } // ii
#ifndef NDEBUG
  for (size_t ii = 0; ii < rows; ++ii)
    for (size_t jj = 0; jj < cols; ++jj)
      if (std::isnan(ret[ii][jj]) || std::isinf(ret[ii][jj]))
        DUNE_THROW(Dune::MathError, "Inversion of triangular matrix failed!");
#endif
  return ret;
} // ... tridiagonal_matrix_inverse(...)

// After each refinement step:
// num_vertices_new = num_vertices_old + num_intersections_old
// num_intersections_new = 2*num_intersections_old + 3*num_faces_old
// num_faces_new = 4*num_faces_old
// Initially, there are 6 vertices, 12 intersections and 8 faces.
template <size_t refinements>
struct OctaederStatistics
{
  static constexpr size_t constexpr_pow(size_t base, size_t exponent)
  {
    return (exponent == 0) ? 1 : (base * constexpr_pow(base, exponent - 1));
  }

  static constexpr size_t num_faces()
  {
    return 8 * constexpr_pow(4, refinements);
  }

  static constexpr size_t num_intersections()
  {
    return 2 * OctaederStatistics<refinements - 1>::num_intersections()
           + 3 * OctaederStatistics<refinements - 1>::num_faces();
  }

  static constexpr size_t num_vertices()
  {
    return OctaederStatistics<refinements - 1>::num_vertices()
           + OctaederStatistics<refinements - 1>::num_intersections();
  }
};

template <>
struct OctaederStatistics<0>
{
  static constexpr size_t num_faces()
  {
    return 8;
  }

  static constexpr size_t num_intersections()
  {
    return 12;
  }

  static constexpr size_t num_vertices()
  {
    return 6;
  }
};


enum class EntropyType
{
  MaxwellBoltzmann,
  BoseEinstein
};


template <class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1,
          size_t fluxDim = domainDim,
          EntropyType entrpy = EntropyType::MaxwellBoltzmann>
class MomentBasisInterface
{
public:
  static const size_t dimDomain = domainDim;
  static const size_t dimRange = rangeDim;
  static const size_t dimRangeCols = rangeDimCols;
  static const size_t dimFlux = fluxDim;
  static const size_t d = domainDim;
  static const size_t r = rangeDim;
  static const size_t rC = rangeDimCols;
  static const size_t d_flux = fluxDim;
  static const EntropyType entropy = entrpy;
  using D = DomainFieldImp;
  using R = RangeFieldImp;

  using DomainFieldType = DomainFieldImp;
  using DomainType = FieldVector<DomainFieldType, dimDomain>;
  using RangeFieldType = RangeFieldImp;
  using MatrixType = DynamicMatrix<RangeFieldType>;
  using RangeType = typename XT::Functions::RangeTypeSelector<RangeFieldType, dimRange, dimRangeCols>::type;
  using DynamicRangeType = DynamicVector<RangeFieldType>;
  using QuadratureType = QuadratureRule<DomainFieldType, dimDomain>;
  using QuadraturesType = std::vector<QuadratureType>;
  using VisualizerType = XT::Functions::VisualizerInterface<dimRange, dimRangeCols, RangeFieldType>;
  using StringifierType = std::function<std::string(const RangeType&)>;
  using Partitioning1dType = std::vector<RangeFieldType>;
  using SphericalTriangulationType = SphericalTriangulation<RangeFieldType>;
  using MergedQuadratureIterator =
      typename XT::Data::MergedQuadrature<RangeFieldType, dimDomain>::MergedQuadratureIterator;

  static size_t default_quad_order()
  {
    DUNE_THROW(Dune::NotImplemented, "Please overwrite this function in derived classes!");
    return 0;
  }

  static size_t default_quad_refinements()
  {
    return 0;
  }

  MomentBasisInterface(const QuadraturesType& quadratures = QuadraturesType())
    : triangulation_(stored_triangulation_)
    , quadratures_(quadratures)
  {}

  MomentBasisInterface(const size_t refinements, const QuadraturesType& quadratures = QuadraturesType())
    : stored_triangulation_(refinements)
    , triangulation_(stored_triangulation_)
    , quadratures_(quadratures)
  {}

  MomentBasisInterface(const SphericalTriangulationType& triangulation, const QuadraturesType& quadratures)
    : triangulation_(triangulation)
    , quadratures_(quadratures)
  {}

  virtual ~MomentBasisInterface() {}

  virtual DynamicRangeType evaluate(const DomainType& v) const = 0;

  virtual DynamicRangeType evaluate(const DomainType& v, const size_t /*face_index*/) const
  {
    return evaluate(v);
  }

  // returns <b>, where b is the basis functions vector
  virtual const DynamicRangeType& integrated() const
  {
    return integrated_;
  }

  virtual MatrixType mass_matrix() const
  {
    MatrixType M(dimRange, dimRange, 0.);
    parallel_quadrature(quadratures_, M, size_t(-1));
    return M;
  } // ... mass_matrix()

  virtual MatrixType mass_matrix_inverse() const
  {
    auto ret = mass_matrix();
    ret.invert();
    return ret;
  }

  // Get moment vector from distribution function psi. psi also takes an array of bools as argument to decide which
  // value to use at the discontinuities. Is currently ignored in all testcases except the Heaviside test in one
  // dimension, so we only have to test if we are on the positive or negative interval touching 0. For more general
  // testcases, this would have to be adapted.
  virtual DynamicRangeType get_u(const std::function<RangeFieldType(DomainType, std::array<int, dimDomain>)>& psi) const
  {
    DynamicRangeType ret(dimRange, 0.);
    const auto merged_quads = XT::Data::merged_quadrature(quadratures());
    for (auto it = merged_quads.begin(); it != merged_quads.end(); ++it) {
      const auto& quad_point = *it;
      const auto& v = quad_point.position();
      ret += evaluate(v, it.first_index()) * psi(v, has_fixed_sign(it.first_index())) * quad_point.weight();
    }
    return ret;
  }

  // Tests (per coordinate direction) whether the positions in the quadrature with index index are all negative or
  // positive. Returns -1 if the entries are all non-positive, 1 if the entries are all non-negative, and 0 else.
  virtual std::array<int, dimDomain> has_fixed_sign(const size_t /*index*/) const
  {
    std::array<int, dimDomain> ret;
    for (size_t dd = 0; dd < dimDomain; ++dd)
      ret[dd] = 0;
    return ret;
  }

  virtual FieldVector<MatrixType, dimFlux> flux_matrix() const
  {
    FieldVector<MatrixType, dimFlux> B(MatrixType(dimRange, dimRange, 0));
    for (size_t dd = 0; dd < dimFlux; ++dd)
      parallel_quadrature(quadratures_, B[dd], dd);
    return B;
  }

  virtual std::unique_ptr<VisualizerType> visualizer() const
  {
    return std::make_unique<XT::Functions::GenericVisualizer<dimRange, dimRangeCols, RangeFieldType>>(
        1, [this](const int /*comp*/, const RangeType& val) { return this->density(val); });
  }

  // returns V M^-1 where the matrix V has entries <v h_i h_j>_- and <v h_i h_j>_+
  virtual FieldVector<FieldVector<MatrixType, 2>, dimFlux> kinetic_flux_matrices() const
  {
    const auto M = std::make_unique<XT::Common::FieldMatrix<RangeFieldType, dimRange, dimRange>>(mass_matrix());
    FieldVector<FieldVector<MatrixType, 2>, dimFlux> B_kinetic(
        FieldVector<MatrixType, 2>(MatrixType(dimRange, dimRange, 0.)));
    MatrixType tmp_mat(dimRange, dimRange, 0.);
    for (size_t dd = 0; dd < dimFlux; ++dd) {
      QuadraturesType neg_quadratures(quadratures_.size());
      QuadraturesType pos_quadratures(quadratures_.size());
      get_pos_and_neg_quadratures(neg_quadratures, pos_quadratures, dd);
      parallel_quadrature(neg_quadratures, tmp_mat, dd);
      for (size_t rr = 0; rr < dimRange; ++rr)
        M->solve(B_kinetic[dd][0][rr], tmp_mat[rr]);
      parallel_quadrature(pos_quadratures, tmp_mat, dd);
      for (size_t rr = 0; rr < dimRange; ++rr)
        M->solve(B_kinetic[dd][1][rr], tmp_mat[rr]);
    } // dd
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
    parallel_quadrature(quadratures_, ret, direction, true);
    ret.rightmultiply(mass_matrix_inverse());
    // We need the exact reflection matrix to guarantee Q-realizability, the matrix should only contain 0, +-1, so just
    // ensure it does
    for (size_t ii = 0; ii < dimRange; ++ii) {
      for (size_t jj = 0; jj < dimRange; ++jj) {
        if (std::abs(ret[ii][jj]) > 0.99 && std::abs(ret[ii][jj]) < 1.01)
          ret[ii][jj] = ret[ii][jj] / std::abs(ret[ii][jj]);
        else if (std::abs(ret[ii][jj]) < 0.01)
          ret[ii][jj] = 0;
        else
          DUNE_THROW(Dune::MathError, "Invalid reflection matrix!");
      }
    }
    return ret;
  }

  // return alpha s.t. alpha_one * b(v) == 1 for all v
  virtual DynamicRangeType alpha_one() const = 0;

  // returns alpha s.t. the distribution is isotropic and has density rho
  virtual DynamicRangeType alpha_iso(const RangeFieldType rho = 1.) const
  {
    const auto scale_factor = rho / density(integrated());
    if (entropy == EntropyType::MaxwellBoltzmann)
      return alpha_one() * std::log(scale_factor);
    else
      return alpha_one() * std::log(scale_factor / (scale_factor + 1));
  }

  virtual RangeFieldType density(const DynamicRangeType& u) const = 0;

  virtual RangeFieldType density(const RangeType& u) const
  {
    return density(XT::Common::convert_to<DynamicRangeType>(u));
  }

  virtual void ensure_min_density(DynamicRangeType& u, const RangeFieldType min_density) const
  {
    if (density(u) < min_density)
      u = u_iso() * min_density;
  }

  virtual void ensure_min_density(RangeType& u, const RangeFieldType min_density) const
  {
    if (density(u) < min_density) {
      u = u_iso();
      u *= min_density;
    }
  }

  virtual bool adjust_alpha_to_ensure_min_density(RangeType& /*alpha*/, const RangeFieldType /*min_density*/) const
  {
    return false;
  }

  // Volume of integration domain. For the Mn models it is important that u_iso has density 1. If the basis is exactly
  // integrated, we thus use the exact unit ball volume. If the basis is only integrated by quadrature, we have to use
  // <1> as volume to get a density of 1.
  virtual RangeFieldType unit_ball_volume() const
  {
    return unit_ball_volume_exact();
  }

  virtual const DynamicRangeType& u_iso() const
  {
    return u_iso_;
  }

  virtual std::string short_id() const = 0;

  virtual std::string mn_name() const = 0;

  virtual std::string pn_name() const = 0;

  static QuadratureRule<RangeFieldType, 2> barycentre_rule()
  {
    Dune::QuadratureRule<RangeFieldType, 2> ret;
    ret.push_back(Dune::QuadraturePoint<RangeFieldType, 2>({1. / 3., 1. / 3.}, 0.5));
    return ret;
  }

  static Partitioning1dType create_1d_partitioning(const size_t num_intervals)
  {
    Partitioning1dType ret(num_intervals + 1);
    for (size_t ii = 0; ii <= num_intervals; ++ii)
      ret[ii] = -1. + (2. * ii) / num_intervals;
    return ret;
  }

  const QuadraturesType& quadratures() const
  {
    return quadratures_;
  }

  // A Gauss-Lobatto quadrature on each interval
  template <size_t dD = dimDomain>
  static std::enable_if_t<dD == 1, QuadraturesType> gauss_lobatto_quadratures(const size_t num_intervals,
                                                                              const size_t quad_order,
                                                                              const size_t additional_refinements = 0)
  {
    QuadraturesType ret(num_intervals);
    const auto quads_per_interval = std::pow(2, additional_refinements);
    const auto quadrature_boundaries = create_1d_partitioning(num_intervals * quads_per_interval);
    // quadrature on reference interval [0, 1]
    const auto reference_quadrature = XT::Data::GaussLobattoQuadrature<DomainFieldType>::get(quad_order);
    // map to quadrature on interval [a, b] by
    // x_i -> (1-x_i) a + x_i b
    // w_i -> w_i * (b-a)
    for (size_t ii = 0; ii < num_intervals; ++ii) {
      for (size_t jj = 0; jj < quads_per_interval; ++jj) {
        for (const auto& quad_point : reference_quadrature) {
          const auto& x = quad_point.position()[0];
          const auto& a = quadrature_boundaries[ii * quads_per_interval + jj];
          const auto& b = quadrature_boundaries[ii * quads_per_interval + jj + 1];
          const auto pos = (1 - x) * a + x * b;
          const auto weight = quad_point.weight() * (b - a);
          ret[ii].emplace_back(pos, weight);
        } // quad_points
      } // jj
    } // quad_cells
    return ret;
  }

  template <size_t dD = dimDomain>
  static std::enable_if_t<dD == 3, QuadraturesType> lebedev_quadrature(const size_t quad_order)
  {
    return QuadraturesType(1, XT::Data::LebedevQuadrature<DomainFieldType, true>::get(quad_order));
  }

  static RangeFieldType unit_ball_volume_exact()
  {
    if (dimDomain == 1)
      return 2;
    else if (dimDomain == 2)
      return 2 * M_PI;
    else if (dimDomain == 3)
      return 4 * M_PI;
    else {
      DUNE_THROW(NotImplemented, "");
      return 0;
    }
  }

  RangeFieldType unit_ball_volume_quad() const
  {
    RangeFieldType ret(0.);
    for (const auto& quad_point : XT::Data::merged_quadrature(quadratures_))
      ret += quad_point.weight();
    return ret;
  }

protected:
  void initialize_base_values()
  {
    integrated_ = integrated_initializer(quadratures_);
    u_iso_ = integrated() / density(integrated());
  }

  std::array<int, dimDomain> interval_has_fixed_sign(const size_t index, const size_t num_intervals) const
  {
    std::array<int, dimDomain> ret;
    const bool odd = num_intervals % 2;
    if (index < num_intervals / 2)
      ret[0] = -1;
    else if (odd && index == num_intervals / 2)
      ret[0] = 0;
    else
      ret[0] = 1;
    return ret;
  }

  // assumes that each spherical triangle is contained in one octand of the sphere
  std::array<int, dimDomain> triangle_has_fixed_sign(const size_t index) const
  {
    std::array<int, dimDomain> ret;
    const auto center = triangulation_.faces()[index].center();
    for (size_t dd = 0; dd < dimDomain; ++dd)
      ret[dd] = center[dd] < 0. ? -1 : 1;
    return ret;
  }

  void
  get_pos_and_neg_quadratures(QuadraturesType& neg_quadratures, QuadraturesType& pos_quadratures, const size_t dd) const
  {
    for (size_t ii = 0; ii < quadratures_.size(); ++ii) {
      for (const auto& quad_point : quadratures_[ii]) {
        const auto& v = quad_point.position();
        const auto& weight = quad_point.weight();
        // if v[dd] = 0 the quad_point does not contribute to the integral
        if (v[dd] > 0.)
          pos_quadratures[ii].emplace_back(v, weight);
        else if (v[dd] < 0.)
          neg_quadratures[ii].emplace_back(v, weight);
      } // quad_points
    } // quadratures
  }

  static std::vector<MergedQuadratureIterator> create_decomposition(const QuadraturesType& quadratures,
                                                                    const size_t num_threads)
  {
    const size_t size = XT::Data::merged_quadrature(quadratures).size();
    std::vector<MergedQuadratureIterator> ret(num_threads + 1);
    for (size_t ii = 0; ii < num_threads; ++ii)
      ret[ii] = XT::Data::merged_quadrature(quadratures).iterator(size / num_threads * ii);
    ret[num_threads] = XT::Data::merged_quadrature(quadratures).iterator(size);
    return ret;
  }

  virtual void parallel_quadrature(const QuadraturesType& quadratures,
                                   MatrixType& matrix,
                                   const size_t v_index,
                                   const bool reflecting = false) const
  {
    size_t num_threads =
        std::min(XT::Common::threadManager().max_threads(), XT::Data::merged_quadrature(quadratures).size());
    auto decomposition = create_decomposition(quadratures, num_threads);
    std::vector<std::thread> threads(num_threads);
    // Launch a group of threads
    std::vector<MatrixType> local_matrices(num_threads, MatrixType(matrix.N(), matrix.M(), 0.));
    for (size_t ii = 0; ii < num_threads; ++ii)
      threads[ii] = std::thread(&MomentBasisInterface::calculate_in_thread,
                                this,
                                std::ref(local_matrices[ii]),
                                v_index,
                                std::cref(decomposition),
                                ii,
                                reflecting);
    // Join the threads with the main thread
    for (size_t ii = 0; ii < num_threads; ++ii)
      threads[ii].join();
    // add local matrices
    matrix *= 0.;
    for (size_t ii = 0; ii < num_threads; ++ii)
      matrix += local_matrices[ii];
  } // void parallel_quadrature(...)

  virtual void calculate_in_thread(MatrixType& local_matrix,
                                   const size_t v_index,
                                   const std::vector<MergedQuadratureIterator>& decomposition,
                                   const size_t ii,
                                   const bool reflecting) const
  {
    const auto& reflected_indices = triangulation_.reflected_face_indices();
    for (auto it = decomposition[ii]; it != decomposition[ii + 1]; ++it) {
      const auto& quad_point = *it;
      const auto& v = quad_point.position();
      const auto basis_evaluated = evaluate(v, it.first_index());
      auto basis_reflected = basis_evaluated;
      if (reflecting) {
        auto v_reflected = v;
        v_reflected[v_index] *= -1.;
        // If the basis functions have a triangulation, get index of reflected triangle. Otherwise set to 0, will be
        // ignored.
        const size_t reflected_index = reflected_indices.size() ? reflected_indices[it.first_index()][v_index] : 0;
        basis_reflected = evaluate(v_reflected, reflected_index);
      }
      const auto& weight = quad_point.weight();
      const auto factor = (reflecting || v_index == size_t(-1)) ? 1. : v[v_index];
      for (size_t nn = 0; nn < local_matrix.N(); ++nn)
        for (size_t mm = 0; mm < local_matrix.M(); ++mm)
          local_matrix[nn][mm] +=
              basis_evaluated[nn] * (reflecting ? basis_reflected[mm] : basis_evaluated[mm]) * factor * weight;
    }
  } // void calculate_in_thread(...)

  virtual DynamicRangeType integrated_initializer(const QuadraturesType& quadratures) const
  {
    size_t num_threads =
        std::min(XT::Common::threadManager().max_threads(), XT::Data::merged_quadrature(quadratures).size());
    auto decomposition = create_decomposition(quadratures, num_threads);
    std::vector<std::thread> threads(num_threads);
    std::vector<DynamicRangeType> local_vectors(num_threads, DynamicRangeType(dimRange, 0.));
    for (size_t ii = 0; ii < num_threads; ++ii)
      threads[ii] = std::thread(&MomentBasisInterface::integrated_initializer_thread,
                                this,
                                std::ref(local_vectors[ii]),
                                std::cref(decomposition),
                                ii);
    // Join the threads with the main thread
    for (size_t ii = 0; ii < num_threads; ++ii)
      threads[ii].join();
    // add local matrices
    DynamicRangeType ret(dimRange, 0.);
    for (size_t ii = 0; ii < num_threads; ++ii)
      ret += local_vectors[ii];
    return ret;
  }

  virtual void integrated_initializer_thread(DynamicRangeType& local_range,
                                             const std::vector<MergedQuadratureIterator>& decomposition,
                                             const size_t ii) const
  {
    for (auto it = decomposition[ii]; it != decomposition[ii + 1]; ++it) {
      const auto& quad_point = *it;
      auto basis_evaluated = evaluate(quad_point.position(), it.first_index());
      basis_evaluated *= quad_point.weight();
      local_range += basis_evaluated;
    } // jj
  } // void calculate_in_thread(...)

  SphericalTriangulationType stored_triangulation_;
  const SphericalTriangulationType& triangulation_;
  QuadraturesType quadratures_;
  DynamicRangeType integrated_;
  DynamicRangeType u_iso_;
};

template <class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols,
          size_t fluxDim,
          EntropyType entropy>
const size_t
    MomentBasisInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols, fluxDim, entropy>::dimRange;


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_MOMENTMODELS_BASISFUNCTIONS_INTERFACE_HH
