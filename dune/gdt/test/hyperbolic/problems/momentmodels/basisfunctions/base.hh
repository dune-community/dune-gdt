// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_BASISFUNCTIONS_BASE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_BASISFUNCTIONS_BASE_HH

#include <memory>
#include <vector>
#include <string>

#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#include <dune/xt/common/math.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/tuple.hh>

#include <dune/xt/functions/affine.hh>

#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/la/container.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/l2.hh>
#include <dune/gdt/spaces/cg.hh>
#include <dune/gdt/test/hyperbolic/problems/momentmodels/triangulation.hh>
#include <dune/gdt/test/hyperbolic/quadratures/gausslobatto.hh>
#include <dune/gdt/test/hyperbolic/spherical_quadratures/fekete.hh>
#include <dune/gdt/test/hyperbolic/spherical_quadratures/lebedev.hh>
#include <dune/gdt/test/hyperbolic/spherical_quadratures/octant.hh>
#include <dune/gdt/test/hyperbolic/spherical_quadratures/wrapper.hh>

namespace Dune {
namespace GDT {


template <class DiscreteFunctionType>
auto get_factor_space(const DiscreteFunctionType& discrete_function)
    -> FvSpace<typename DiscreteFunctionType::SpaceType::GridLayerType,
               typename DiscreteFunctionType::RangeFieldType,
               1,
               1>
{
  using FactorSpaceType = FvSpace<typename DiscreteFunctionType::SpaceType::GridLayerType,
                                  typename DiscreteFunctionType::RangeFieldType,
                                  1,
                                  1>;
  return FactorSpaceType(discrete_function.space().grid_layer());
}

// take a FiniteVolume-DiscreteFunction and return a DiscreteFunction corresponding to component ii
template <class FactorSpaceType, class DiscreteFunctionType>
auto get_factor_discrete_function(const size_t ii,
                                  const FactorSpaceType& factor_space,
                                  const DiscreteFunctionType& discrete_function)
    -> DiscreteFunction<FactorSpaceType, typename DiscreteFunctionType::VectorType>
{
  assert(ii < DiscreteFunctionType::dimRange);
  using VectorType = typename DiscreteFunctionType::VectorType;
  using FactorDiscreteFunctionType = typename Dune::GDT::DiscreteFunction<FactorSpaceType, VectorType>;
  const auto& space = discrete_function.space();
  FactorDiscreteFunctionType factor_discrete_function(factor_space);
  auto& factor_vector = factor_discrete_function.vector();
  const auto it_end = space.grid_layer().template end<0>();
  for (auto it = space.grid_layer().template begin<0>(); it != it_end; ++it) {
    const auto& entity = *it;
    factor_vector.set_entry(factor_space.mapper().mapToGlobal(entity, 0),
                            discrete_function.vector().get_entry(space.mapper().mapToGlobal(entity, ii)));
  }
  return factor_discrete_function;
}

// visualizes sum of components of discrete_function
template <class DiscreteFunctionType>
void sum_visualizer(const DiscreteFunctionType& u_n, const std::string& filename_prefix, const size_t ii)
{
  auto factor_space = get_factor_space(u_n);
  auto sum_function = get_factor_discrete_function(0, factor_space, u_n);
  for (size_t jj = 1; jj < DiscreteFunctionType::dimRange; ++jj)
    sum_function.vector() += get_factor_discrete_function(jj, factor_space, u_n).vector();
  sum_function.visualize(filename_prefix + "_" + Dune::XT::Common::to_string(ii));
}

// visualizes sum of components with index divisible by divisor
template <class DiscreteFunctionType, size_t dimRange>
void sum_divisible_by_visualizer(const DiscreteFunctionType& u_n,
                                 const std::string& filename_prefix,
                                 const size_t ii,
                                 const size_t divisor)
{
  auto factor_space = get_factor_space(u_n);
  auto sum_function = get_factor_discrete_function(0, factor_space, u_n);
  for (size_t jj = divisor; jj < DiscreteFunctionType::dimRange; jj += divisor)
    sum_function.vector() += get_factor_discrete_function(jj, factor_space, u_n).vector();
  sum_function.visualize(filename_prefix + "_" + Dune::XT::Common::to_string(ii));
}

// visualizes factor * component of discrete function
template <class DiscreteFunctionType>
void component_visualizer(const size_t component,
                          const DiscreteFunctionType& u_n,
                          const std::string& filename_prefix,
                          const size_t ii,
                          const double factor = 1.)
{
  auto factor_space = get_factor_space(u_n);
  auto u_n_comp = get_factor_discrete_function(component, factor_space, u_n);
  u_n_comp.vector() *= factor;
  u_n_comp.visualize(filename_prefix + "_" + Dune::XT::Common::to_string(ii));
}

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


template <class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1,
          size_t fluxDim = domainDim>
class BasisfunctionsInterface
{
public:
  static const size_t dimDomain = domainDim;
  static const size_t dimRange = rangeDim;
  static const size_t dimRangeCols = rangeDimCols;
  static const size_t dimFlux = fluxDim;
  using DomainFieldType = DomainFieldImp;
  using DomainType = XT::Common::FieldVector<DomainFieldType, dimDomain>;
  using RangeFieldType = RangeFieldImp;
  using MatrixType = DynamicMatrix<RangeFieldType>;
  using RangeType = typename XT::Functions::RangeTypeSelector<RangeFieldType, dimRange, dimRangeCols>::type;
  using DynamicRangeType = DynamicVector<RangeFieldType>;
  using QuadraturesType = QuadraturesWrapper<DomainFieldType, dimDomain>;
  using QuadratureType = typename QuadraturesType::QuadratureType;
  template <class DiscreteFunctionType>
  using VisualizerType = std::function<void(const DiscreteFunctionType&, const std::string&, const size_t)>;
  using StringifierType = std::function<std::string(const RangeType&)>;
  using Triangulation1dType = std::vector<RangeFieldType>;
  using SphericalTriangulationType = SphericalTriangulation<RangeFieldType>;
  using MergedQuadratureIterator = typename QuadraturesType::MergedQuadratureType::ConstIteratorType;

  BasisfunctionsInterface(const QuadraturesType& quadratures = QuadraturesType(),
                          const QuadraturesType& fine_quadratures = QuadraturesType())
    : quadratures_(quadratures)
    , fine_quadratures_(fine_quadratures)
  {}

  BasisfunctionsInterface(const size_t refinements,
                          const QuadraturesType& quadratures = QuadraturesType(),
                          const QuadraturesType& fine_quadratures = QuadraturesType())
    : quadratures_(quadratures)
    , fine_quadratures_(fine_quadratures)
    , triangulation_(refinements)
  {}

  virtual ~BasisfunctionsInterface() {}

  virtual DynamicRangeType evaluate(const DomainType& v) const = 0;

  virtual DynamicRangeType evaluate(const DomainType& v, const size_t /*face_index*/) const
  {
    return evaluate(v);
  }

  // returns <b>, where b is the basis functions vector
  virtual DynamicRangeType integrated(const bool use_fine_quadratures = false) const
  {
    if (use_fine_quadratures)
      return integrated_initializer(fine_quadratures_);
    return integrated_;
  }

  virtual MatrixType mass_matrix(const bool use_fine_quadratures = false) const
  {
    MatrixType M(dimRange, dimRange, 0.);
    parallel_quadrature(use_fine_quadratures ? fine_quadratures_ : quadratures_, M, size_t(-1));
    return M;
  } // ... mass_matrix()

  virtual MatrixType mass_matrix_inverse(bool use_fine_quadratures = false) const
  {
    auto ret = mass_matrix(use_fine_quadratures);
    ret.invert();
    return ret;
  }

  virtual DynamicRangeType get_moment_vector(const std::function<RangeFieldType(DomainType, bool)>& psi) const
  {
    DynamicRangeType ret(dimRange, 0.);
    const auto merged_quads = quadratures().merged();
    for (auto it = merged_quads.begin(); it != merged_quads.end(); ++it) {
      const auto& quad_point = *it;
      const auto& v = quad_point.position();
      ret += evaluate(v, it.first_index()) * psi(v, is_negative(it)) * quad_point.weight();
    }
    return ret;
  }

  virtual bool
  is_negative(const typename MergedQuadrature<RangeFieldType, dimDomain>::MergedQuadratureIterator& /*it*/) const
  {
    return false;
  }

  virtual FieldVector<MatrixType, dimFlux> mass_matrix_with_v() const
  {
    FieldVector<MatrixType, dimFlux> B(MatrixType(dimRange, dimRange, 0));
    for (size_t dd = 0; dd < dimFlux; ++dd)
      parallel_quadrature(quadratures_, B[dd], dd);
    return B;
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

  virtual DynamicRangeType alpha_iso() const = 0;

  virtual RangeFieldType density(const DynamicRangeType& u) const = 0;

  virtual void ensure_min_density(DynamicRangeType& u, const RangeFieldType min_density) const
  {
    if (density(u) < min_density)
      u = u_iso() * min_density;
  }

  // Volume of integration domain. For the Mn models it is important that u_iso has density 1. If the basis is exactly
  // integrated, we thus use the exact unit ball volume. If the basis is only integrated by quadrature, we have to use
  // <1> as volume to get a density of 1.
  virtual RangeFieldType unit_ball_volume(const bool /*use_fine_quadratures*/ = false) const
  {
    return unit_ball_volume_exact();
  }

  virtual DynamicRangeType u_iso() const
  {
    return u_iso_;
  }

  virtual std::string short_id() const = 0;

  static QuadratureRule<RangeFieldType, 2> barycentre_rule()
  {
    Dune::QuadratureRule<RangeFieldType, 2> ret;
    ret.push_back(Dune::QuadraturePoint<RangeFieldType, 2>({1. / 3., 1. / 3.}, 0.5));
    return ret;
  }

  static Triangulation1dType create_1d_triangulation(const size_t num_intervals)
  {
    Triangulation1dType ret(num_intervals + 1);
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
                                                                              const size_t quad_order)
  {
    QuadraturesType ret(num_intervals);
    auto interval_boundaries = create_1d_triangulation(num_intervals);
    // quadrature on reference interval [0, 1]
    const auto reference_quadrature = GaussLobattoQuadrature<DomainFieldType>::get(quad_order);
    // map to quadrature on interval [a, b] by
    // x_i -> (1-x_i) a + x_i b
    // w_i -> w_i * (b-a)
    for (size_t ii = 0; ii < num_intervals; ++ii) {
      for (const auto& quad_point : reference_quadrature) {
        const auto& x = quad_point.position()[0];
        const auto& a = interval_boundaries[ii];
        const auto& b = interval_boundaries[ii + 1];
        const auto pos = (1 - x) * a + x * b;
        const auto weight = quad_point.weight() * (b - a);
        ret[ii].emplace_back(pos, weight);
      } // quad_points
    } // quad_cells
    return ret;
  }

  template <size_t dD = dimDomain>
  static std::enable_if_t<dD == 3, QuadraturesType> lebedev_quadrature(const size_t quad_order)
  {
    return QuadraturesType(1, LebedevQuadrature<DomainFieldType, true>::get(quad_order));
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

  RangeFieldType unit_ball_volume_quad(const bool use_fine_quadratures = false) const
  {
    RangeFieldType ret(0.);
    for (const auto& quad_point : (use_fine_quadratures ? fine_quadratures_ : quadratures_).merged())
      ret += quad_point.weight();
    return ret;
  }

protected:
  void initialize_base_values()
  {
    integrated_ = integrated_initializer(quadratures_);
    u_iso_ = integrated() / density(integrated());
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
    const size_t size = quadratures.merged().size();
    std::vector<MergedQuadratureIterator> ret(num_threads + 1);
    for (size_t ii = 0; ii < num_threads; ++ii)
      ret[ii] = quadratures.merged().iterator(size / num_threads * ii);
    ret[num_threads] = quadratures.merged().iterator(size);
    return ret;
  }

  virtual void parallel_quadrature(const QuadraturesType& quadratures,
                                   MatrixType& matrix,
                                   const size_t v_index,
                                   const bool reflecting = false) const
  {
    size_t num_threads = std::min(XT::Common::threadManager().max_threads(), quadratures.merged().size());
    auto decomposition = create_decomposition(quadratures, num_threads);
    std::vector<std::thread> threads(num_threads);
    // Launch a group of threads
    std::vector<MatrixType> local_matrices(num_threads, MatrixType(matrix.N(), matrix.M(), 0.));
    for (size_t ii = 0; ii < num_threads; ++ii)
      threads[ii] = std::thread(&BasisfunctionsInterface::calculate_in_thread,
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

  DynamicRangeType integrated_initializer(const QuadraturesType& quadratures) const
  {
    size_t num_threads = std::min(XT::Common::threadManager().max_threads(), quadratures.merged().size());
    auto decomposition = create_decomposition(quadratures, num_threads);
    std::vector<std::thread> threads(num_threads);
    std::vector<DynamicRangeType> local_vectors(num_threads, DynamicRangeType(dimRange, 0.));
    for (size_t ii = 0; ii < num_threads; ++ii)
      threads[ii] = std::thread(&BasisfunctionsInterface::integrated_initializer_thread,
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

  void integrated_initializer_thread(DynamicRangeType& local_range,
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

  QuadraturesType quadratures_;
  QuadraturesType fine_quadratures_;
  SphericalTriangulationType triangulation_;
  DynamicRangeType integrated_;
  DynamicRangeType u_iso_;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_BASISFUNCTIONS_BASE_HH
