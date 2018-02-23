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

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


// take a DiscreteFunction and return a DiscreteFunction corresponding to component ii
template <size_t ii, class DiscreteFunctionType>
auto get_factor_discrete_function(const DiscreteFunctionType& discrete_function) ->
    typename Dune::GDT::DiscreteFunction<
        typename XT::Common::tuple_element<ii, typename DiscreteFunctionType::SpaceType::SpaceTupleType>::type,
        typename DiscreteFunctionType::VectorType>
{
  typedef typename Dune::GDT::DiscreteFunction<
      typename XT::Common::tuple_element<ii, typename DiscreteFunctionType::SpaceType::SpaceTupleType>::type,
      typename DiscreteFunctionType::VectorType>
      FactorDiscreteFunctionType;
  static_assert(ii < DiscreteFunctionType::SpaceType::num_factors, "This factor does not exist.");
  const auto& space = discrete_function.space();
  const auto& factor_space = space.template factor<ii>();
  typename DiscreteFunctionType::VectorType factor_vector(factor_space.mapper().size());
  const auto it_end = space.grid_layer().template end<0>();
  for (auto it = space.grid_layer().template begin<0>(); it != it_end; ++it) {
    const auto& entity = *it;
    for (size_t jj = 0; jj < factor_space.mapper().numDofs(entity); ++jj)
      factor_vector.set_entry(factor_space.mapper().mapToGlobal(entity, jj),
                              discrete_function.vector().get_entry(space.mapper().mapToGlobal(ii, entity, jj)));
  }
  FactorDiscreteFunctionType factor_discrete_function(factor_space);
  factor_discrete_function.vector() = factor_vector;
  //  typedef Dune::GDT::DiscreteFunctionDataHandle<FactorDiscreteFunctionType> DataHandleType;
  //  DataHandleType handle(factor_discrete_function);
  //  factor_space.grid_layer().template communicate<DataHandleType>(
  //      handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
  return factor_discrete_function;
}

// static for loop to sum components of a DiscreteFunction
template <size_t index, size_t N>
struct static_discrete_function_loop
{
  template <class DiscreteFunctionType>
  static typename DiscreteFunctionType::VectorType sum_vectors(const DiscreteFunctionType& discrete_function)
  {
    return static_discrete_function_loop<index, N / 2>::sum_vectors(discrete_function)
           + static_discrete_function_loop<index + N / 2, N - N / 2>::sum_vectors(discrete_function);
  }

  template <class DiscreteFunctionType>
  static typename DiscreteFunctionType::VectorType
  sum_vectors_divisible_by(const DiscreteFunctionType& discrete_function, const size_t divisor)
  {
    return static_discrete_function_loop<index, N / 2>::sum_vectors_divisible_by(discrete_function, divisor)
           + static_discrete_function_loop<index + N / 2, N - N / 2>::sum_vectors_divisible_by(discrete_function,
                                                                                               divisor);
  }
};

// specialization to end the loop
template <size_t index>
struct static_discrete_function_loop<index, 1>
{
  template <class DiscreteFunctionType>
  static typename DiscreteFunctionType::VectorType sum_vectors(const DiscreteFunctionType& discrete_function)
  {
    return get_factor_discrete_function<index, DiscreteFunctionType>(discrete_function).vector();
  }

  template <class DiscreteFunctionType>
  static typename DiscreteFunctionType::VectorType
  sum_vectors_divisible_by(const DiscreteFunctionType& discrete_function, const size_t divisor)
  {
    if (!(index % divisor))
      return get_factor_discrete_function<index, DiscreteFunctionType>(discrete_function).vector();
    else
      return typename DiscreteFunctionType::VectorType(
          discrete_function.space().template factor<index>().mapper().size(), 0.);
  }
};

// visualizes sum of components of discrete_function
template <class DiscreteFunctionType, size_t dimRange>
void sum_visualizer(const DiscreteFunctionType& u_n, const std::string& filename_prefix, const size_t ii)
{
  auto sum_function = get_factor_discrete_function<0, DiscreteFunctionType>(u_n);
  sum_function.vector() = static_discrete_function_loop<0, dimRange>::sum_vectors(u_n);
  sum_function.visualize(filename_prefix + "_" + Dune::XT::Common::to_string(ii));
}

// visualizes sum of components with index divisible by divisor
template <class DiscreteFunctionType, size_t dimRange>
void sum_divisible_by_visualizer(const DiscreteFunctionType& u_n,
                                 const std::string& filename_prefix,
                                 const size_t ii,
                                 const size_t divisor)
{
  auto sum_function = get_factor_discrete_function<0, DiscreteFunctionType>(u_n);
  sum_function.vector() = static_discrete_function_loop<0, dimRange>::sum_vectors_divisible_by(u_n, divisor);
  sum_function.visualize(filename_prefix + "_" + Dune::XT::Common::to_string(ii));
}

// visualizes factor * component of discrete function
template <class DiscreteFunctionType, size_t dimRange, size_t component>
void component_visualizer(const DiscreteFunctionType& u_n,
                          const std::string& filename_prefix,
                          const size_t ii,
                          const double factor = 1.)
{
  auto u_n_comp = get_factor_discrete_function<component, DiscreteFunctionType>(u_n);
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
  static constexpr size_t num_faces()
  {
    return 8 * (1 << 2 * refinements);
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
  typedef DomainFieldImp DomainFieldType;
  typedef XT::Common::FieldVector<DomainFieldType, dimDomain> DomainType;
  typedef RangeFieldImp RangeFieldType;
  typedef DynamicMatrix<RangeFieldType> MatrixType;
  typedef typename XT::Functions::RangeTypeSelector<RangeFieldType, dimRange, dimRangeCols>::type RangeType;
  template <class DiscreteFunctionType>
  using VisualizerType = typename std::function<void(const DiscreteFunctionType&, const std::string&, const size_t)>;
  using StringifierType = std::function<std::string(const RangeType&)>;
  typedef typename Dune::QuadratureRule<RangeFieldType, dimDomain> QuadratureType;

  virtual ~BasisfunctionsInterface(){};

  virtual RangeType evaluate(const DomainType& v) const = 0;

  virtual RangeType integrated() const = 0;

  virtual MatrixType mass_matrix() const = 0;

  virtual MatrixType mass_matrix_inverse() const = 0;

  virtual FieldVector<MatrixType, dimFlux> mass_matrix_with_v() const = 0;

  static QuadratureRule<RangeFieldType, 2> barycentre_rule()
  {
    Dune::QuadratureRule<RangeFieldType, 2> ret;
    ret.push_back(Dune::QuadraturePoint<RangeFieldType, 2>({1. / 3., 1. / 3.}, 0.5));
    return ret;
  }

protected:
  std::vector<std::vector<size_t>> create_decomposition(const size_t num_threads, const size_t size) const
  {
    std::vector<std::vector<size_t>> decomposition(num_threads);
    for (size_t ii = 0; ii < num_threads - 1; ++ii) {
      decomposition[ii].reserve(size / num_threads * (ii + 1) - size / num_threads * ii);
      for (size_t jj = size / num_threads * ii; jj < size / num_threads * (ii + 1); ++jj)
        decomposition[ii].push_back(jj);
    }
    decomposition[num_threads - 1].reserve(size - (size / num_threads) * (num_threads - 1));
    for (size_t jj = size / num_threads * (num_threads - 1); jj < size; ++jj)
      decomposition[num_threads - 1].push_back(jj);
    return decomposition;
  }

  virtual void parallel_quadrature(const QuadratureType& quadrature,
                                   MatrixType& matrix,
                                   const size_t v_index,
                                   const bool reflecting = false) const
  {
    size_t num_threads = std::min(XT::Common::threadManager().max_threads(), quadrature.size());
    auto decomposition = create_decomposition(num_threads, quadrature.size());
    std::vector<std::thread> threads(num_threads);
    // Launch a group of threads
    std::vector<MatrixType> local_matrices(num_threads, MatrixType(matrix.N(), matrix.M(), 0.));
    for (size_t ii = 0; ii < num_threads; ++ii)
      threads[ii] = std::thread(&BasisfunctionsInterface::calculate_in_thread,
                                this,
                                std::cref(quadrature),
                                std::ref(local_matrices[ii]),
                                v_index,
                                std::cref(decomposition[ii]),
                                reflecting);
    // Join the threads with the main thread
    for (size_t ii = 0; ii < num_threads; ++ii)
      threads[ii].join();
    // add local matrices
    matrix *= 0.;
    for (size_t ii = 0; ii < num_threads; ++ii)
      matrix += local_matrices[ii];
  } // void parallel_quadrature(...)

  virtual void calculate_in_thread(const QuadratureType& quadrature,
                                   MatrixType& local_matrix,
                                   const size_t v_index,
                                   const std::vector<size_t>& indices,
                                   const bool reflecting) const
  {
    for (const auto& jj : indices) {
      const auto& quad_point = quadrature[jj];
      const auto& v = quad_point.position();
      auto v_reflected = v;
      if (reflecting)
        v_reflected[v_index] *= -1.;
      const auto basis_evaluated = evaluate(v);
      const auto basis_reflected = evaluate(v_reflected);
      const auto& weight = quad_point.weight();
      const auto factor = (reflecting || v_index == size_t(-1)) ? 1. : v[v_index];
      for (size_t nn = 0; nn < local_matrix.N(); ++nn)
        for (size_t mm = 0; mm < local_matrix.M(); ++mm)
          local_matrix[nn][mm] +=
              basis_evaluated[nn] * (reflecting ? basis_reflected[mm] : basis_evaluated[mm]) * factor * weight;
    } // ii
  } // void calculate_in_thread(...)

  RangeType integrated_initializer(const QuadratureType& quadrature) const
  {
    size_t num_threads = std::min(XT::Common::threadManager().max_threads(), quadrature.size());
    auto decomposition = create_decomposition(num_threads, quadrature.size());
    std::vector<std::thread> threads(num_threads);
    std::vector<RangeType> local_vectors(num_threads, RangeType(0.));
    for (size_t ii = 0; ii < num_threads; ++ii)
      threads[ii] = std::thread(&BasisfunctionsInterface::integrated_initializer_thread,
                                this,
                                std::cref(quadrature),
                                std::ref(local_vectors[ii]),
                                std::cref(decomposition[ii]));
    // Join the threads with the main thread
    for (size_t ii = 0; ii < num_threads; ++ii)
      threads[ii].join();
    // add local matrices
    RangeType ret(0.);
    for (size_t ii = 0; ii < num_threads; ++ii)
      ret += local_vectors[ii];
    return ret;
  }

  void integrated_initializer_thread(const QuadratureType& quadrature,
                                     RangeType& local_range,
                                     const std::vector<size_t>& indices) const
  {
    for (const auto& jj : indices) {
      const auto& quad_point = quadrature[jj];
      auto basis_evaluated = evaluate(quad_point.position());
      basis_evaluated *= quad_point.weight();
      local_range += basis_evaluated;
    } // jj
  } // void calculate_in_thread(...)
};


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_BASISFUNCTIONS_BASE_HH
