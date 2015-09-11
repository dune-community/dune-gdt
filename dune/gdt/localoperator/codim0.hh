// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#warning This header is deprecated, include <dune/gdt/localoperator/integrals.hh> instead (08.09.2015)!

#ifndef DUNE_GDT_LOCALOPERATOR_INTEGRAL_HH
#define DUNE_GDT_LOCALOPERATOR_INTEGRAL_HH

#include <vector>
#include <utility>
#include <type_traits>
#include <limits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/deprecated.hh>
#include <dune/common/densematrix.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/stuff/functions/interfaces.hh>

#include "../localevaluation/interface.hh"
#include "interface.hh"

namespace Dune {
namespace GDT {
namespace LocalOperator {


// forward, to be used in the traits
template< class BinaryEvaluationImp >
class Codim0Integral;


namespace internal {


template< class BinaryEvaluationImp >
class Codim0IntegralTraits
{
  static_assert(std::is_base_of< LocalEvaluation::Codim0Interface< typename BinaryEvaluationImp::Traits, 2 >,
                                 BinaryEvaluationImp >::value,
                "BinaryEvaluationImp has to be derived from LocalEvaluation::Codim0Interface< ..., 2 >!");
public:
  typedef Codim0Integral< BinaryEvaluationImp > derived_type;
};


} // namespace internal


template< class BinaryEvaluationType >
class
  DUNE_DEPRECATED_MSG("Use LocalVolumeIntegralOperator instead (08.09.2015)!")
      Codim0Integral
  : public LocalOperator::Codim0Interface< internal::Codim0IntegralTraits< BinaryEvaluationType > >
{
  static const size_t numTmpObjectsRequired_ = 1;
public:
  typedef internal::Codim0IntegralTraits< BinaryEvaluationType > Traits;

  template< class... Args >
  explicit Codim0Integral(Args&& ...args)
    : integrand_(std::forward< Args >(args)...)
    , over_integrate_(0)
  {}

  template< class... Args >
  explicit Codim0Integral(const int over_integrate, Args&& ...args)
    : integrand_(std::forward< Args >(args)...)
    , over_integrate_(boost::numeric_cast< size_t >(over_integrate))
  {}

  template< class... Args >
  explicit Codim0Integral(const size_t over_integrate, Args&& ...args)
    : integrand_(std::forward< Args >(args)...)
    , over_integrate_(over_integrate)
  {}

  size_t numTmpObjectsRequired() const
  {
    return numTmpObjectsRequired_;
  }

  template< class E, class D, size_t d, class R, size_t rT, size_t rCT, size_t rA, size_t rCA >
  void apply(const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
             const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase,
             Dune::DynamicMatrix< R >& ret,
             std::vector< Dune::DynamicMatrix< R > >& tmpLocalMatrices) const
  {
    const auto& entity = ansatzBase.entity();
    const auto localFunctions = integrand_.localFunctions(entity);
    // quadrature
    const size_t integrand_order = integrand_.order(localFunctions, ansatzBase, testBase) + over_integrate_;
    const auto& volumeQuadrature = QuadratureRules< D, d >::rule(entity.type(),
                                                                 boost::numeric_cast< int >(integrand_order));
    // check matrix and tmp storage
    const size_t rows = testBase.size();
    const size_t cols = ansatzBase.size();
    ret *= 0.0;
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    assert(tmpLocalMatrices.size() >= numTmpObjectsRequired_);
    auto& evaluationResult = tmpLocalMatrices[0];
    // loop over all quadrature points
    for (const auto& quadPoint : volumeQuadrature) {
      const auto x = quadPoint.position();
      // integration factors
      const auto integrationFactor = entity.geometry().integrationElement(x);
      const auto quadratureWeight = quadPoint.weight();
      // evaluate the integrand
      integrand_.evaluate(localFunctions, testBase, ansatzBase, x, evaluationResult);
      // compute integral
      for (size_t ii = 0; ii < rows; ++ii) {
        auto& retRow = ret[ii];
        const auto& evaluationResultRow = evaluationResult[ii];
        for (size_t jj = 0; jj < cols; ++jj)
          retRow[jj] += evaluationResultRow[jj] * integrationFactor * quadratureWeight;
      } // compute integral
    } // loop over all quadrature points
  } // ... apply(...)

private:
  const BinaryEvaluationType integrand_;
  const size_t over_integrate_;
}; // class Codim0Integral


} // namespace LocalOperator
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALOPERATOR_INTEGRAL_HH
