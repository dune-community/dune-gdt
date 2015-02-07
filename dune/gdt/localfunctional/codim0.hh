// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LOCALFUNCTIONAL_CODIM0_HH
#define DUNE_GDT_LOCALFUNCTIONAL_CODIM0_HH

#include <vector>
#include <utility>
#include <type_traits>
#include <limits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/dynmatrix.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/stuff/functions/interfaces.hh>

#include "../localevaluation/interface.hh"
#include "interface.hh"

namespace Dune {
namespace GDT {
namespace LocalFunctional {


// forward
template< class UnaryEvaluationImp >
class Codim0Integral;


namespace internal {


template< class UnaryEvaluationImp >
class Codim0IntegralTraits
{
  static_assert(std::is_base_of< LocalEvaluation::Codim0Interface< typename UnaryEvaluationImp::Traits, 1 >,
                                 UnaryEvaluationImp >::value,
                "UnaryEvaluationImp has to be derived from LocalEvaluation::Codim0Interface< ..., 1 >!");
public:
  typedef Codim0Integral< UnaryEvaluationImp > derived_type;
};


} // namespace internal


template< class UnaryEvaluationImp >
class Codim0Integral
  : public LocalFunctional::Codim0Interface< internal::Codim0IntegralTraits< UnaryEvaluationImp > >
{
public:
  typedef internal::Codim0IntegralTraits< UnaryEvaluationImp > Traits;
  typedef UnaryEvaluationImp                                   UnaryEvaluationType;

private:
  static const size_t numTmpObjectsRequired_ = 1;

public:
  template< class... Args >
  explicit Codim0Integral(Args&& ...args)
    : evaluation_(std::forward< Args >(args)...)
    , over_integrate_(0)
  {}

  template< class... Args >
  explicit Codim0Integral(const int over_integrate, Args&& ...args)
    : evaluation_(std::forward< Args >(args)...)
    , over_integrate_(boost::numeric_cast< size_t >(over_integrate))
  {}

  template< class... Args >
  explicit Codim0Integral(const size_t over_integrate, Args&& ...args)
    : evaluation_(std::forward< Args >(args)...)
    , over_integrate_(over_integrate)
  {}

  size_t numTmpObjectsRequired() const
  {
    return numTmpObjectsRequired_;
  }

  /**
   *  \brief      Applies the local functional.
   *  \tparam T   Traits of the BaseFunctionSetInterface implementation, representing the type of the testBase
   *  \attention  ret is assumed to be zero!
   */
  template< class E, class D, size_t d, class R, size_t r, size_t rC >
  void apply(const Stuff::LocalfunctionSetInterface< E, D, d, R, r, rC >& testBase,
             Dune::DynamicVector< R >& ret,
             std::vector< Dune::DynamicVector< R > >& tmpLocalVectors) const
  {
    // local inducing function
    const auto& entity = testBase.entity();
    const auto localFunctions = evaluation_.localFunctions(entity);
    // quadrature
    const auto integrand_order = evaluation_.order(localFunctions, testBase) + over_integrate_;
    const auto& volumeQuadrature = QuadratureRules< D, d >::rule(entity.type(),
                                                                 boost::numeric_cast< int >(integrand_order));
    // check vector and tmp storage
    ret *= 0.0;
    const size_t size = testBase.size();
    assert(ret.size() >= size);
    assert(tmpLocalVectors.size() >= numTmpObjectsRequired_);
    auto& localVector = tmpLocalVectors[0];
    // loop over all quadrature points
    const auto quadPointEndIt = volumeQuadrature.end();
    for (auto quadPointIt = volumeQuadrature.begin(); quadPointIt != quadPointEndIt; ++quadPointIt) {
      const Dune::FieldVector< D, d > x = quadPointIt->position();
      // integration factors
      const auto integrationFactor = entity.geometry().integrationElement(x);
      const auto quadratureWeight = quadPointIt->weight();
      // evaluate the local operation
      evaluation_.evaluate(localFunctions, testBase, x, localVector);
      // compute integral
      for (size_t ii = 0; ii < size; ++ii)
        ret[ii] += localVector[ii] * integrationFactor * quadratureWeight;
    } // loop over all quadrature points
  } // ... apply(...)

private:
  const UnaryEvaluationType evaluation_;
  const size_t over_integrate_;
}; // class Codim0Integral

} // namespace LocalFunctional
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALFUNCTIONAL_CODIM0_HH
