// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LOCALFUNCTIONAL_CODIM1_HH
#define DUNE_GDT_LOCALFUNCTIONAL_CODIM1_HH

#include <utility>
#include <vector>
#include <type_traits>
#include <limits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/densevector.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/stuff/functions/interfaces.hh>

#include "../localevaluation/interface.hh"
#include "interface.hh"

namespace Dune {
namespace GDT {
namespace LocalFunctional {


// forward
template< class UnaryEvaluationImp >
class Codim1Integral;


namespace internal {


template< class UnaryEvaluationImp >
class Codim1IntegralTraits
{
  static_assert(std::is_base_of<  LocalEvaluation::Codim1Interface< typename UnaryEvaluationImp::Traits, 1 >,
                                  UnaryEvaluationImp >::value,
                "UnaryEvaluationImp has to be derived from LocalEvaluation::Codim1Interface< ..., 1 >!");
public:
  typedef Codim1Integral< UnaryEvaluationImp > derived_type;
};


} // namespace internal


template< class UnaryEvaluationImp >
class Codim1Integral
    : public LocalFunctional::Codim1Interface< internal::Codim1IntegralTraits< UnaryEvaluationImp > >
{
public:
  typedef internal::Codim1IntegralTraits< UnaryEvaluationImp > Traits;
  typedef UnaryEvaluationImp                                   UnaryEvaluationType;

private:
  static const size_t numTmpObjectsRequired_ = 1;

public:
  template< class... Args >
  explicit Codim1Integral(Args&& ...args)
    : evaluation_(std::forward< Args >(args)...)
    , over_integrate_(0)
  {}

  template< class... Args >
  explicit Codim1Integral(const int over_integrate, Args&& ...args)
    : evaluation_(std::forward< Args >(args)...)
    , over_integrate_(boost::numeric_cast< size_t >(over_integrate))
  {}

  template< class... Args >
  explicit Codim1Integral(const size_t over_integrate, Args&& ...args)
    : evaluation_(std::forward< Args >(args)...)
    , over_integrate_(over_integrate)
  {}

  size_t numTmpObjectsRequired() const
  {
    return numTmpObjectsRequired_;
  }

  template< class E, class IntersectionType, class D, int d, class R, int r, int rC >
  void apply(const Stuff::LocalfunctionSetInterface< E, D, d, R, r, rC >& testBase,
             const IntersectionType& intersection,
             Dune::DynamicVector< R >& ret,
             std::vector< Dune::DynamicVector< R > >& tmpLocalVectors) const
  {
    // local inducing function
    const auto& entity = testBase.entity();
    const auto localFunctions = evaluation_.localFunctions(entity);
    // quadrature
    const auto integrand_order = evaluation_.order(localFunctions, testBase) + over_integrate_;
    const auto& faceQuadrature = QuadratureRules< D, d - 1 >::rule(intersection.type(),
                                                                   boost::numeric_cast< int >(integrand_order));
    // check vector and tmp storage
    ret *= 0.0;
    const size_t size = testBase.size();
    assert(ret.size() >= size);
    assert(tmpLocalVectors.size() >= numTmpObjectsRequired_);
    auto& localVector = tmpLocalVectors[0];
    // loop over all quadrature points
    for (auto quadPoint = faceQuadrature.begin(); quadPoint != faceQuadrature.end(); ++quadPoint) {
      const Dune::FieldVector< D, d - 1 > localPoint = quadPoint->position();
      const auto integrationFactor = intersection.geometry().integrationElement(localPoint);
      const auto quadratureWeight = quadPoint->weight();
      // evaluate local
      evaluation_.evaluate(localFunctions, testBase, intersection, localPoint, localVector);
      // compute integral
      assert(localVector.size() >= size);
      // loop over all test basis functions
      for (size_t ii = 0; ii < size; ++ii)
        ret[ii] += localVector[ii] * integrationFactor * quadratureWeight;
    } // loop over all quadrature points
  } // void apply(...) const

private:
  const UnaryEvaluationType evaluation_;
  const size_t over_integrate_;
}; // class Codim1Integral


} // namespace LocalFunctional
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALFUNCTIONAL_CODIM1_HH
