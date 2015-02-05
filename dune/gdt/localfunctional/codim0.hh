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


#include <dune/common/dynmatrix.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/stuff/functions/interfaces.hh>

#include "../localevaluation/interface.hh"
#include "interface.hh"

namespace Dune {
namespace GDT {
namespace LocalFunctional {


// forward, to be used in the traits
template< class UnaryEvaluationImp >
class Codim0Integral;


template< class UnaryEvaluationImp >
class Codim0IntegralTraits
{
  static_assert(std::is_base_of<  LocalEvaluation::Codim0Interface< typename UnaryEvaluationImp::Traits, 1 >,
                                  UnaryEvaluationImp >::value,
                "UnaryEvaluationImp has to be derived from LocalEvaluation::Codim0Interface< ..., 1 >!");
public:
  typedef Codim0Integral< UnaryEvaluationImp >                                        derived_type;
  typedef LocalEvaluation::Codim0Interface< typename UnaryEvaluationImp::Traits, 1 >  UnaryEvaluationType;
};


template< class UnaryEvaluationImp >
class Codim0Integral
  : public LocalFunctional::Codim0Interface< Codim0IntegralTraits< UnaryEvaluationImp > >
{
public:
  typedef Codim0IntegralTraits< UnaryEvaluationImp >  Traits;
  typedef typename Traits::UnaryEvaluationType        UnaryEvaluationType;

private:
  static const size_t numTmpObjectsRequired_ = 1;

public:
  template< class... Args >
  explicit Codim0Integral(Args&& ...args)
    : evaluation_(std::forward< Args >(args)...)
  {}

  Codim0Integral(const UnaryEvaluationImp eval)
    : evaluation_(eval)
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
  template< class E, class D, int d, class R, int r, int rC >
  void apply(const Stuff::LocalfunctionSetInterface< E, D, d, R, r, rC >& testBase,
             Dune::DynamicVector< R >& ret,
             std::vector< Dune::DynamicVector< R > >& tmpLocalVectors) const
  {
    // local inducing function
    const auto& entity = testBase.entity();
    const auto localFunctions = evaluation_.localFunctions(entity);
    // quadrature
    typedef Dune::QuadratureRules< D, d > VolumeQuadratureRules;
    typedef Dune::QuadratureRule< D, d > VolumeQuadratureType;
    const size_t integrand_order = evaluation().order(localFunctions, testBase);
    assert(integrand_order < std::numeric_limits< int >::max());
    const VolumeQuadratureType& volumeQuadrature = VolumeQuadratureRules::rule(entity.type(), int(integrand_order));
    // check vector and tmp storage
    ret *= 0.0;
    const size_t size = testBase.size();
    assert(ret.size() >= size);
    assert(tmpLocalVectors.size() >= numTmpObjectsRequired_);
    Dune::DynamicVector< R >& localVector = tmpLocalVectors[0];
    // loop over all quadrature points
    const auto quadPointEndIt = volumeQuadrature.end();
    for (auto quadPointIt = volumeQuadrature.begin(); quadPointIt != quadPointEndIt; ++quadPointIt) {
      const Dune::FieldVector< D, d > x = quadPointIt->position();
      // integration factors
      const double integrationFactor = entity.geometry().integrationElement(x);
      const double quadratureWeight = quadPointIt->weight();
      // evaluate the local operation
      evaluation().evaluate(localFunctions, testBase, x, localVector);
      // compute integral
      for (size_t ii = 0; ii < size; ++ii) {
        ret[ii] += localVector[ii] * integrationFactor * quadratureWeight;
      } // compute integral
    } // loop over all quadrature points
  } // ... apply(...)

private:
  const UnaryEvaluationType& evaluation() const
  {
    return evaluation_;
  }

  const UnaryEvaluationImp evaluation_;
}; // class Codim0Integral

} // namespace LocalFunctional
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALFUNCTIONAL_CODIM0_HH
