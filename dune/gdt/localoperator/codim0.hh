// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LOCALOPERATOR_INTEGRAL_HH
#define DUNE_GDT_LOCALOPERATOR_INTEGRAL_HH

#include <vector>
#include <utility>
#include <type_traits>
#include <limits>

#include <dune/stuff/common/disable_warnings.hh>
#include <dune/common/densematrix.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/common/disable_warnings.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/common/matrix.hh>
#include <dune/stuff/functions/interfaces.hh>

#include "../localevaluation/interface.hh"
#include "interface.hh"

namespace Dune {
namespace GDT {
namespace LocalOperator {


// forward, to be used in the traits
template <class BinaryEvaluationImp>
class Codim0Integral;


template <class BinaryEvaluationImp>
class Codim0IntegralTraits
{
  static_assert(std::is_base_of<LocalEvaluation::Codim0Interface<typename BinaryEvaluationImp::Traits, 2>,
                                BinaryEvaluationImp>::value,
                "BinaryEvaluationImp has to be derived from LocalEvaluation::Codim0Interface< ..., 2 >!");

public:
  typedef Codim0Integral<BinaryEvaluationImp> derived_type;
  typedef LocalEvaluation::Codim0Interface<typename BinaryEvaluationImp::Traits, 2> BinaryEvaluationType;
};


template <class BinaryEvaluationImp>
class Codim0Integral : public LocalOperator::Codim0Interface<Codim0IntegralTraits<BinaryEvaluationImp>>
{
public:
  typedef Codim0IntegralTraits<BinaryEvaluationImp> Traits;
  typedef typename Traits::BinaryEvaluationType BinaryEvaluationType;

private:
  static const size_t numTmpObjectsRequired_ = 1;

public:
  Codim0Integral(const BinaryEvaluationImp eval, const size_t over_integrate = 0)
    : evaluation_(eval)
    , over_integrate_(over_integrate)
  {
  }

  Codim0Integral(const size_t over_integrate, const BinaryEvaluationImp eval)
    : evaluation_(eval)
    , over_integrate_(over_integrate)
  {
  }

  template <class... Args>
  Codim0Integral(Args&&... args)
    : evaluation_(std::forward<Args>(args)...)
    , over_integrate_(0)
  {
  }

  template <class... Args>
  Codim0Integral(const int over_integrate, Args&&... args)
    : evaluation_(std::forward<Args>(args)...)
    , over_integrate_(over_integrate)
  {
  }

  template <class... Args>
  Codim0Integral(const size_t over_integrate, Args&&... args)
    : evaluation_(std::forward<Args>(args)...)
    , over_integrate_(over_integrate)
  {
  }

  const BinaryEvaluationType& inducingEvaluation() const
  {
    return evaluation_;
  }

  size_t numTmpObjectsRequired() const
  {
    return numTmpObjectsRequired_;
  }

  template <class E, class D, int d, class R, int rT, int rCT, int rA, int rCA>
  void apply(const Stuff::LocalfunctionSetInterface<E, D, d, R, rT, rCT>& testBase,
             const Stuff::LocalfunctionSetInterface<E, D, d, R, rA, rCA>& ansatzBase, Dune::DynamicMatrix<R>& ret,
             std::vector<Dune::DynamicMatrix<R>>& tmpLocalMatrices) const
  {
    const auto& entity        = ansatzBase.entity();
    const auto localFunctions = evaluation_.localFunctions(entity);
    // quadrature
    typedef Dune::QuadratureRules<D, d> VolumeQuadratureRules;
    typedef Dune::QuadratureRule<D, d> VolumeQuadratureType;
    const size_t integrand_order = evaluation().order(localFunctions, ansatzBase, testBase) + over_integrate_;
    assert(integrand_order < std::numeric_limits<int>::max());
    const VolumeQuadratureType& volumeQuadrature = VolumeQuadratureRules::rule(entity.type(), int(integrand_order));
    // check matrix and tmp storage
    const size_t rows = testBase.size();
    const size_t cols = ansatzBase.size();
    Dune::Stuff::Common::clear(ret);
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    assert(tmpLocalMatrices.size() >= numTmpObjectsRequired_);
    auto& evaluationResult = tmpLocalMatrices[0];
    // loop over all quadrature points
    const auto quadPointEndIt = volumeQuadrature.end();
    for (auto quadPointIt = volumeQuadrature.begin(); quadPointIt != quadPointEndIt; ++quadPointIt) {
      const Dune::FieldVector<D, d> x = quadPointIt->position();
      // integration factors
      const double integrationFactor = entity.geometry().integrationElement(x);
      const double quadratureWeight  = quadPointIt->weight();
      // evaluate the local operation
      evaluation().evaluate(localFunctions, ansatzBase, testBase, x, evaluationResult);
      // compute integral
      for (size_t ii = 0; ii < rows; ++ii) {
        auto& retRow                    = ret[ii];
        const auto& evaluationResultRow = evaluationResult[ii];
        for (size_t jj = 0; jj < cols; ++jj)
          retRow[jj] += evaluationResultRow[jj] * integrationFactor * quadratureWeight;
      } // compute integral
    } // loop over all quadrature points
  } // ... apply(...)

private:
  const BinaryEvaluationType& evaluation() const
  {
    return static_cast<const BinaryEvaluationType&>(evaluation_);
  }

  const BinaryEvaluationImp evaluation_;
  const size_t over_integrate_;
}; // class Codim0Integral


} // namespace LocalOperator
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALOPERATOR_INTEGRAL_HH
