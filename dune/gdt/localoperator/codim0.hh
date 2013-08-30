// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LOCALOPERATOR_INTEGRAL_HH
#define DUNE_GDT_LOCALOPERATOR_INTEGRAL_HH

#include <vector>
#include <utility>
#include <type_traits>

#include <dune/common/densematrix.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/stuff/common/matrix.hh>

#include "../basefunctionset/interface.hh"
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
  Codim0Integral(const BinaryEvaluationImp evaluation)
    : evaluation_(evaluation)
  {
  }

  template <class... Args>
  Codim0Integral(Args&&... args)
    : evaluation_(std::forward<Args>(args)...)
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

  template <class T, class A, class D, int d, class R, int rT, int rCT, int rA, int rCA>
  void apply(const BaseFunctionSetInterface<T, D, d, R, rT, rCT>& testBase,
             const BaseFunctionSetInterface<A, D, d, R, rA, rCA>& ansatzBase, Dune::DynamicMatrix<R>& ret,
             std::vector<Dune::DynamicMatrix<R>>& tmpLocalMatrices) const
  {
    const auto& entity        = ansatzBase.entity();
    const auto localFunctions = evaluation_.localFunctions(entity);
    // quadrature
    typedef Dune::QuadratureRules<D, d> VolumeQuadratureRules;
    typedef Dune::QuadratureRule<D, d> VolumeQuadratureType;
    const int quadratureOrder = evaluation().order(localFunctions, ansatzBase, testBase);
    assert(quadratureOrder >= 0 && "Not implemented for negative integration orders!");
    const VolumeQuadratureType& volumeQuadrature = VolumeQuadratureRules::rule(entity.type(), 2 * quadratureOrder + 1);
    // check matrix and tmp storage
    const size_t rows = testBase.size();
    const size_t cols = ansatzBase.size();
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    assert(tmpLocalMatrices.size() >= numTmpObjectsRequired_);
    // loop over all quadrature points
    const auto quadPointEndIt = volumeQuadrature.end();
    for (auto quadPointIt = volumeQuadrature.begin(); quadPointIt != quadPointEndIt; ++quadPointIt) {
      const Dune::FieldVector<D, d> x = quadPointIt->position();
      // integration factors
      const double integrationFactor = entity.geometry().integrationElement(x);
      const double quadratureWeight  = quadPointIt->weight();
      // clear tmp matrix
      Dune::Stuff::Common::clear(tmpLocalMatrices[0]);
      // evaluate the local operation
      evaluation().evaluate(localFunctions, ansatzBase, testBase, x, tmpLocalMatrices[0]);
      // compute integral
      for (size_t ii = 0; ii < rows; ++ii) {
        auto& retRow       = ret[ii];
        const auto& tmpRow = tmpLocalMatrices[0][ii];
        for (size_t jj = 0; jj < cols; ++jj)
          retRow[jj] += tmpRow[jj] * integrationFactor * quadratureWeight;
      } // compute integral
    } // loop over all quadrature points
  } // ... apply(...)

private:
  const BinaryEvaluationType& evaluation() const
  {
    return static_cast<const BinaryEvaluationType&>(evaluation_);
  }

  const BinaryEvaluationImp evaluation_;
}; // class Codim0Integral


} // namespace LocalOperator
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALOPERATOR_INTEGRAL_HH
