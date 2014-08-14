// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LOCALOPERATOR_CODIM1_HH
#define DUNE_GDT_LOCALOPERATOR_CODIM1_HH

#include <utility>
#include <vector>
#include <type_traits>
#include <limits>

#include <dune/stuff/common/disable_warnings.hh>
# include <dune/common/densematrix.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/common/disable_warnings.hh>
# include <dune/geometry/quadraturerules.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/common/matrix.hh>
#include <dune/stuff/functions/interfaces.hh>

#include "../localevaluation/interface.hh"
#include "interface.hh"

namespace Dune {
namespace GDT {
namespace LocalOperator {


// forward, to be used in the traits
template< class QuaternaryEvaluationImp >
class Codim1CouplingIntegral;


template< class QuaternaryEvaluationImp >
class Codim1CouplingIntegralTraits
{
  static_assert(std::is_base_of<  LocalEvaluation::Codim1Interface< typename QuaternaryEvaluationImp::Traits, 4 >,
                                  QuaternaryEvaluationImp >::value,
                "QuaternaryEvaluationImp has to be derived from LocalEvaluation::Codim1Interface< ..., 4 >!");
public:
  typedef Codim1CouplingIntegral< QuaternaryEvaluationImp > derived_type;
  typedef LocalEvaluation::Codim1Interface< typename QuaternaryEvaluationImp::Traits, 4 > QuaternaryEvaluationType;
};


template< class QuaternaryEvaluationImp >
class Codim1CouplingIntegral
    : public LocalOperator::Codim1CouplingInterface< Codim1CouplingIntegralTraits< QuaternaryEvaluationImp > >
{
public:
  typedef Codim1CouplingIntegralTraits< QuaternaryEvaluationImp > Traits;
  typedef typename Traits::QuaternaryEvaluationType QuaternaryEvaluationType;

private:
  static const size_t numTmpObjectsRequired_ = 4;

public:
  Codim1CouplingIntegral(const QuaternaryEvaluationImp eval)
    : evaluation_(eval)
  {}

  template< class... Args >
  explicit Codim1CouplingIntegral(Args&& ...args)
    : evaluation_(std::forward< Args >(args)...)
  {}

  size_t numTmpObjectsRequired() const
  {
    return numTmpObjectsRequired_;
  }

  template< class E, class N, class IntersectionType, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  void apply(const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& entityTestBase,
             const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& entityAnsatzBase,
             const Stuff::LocalfunctionSetInterface< N, D, d, R, rT, rCT >& neighborTestBase,
             const Stuff::LocalfunctionSetInterface< N, D, d, R, rA, rCA >& neighborAnsatzBase,
             const IntersectionType& intersection,
             Dune::DynamicMatrix< R >& entityEntityRet,
             Dune::DynamicMatrix< R >& neighborNeighborRet,
             Dune::DynamicMatrix< R >& entityNeighborRet,
             Dune::DynamicMatrix< R >& neighborEntityRet,
             std::vector< Dune::DynamicMatrix< R > >& tmpLocalMatrices) const
  {
    // local inducing function
    const auto& entity = entityTestBase.entity();
    const auto localFunctionsEn = evaluation_.localFunctions(entity);
    const auto& neighbor = neighborTestBase.entity();
    const auto localFunctionsNe = evaluation_.localFunctions(neighbor);
    // quadrature
    const size_t integrand_order = evaluation().order(localFunctionsEn, localFunctionsNe,
                                                       entityTestBase, entityAnsatzBase,
                                                       neighborTestBase, neighborAnsatzBase);
    assert(integrand_order < std::numeric_limits< int >::max());
    const auto& faceQuadrature = QuadratureRules< D, d - 1 >::rule(intersection.type(), int(integrand_order));
    // check matrices
    entityEntityRet *= 0.0;
    neighborNeighborRet *= 0.0;
    entityNeighborRet *= 0.0;
    neighborEntityRet *= 0.0;
    const size_t rowsEn = entityTestBase.size();
    const size_t colsEn = entityAnsatzBase.size();
    const size_t rowsNe = neighborTestBase.size();
    const size_t colsNe = neighborAnsatzBase.size();
    assert(entityEntityRet.rows() >= rowsEn);
    assert(entityEntityRet.cols() >= colsEn);
    assert(neighborNeighborRet.rows() >= rowsNe);
    assert(neighborNeighborRet.cols() >= colsNe);
    assert(entityNeighborRet.rows() >= rowsEn);
    assert(entityNeighborRet.cols() >= colsNe);
    assert(neighborEntityRet.rows() >= rowsEn);
    assert(neighborEntityRet.cols() >= colsEn);
    assert(tmpLocalMatrices.size() >= numTmpObjectsRequired_);
    Dune::DynamicMatrix< R >& entityEntityVals = tmpLocalMatrices[0];
    Dune::DynamicMatrix< R >& neighborNeighborVals = tmpLocalMatrices[1];
    Dune::DynamicMatrix< R >& entityNeighborVals = tmpLocalMatrices[2];
    Dune::DynamicMatrix< R >& neighborEntityVals = tmpLocalMatrices[3];
    // loop over all quadrature points
    for (auto quadPoint = faceQuadrature.begin(); quadPoint != faceQuadrature.end(); ++quadPoint) {
      const Dune::FieldVector< D, d - 1 > localPoint = quadPoint->position();
      const R integrationFactor = intersection.geometry().integrationElement(localPoint);
      const R quadratureWeight = quadPoint->weight();
      // evaluate local
      evaluation().evaluate(localFunctionsEn, localFunctionsNe,
                            entityTestBase, entityAnsatzBase,
                            neighborTestBase, neighborAnsatzBase,
                            intersection, localPoint,
                            entityEntityVals,
                            neighborNeighborVals,
                            entityNeighborVals,
                            neighborEntityVals);
      // compute integral
      assert(entityEntityVals.rows() >= rowsEn);
      assert(entityEntityVals.cols() >= colsEn);
      assert(neighborNeighborVals.rows() >= rowsNe);
      assert(neighborNeighborVals.cols() >= colsNe);
      assert(entityNeighborVals.rows() >= rowsEn);
      assert(entityNeighborVals.cols() >= colsNe);
      assert(neighborEntityVals.rows() >= rowsEn);
      assert(neighborEntityVals.cols() >= colsEn);
      // loop over all entity test basis functions
      for (size_t ii = 0; ii < rowsEn; ++ii) {
        auto& entityEntityRetRow = entityEntityRet[ii];
        const auto& entityEntityValsRow = entityEntityVals[ii];
        auto& entityNeighborRetRow = entityNeighborRet[ii];
        const auto& entityNeighborValsRow = entityNeighborVals[ii];
        // loop over all entity ansatz basis functions
        for (size_t jj = 0; jj < colsEn; ++jj) {
          entityEntityRetRow[jj] += entityEntityValsRow[jj] * integrationFactor * quadratureWeight;
        } // loop over all entity ansatz basis functions
        // loop over all neighbor ansatz basis functions
        for (size_t jj = 0; jj < colsNe; ++jj) {
          entityNeighborRetRow[jj] += entityNeighborValsRow[jj] * integrationFactor * quadratureWeight;
        } // loop over all neighbor ansatz basis functions
      } // loop over all entity test basis functions
      // loop over all neighbor test basis functions
      for (size_t ii = 0; ii < rowsNe; ++ii) {
        auto& neighborNeighborRetRow = neighborNeighborRet[ii];
        const auto& neighborNeighborValsRow = neighborNeighborVals[ii];
        auto& neighborEntityRetRow = neighborEntityRet[ii];
        const auto& neighborEntityValsRow = neighborEntityVals[ii];
        // loop over all neighbor ansatz basis functions
        for (size_t jj = 0; jj < colsNe; ++jj) {
          neighborNeighborRetRow[jj] += neighborNeighborValsRow[jj] * integrationFactor * quadratureWeight;
        } // loop over all neighbor ansatz basis functions
        // loop over all entity ansatz basis functions
        for (size_t jj = 0; jj < colsEn; ++jj) {
          neighborEntityRetRow[jj] += neighborEntityValsRow[jj] * integrationFactor * quadratureWeight;
        } // loop over all entity ansatz basis functions
      } // loop over all neighbor test basis functions
    } // loop over all quadrature points
  } // void apply(...) const

private:
  const QuaternaryEvaluationType& evaluation() const
  {
    return static_cast< const QuaternaryEvaluationType& >(evaluation_);
  }

  const QuaternaryEvaluationImp evaluation_;
}; // class Codim1CouplingIntegral


// forward, to be used in the traits
template< class BinaryEvaluationImp >
class Codim1BoundaryIntegral;


template< class BinaryEvaluationImp >
class Codim1BoundaryIntegralTraits
{
  static_assert(std::is_base_of<  LocalEvaluation::Codim1Interface< typename BinaryEvaluationImp::Traits, 2 >,
                                  BinaryEvaluationImp >::value,
                "BinaryEvaluationImp has to be derived from LocalEvaluation::Codim1Interface< ..., 2 >!");
public:
  typedef Codim1BoundaryIntegral< BinaryEvaluationImp > derived_type;
  typedef LocalEvaluation::Codim1Interface< typename BinaryEvaluationImp::Traits, 2 > BinaryEvaluationType;
};


template< class BinaryEvaluationImp >
class Codim1BoundaryIntegral
  : public LocalOperator::Codim1BoundaryInterface< Codim1BoundaryIntegralTraits< BinaryEvaluationImp > >
{
public:
  typedef Codim1BoundaryIntegralTraits< BinaryEvaluationImp > Traits;
  typedef typename Traits::BinaryEvaluationType     BinaryEvaluationType;

private:
  static const size_t numTmpObjectsRequired_ = 1;

public:
  Codim1BoundaryIntegral(const BinaryEvaluationImp eval)
    : evaluation_(eval)
  {}

  template< class... Args >
  explicit Codim1BoundaryIntegral(Args&& ...args)
    : evaluation_(std::forward< Args >(args)...)
  {}

  size_t numTmpObjectsRequired() const
  {
    return numTmpObjectsRequired_;
  }

  template< class E, class IntersectionType, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  void apply(const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
             const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase,
             const IntersectionType& intersection,
             Dune::DynamicMatrix< R >& ret,
             std::vector< Dune::DynamicMatrix< R > >& tmpLocalMatrices) const
  {
    // local inducing function
    const auto& entity = testBase.entity();
    const auto localFunctions = evaluation_.localFunctions(entity);
    // quadrature
    typedef Dune::QuadratureRules< D, d - 1 > FaceQuadratureRules;
    typedef Dune::QuadratureRule< D, d - 1 > FaceQuadratureType;
    const size_t integrand_order = evaluation().order(localFunctions, testBase, ansatzBase);
    assert(integrand_order < std::numeric_limits< int >::max());
    const FaceQuadratureType& faceQuadrature = FaceQuadratureRules::rule(intersection.type(), int(integrand_order));
    // check matrix and tmp storage
    ret *= 0.0;
    const size_t rows = testBase.size();
    const size_t cols = ansatzBase.size();
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    assert(tmpLocalMatrices.size() >= numTmpObjectsRequired_);
    Dune::DynamicMatrix< R >& localMatrix = tmpLocalMatrices[0];
    // loop over all quadrature points
    for (auto quadPoint = faceQuadrature.begin(); quadPoint != faceQuadrature.end(); ++quadPoint) {
      const Dune::FieldVector< D, d - 1 > localPoint = quadPoint->position();
      const R integrationFactor = intersection.geometry().integrationElement(localPoint);
      const R quadratureWeight = quadPoint->weight();
      // evaluate local
      evaluation().evaluate(localFunctions, testBase, ansatzBase, intersection, localPoint, localMatrix);
      // compute integral
      assert(localMatrix.rows() >= rows);
      assert(localMatrix.cols() >= cols);
      // loop over all test basis functions
      for (size_t ii = 0; ii < rows; ++ii) {
        auto& retRow = ret[ii];
        const auto& localMatrixRow = localMatrix[ii];
        // loop over all ansatz basis functions
        for (size_t jj = 0; jj < cols; ++jj) {
          retRow[jj] += localMatrixRow[jj] * integrationFactor * quadratureWeight;
        } // loop over all ansatz basis functions
      } // loop over all test basis functions
    } // loop over all quadrature points
  } // void apply(...) const

private:
  const BinaryEvaluationType& evaluation() const
  {
    return static_cast< const BinaryEvaluationType& >(evaluation_);
  }

  const BinaryEvaluationImp evaluation_;
}; // class Codim1BoundaryIntegral


} // namespace LocalOperator
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALOPERATOR_CODIM1_HH
