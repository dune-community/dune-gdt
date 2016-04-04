// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LOCALOPERATOR_INTEGRALS_HH
#define DUNE_GDT_LOCALOPERATOR_INTEGRALS_HH

#include <type_traits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/densematrix.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/common/matrix.hh>

#include <dune/gdt/localevaluation/interface.hh>
#include <dune/gdt/type_traits.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forwards, to be used in the traits
template< class BinaryEvaluationType >
class LocalVolumeIntegralOperator;

template< class QuaternaryFaceIntegrandTypeType >
class LocalCouplingIntegralOperator;

template< class BinaryEvaluationType >
class LocalBoundaryIntegralOperator;


namespace internal {


template< class BinaryEvaluationType >
class LocalVolumeIntegralOperatorTraits
{
  static_assert(is_binary_volume_integrand< BinaryEvaluationType >::value,
                "BinaryEvaluationType has to be derived from LocalEvaluation::Codim0Interface< ..., 2 >!");
public:
  typedef LocalVolumeIntegralOperator< BinaryEvaluationType > derived_type;
};


template< class QuaternaryFaceIntegrandTypeType >
class LocalCouplingIntegralOperatorTraits
{
  static_assert(is_quaternary_face_integrand< QuaternaryFaceIntegrandTypeType >::value,
                "QuaternaryFaceIntegrandTypeType has to be derived from LocalEvaluation::Codim1Interface< ..., 4 >!");
public:
  typedef LocalCouplingIntegralOperator< QuaternaryFaceIntegrandTypeType > derived_type;
};


template< class BinaryEvaluationType >
class LocalBoundaryIntegralOperatorTraits
{
  static_assert(is_binary_face_integrand< BinaryEvaluationType >::value,
                "BinaryEvaluationType has to be derived from LocalEvaluation::Codim1Interface< ..., 2 >!");
public:
  typedef LocalBoundaryIntegralOperator< BinaryEvaluationType > derived_type;
};


} // namespace internal


template< class BinaryEvaluationType >
class LocalVolumeIntegralOperator
  : public LocalVolumeTwoFormInterface< internal::LocalVolumeIntegralOperatorTraits< BinaryEvaluationType > >
{
  typedef LocalVolumeIntegralOperator< BinaryEvaluationType >                 ThisType;
public:
  typedef internal::LocalVolumeIntegralOperatorTraits< BinaryEvaluationType > Traits;

  template< class... Args >
  explicit LocalVolumeIntegralOperator(Args&& ...args)
    : integrand_(std::forward< Args >(args)...)
    , over_integrate_(0)
  {}

  template< class... Args >
  explicit LocalVolumeIntegralOperator(const int over_integrate, Args&& ...args)
    : integrand_(std::forward< Args >(args)...)
    , over_integrate_(boost::numeric_cast< size_t >(over_integrate))
  {}

  template< class... Args >
  explicit LocalVolumeIntegralOperator(const size_t over_integrate, Args&& ...args)
    : integrand_(std::forward< Args >(args)...)
    , over_integrate_(over_integrate)
  {}

  LocalVolumeIntegralOperator(const ThisType& other) = default;
  LocalVolumeIntegralOperator(ThisType&& source) = default;

  template< class E, class D, size_t d, class R, size_t rT, size_t rCT, size_t rA, size_t rCA >
  void apply2(const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& test_base,
              const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatz_base,
              Dune::DynamicMatrix< R >& ret) const
  {
    const auto& entity = ansatz_base.entity();
    const auto local_functions = integrand_.localFunctions(entity);
    // create quadrature
    const size_t integrand_order = integrand_.order(local_functions, ansatz_base, test_base) + over_integrate_;
    const auto& quadrature = QuadratureRules< D, d >::rule(entity.type(),
                                                           boost::numeric_cast< int >(integrand_order));
    // prepare storage
    const size_t rows = test_base.size();
    const size_t cols = ansatz_base.size();
    ret *= 0.0;
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    Dune::DynamicMatrix< R > evaluation_result(rows, cols, 0.);  // \todo: make mutable member, after SMP refactor
    // loop over all quadrature points
    for (const auto& quadrature_point : quadrature) {
      const auto xx = quadrature_point.position();
      // integration factors
      const auto integration_factor = entity.geometry().integrationElement(xx);
      const auto quadrature_weight = quadrature_point.weight();
      // evaluate the integrand
      integrand_.evaluate(local_functions, test_base, ansatz_base, xx, evaluation_result);
      // compute integral
      for (size_t ii = 0; ii < rows; ++ii) {
        const auto& evaluation_result_row = evaluation_result[ii];
        auto& ret_row = ret[ii];
        for (size_t jj = 0; jj < cols; ++jj)
          ret_row[jj] += evaluation_result_row[jj] * integration_factor * quadrature_weight;
      } // compute integral
    } // loop over all quadrature points
  } // ... apply2(...)

private:
  const BinaryEvaluationType integrand_;
  const size_t over_integrate_;
}; // class LocalVolumeIntegralOperator


#if 0
template< class QuaternaryEvaluationType >
class LocalCouplingIntegralOperator
    : public LocalOperator::Codim1CouplingInterface< internal::LocalCouplingIntegralOperatorTraits< QuaternaryEvaluationType > >
{
public:
  typedef internal::LocalCouplingIntegralOperatorTraits< QuaternaryEvaluationType > Traits;

  template< class... Args >
  explicit LocalCouplingIntegralOperator(Args&& ...args)
    : evaluation_(std::forward< Args >(args)...)
    , over_integrate_(0)
  {}

  template< class... Args >
  explicit LocalCouplingIntegralOperator(const int over_integrate, Args&& ...args)
    : evaluation_(std::forward< Args >(args)...)
    , over_integrate_(over_integrate)
  {}

  template< class... Args >
  explicit LocalCouplingIntegralOperator(const size_t over_integrate, Args&& ...args)
    : evaluation_(std::forward< Args >(args)...)
    , over_integrate_(boost::numeric_cast< size_t >(over_integrate))
  {}

  template< class E, class N, class IntersectionType, class D, size_t d, class R, size_t rT, size_t rCT, size_t rA, size_t rCA >
  void apply2(const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& entityTestBase,
             const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& entityAnsatzBase,
             const Stuff::LocalfunctionSetInterface< N, D, d, R, rT, rCT >& neighborTestBase,
             const Stuff::LocalfunctionSetInterface< N, D, d, R, rA, rCA >& neighborAnsatzBase,
             const IntersectionType& intersection,
             Dune::DynamicMatrix< R >& entityEntityRet,
             Dune::DynamicMatrix< R >& neighborNeighborRet,
             Dune::DynamicMatrix< R >& entityNeighborRet,
             Dune::DynamicMatrix< R >& neighborEntityRet) const
  {
    // local inducing function
    const auto& entity = entityTestBase.entity();
    const auto localFunctionsEn = evaluation_.localFunctions(entity);
    const auto& neighbor = neighborTestBase.entity();
    const auto localFunctionsNe = evaluation_.localFunctions(neighbor);
    // quadrature
    const size_t integrand_order = evaluation_.order(localFunctionsEn, localFunctionsNe,
                                                     entityTestBase, entityAnsatzBase,
                                                     neighborTestBase, neighborAnsatzBase) + over_integrate_;
    const auto& faceQuadrature = QuadratureRules< D, d - 1 >::rule(intersection.type(),
                                                                   boost::numeric_cast< int >(integrand_order));
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
    auto& entityEntityVals = tmpLocalMatrices[0];
    auto& neighborNeighborVals = tmpLocalMatrices[1];
    auto& entityNeighborVals = tmpLocalMatrices[2];
    auto& neighborEntityVals = tmpLocalMatrices[3];
    // loop over all quadrature points
    for (auto quadPoint = faceQuadrature.begin(); quadPoint != faceQuadrature.end(); ++quadPoint) {
      const Dune::FieldVector< D, d - 1 > localPoint = quadPoint->position();
      const auto integrationFactor = intersection.geometry().integrationElement(localPoint);
      const auto quadratureWeight = quadPoint->weight();
      // evaluate local
      evaluation_.evaluate(localFunctionsEn, localFunctionsNe,
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
        auto& entityEntityret_row = entityEntityRet[ii];
        const auto& entityEntityValsRow = entityEntityVals[ii];
        auto& entityNeighborret_row = entityNeighborRet[ii];
        const auto& entityNeighborValsRow = entityNeighborVals[ii];
        // loop over all entity ansatz basis functions
        for (size_t jj = 0; jj < colsEn; ++jj) {
          entityEntityret_row[jj] += entityEntityValsRow[jj] * integrationFactor * quadratureWeight;
        } // loop over all entity ansatz basis functions
        // loop over all neighbor ansatz basis functions
        for (size_t jj = 0; jj < colsNe; ++jj) {
          entityNeighborret_row[jj] += entityNeighborValsRow[jj] * integrationFactor * quadratureWeight;
        } // loop over all neighbor ansatz basis functions
      } // loop over all entity test basis functions
      // loop over all neighbor test basis functions
      for (size_t ii = 0; ii < rowsNe; ++ii) {
        auto& neighborNeighborret_row = neighborNeighborRet[ii];
        const auto& neighborNeighborValsRow = neighborNeighborVals[ii];
        auto& neighborEntityret_row = neighborEntityRet[ii];
        const auto& neighborEntityValsRow = neighborEntityVals[ii];
        // loop over all neighbor ansatz basis functions
        for (size_t jj = 0; jj < colsNe; ++jj) {
          neighborNeighborret_row[jj] += neighborNeighborValsRow[jj] * integrationFactor * quadratureWeight;
        } // loop over all neighbor ansatz basis functions
        // loop over all entity ansatz basis functions
        for (size_t jj = 0; jj < colsEn; ++jj) {
          neighborEntityret_row[jj] += neighborEntityValsRow[jj] * integrationFactor * quadratureWeight;
        } // loop over all entity ansatz basis functions
      } // loop over all neighbor test basis functions
    } // loop over all quadrature points
  } // void apply(...) const

private:
  const QuaternaryEvaluationType evaluation_;
  const size_t over_integrate_;
}; // class LocalCouplingIntegralOperator


template< class BinaryEvaluationType >
class LocalBoundaryIntegralOperator
  : public LocalOperator::Codim1BoundaryInterface< internal::LocalBoundaryIntegralOperatorTraits< BinaryEvaluationType > >
{
public:
  typedef internal::LocalBoundaryIntegralOperatorTraits< BinaryEvaluationType > Traits;

  template< class... Args >
  LocalBoundaryIntegralOperator(Args&& ...args)
    : evaluation_(std::forward< Args >(args)...)
    , over_integrate_(0)
  {}

  template< class... Args >
  LocalBoundaryIntegralOperator(const int over_integrate, Args&& ...args)
    : evaluation_(std::forward< Args >(args)...)
    , over_integrate_(over_integrate)
  {}

  template< class... Args >
  LocalBoundaryIntegralOperator(const size_t over_integrate, Args&& ...args)
    : evaluation_(std::forward< Args >(args)...)
    , over_integrate_(boost::numeric_cast< size_t >(over_integrate))
  {}

  template< class E, class IntersectionType, class D, size_t d, class R, size_t rT, size_t rCT, size_t rA, size_t rCA >
  void apply2(const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
             const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase,
             const IntersectionType& intersection,
             Dune::DynamicMatrix< R >& ret) const
  {
    // local inducing function
    const auto& entity = testBase.entity();
    const auto localFunctions = evaluation_.localFunctions(entity);
    // quadrature
    typedef Dune::QuadratureRules< D, d - 1 > FaceQuadratureRules;
    typedef Dune::QuadratureRule< D, d - 1 > FaceQuadratureType;
    const auto integrand_order = evaluation_.order(localFunctions, testBase, ansatzBase) + over_integrate_;
    const FaceQuadratureType& faceQuadrature = FaceQuadratureRules::rule(intersection.type(),
                                                                         boost::numeric_cast< int >(integrand_order));
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
      evaluation_.evaluate(localFunctions, testBase, ansatzBase, intersection, localPoint, localMatrix);
      // compute integral
      assert(localMatrix.rows() >= rows);
      assert(localMatrix.cols() >= cols);
      // loop over all test basis functions
      for (size_t ii = 0; ii < rows; ++ii) {
        auto& ret_row = ret[ii];
        const auto& localMatrixRow = localMatrix[ii];
        // loop over all ansatz basis functions
        for (size_t jj = 0; jj < cols; ++jj) {
          ret_row[jj] += localMatrixRow[jj] * integrationFactor * quadratureWeight;
        } // loop over all ansatz basis functions
      } // loop over all test basis functions
    } // loop over all quadrature points
  } // void apply(...) const

private:
  const BinaryEvaluationType evaluation_;
  const size_t over_integrate_;
}; // class LocalBoundaryIntegralOperator
#endif // 0


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALOPERATOR_INTEGRALS_HH
