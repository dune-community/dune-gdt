// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_EVALUATION_IPDG_FLUXES_HH
#define DUNE_GDT_EVALUATION_IPDG_FLUXES_HH

#include <tuple>

#include <dune/common/typetraits.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/densematrix.hh>

#include <dune/stuff/functions/interfaces.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace LocalEvaluation {
namespace IPDGFluxes {


// forward, to be used in the traits
template< class LocalizableFunctionImp >
class CouplingPrimal;


/**
 * \brief Traits for the IPDGFluxes::CouplingPrimal evaluation.
 */
template< class LocalizableFunctionImp >
class CouplingPrimalTraits
{
public:
  typedef CouplingPrimal< LocalizableFunctionImp >  derived_type;
  typedef LocalizableFunctionImp                    LocalizableFunctionType;
  dune_static_assert((Dune::IsBaseOf< Dune::Stuff::LocalizableFunction, LocalizableFunctionImp >::value),
                     "ERROR: LocalizableFunctionImp is not a Dune::Stuff::LocalizableFunction.");
};


template< class LocalizableFunctionImp >
class CouplingPrimal
  : public LocalEvaluation::Codim1Interface< CouplingPrimalTraits< LocalizableFunctionImp >, 4 >
{
public:
  typedef CouplingPrimalTraits< LocalizableFunctionImp > Traits;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;

  CouplingPrimal(const LocalizableFunctionType& inducingFunction, const double penaltyFactor = 1.0)
    : inducingFunction_(inducingFunction)
    , penaltyFactor_(penaltyFactor)
  {}

  template< class EntityType >
  std::tuple< typename LocalizableFunctionType::template LocalFunction< EntityType >::Type > localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(inducingFunction_.localFunction(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template< class LE, class LN, class TE, class AE, class TN, class AN,
            class D, int d, class R, int rT, int rCT, int rA, int rCA >
  static int order(const std::tuple< LE >& localFunctionsEntity,
                   const std::tuple< LN >& localFunctionsNeighbor,
                   const BaseFunctionSetInterface< TE, D, d, R, rT, rCT >& testBaseEntity,
                   const BaseFunctionSetInterface< AE, D, d, R, rA, rCA >& ansatzBaseEntity,
                   const BaseFunctionSetInterface< TN, D, d, R, rT, rCT >& testBaseNeighbor,
                   const BaseFunctionSetInterface< AN, D, d, R, rA, rCA >& ansatzBaseNeighbor)
  {
    const auto& localFunctionEntity = std::get< 0 >(localFunctionsEntity);
    const auto& localFunctionNeighbor = std::get< 0 >(localFunctionsNeighbor);
    return order(localFunctionEntity, localFunctionNeighbor,
                 testBaseEntity, ansatzBaseEntity,
                 testBaseNeighbor, ansatzBaseNeighbor);
  }

  template< class LE, class LN, class TE, class AE, class TN, class AN,
            class D, int d, class R, int rL, int rCL, int rT, int rCT, int rA, int rCA >
  static int order(const Dune::Stuff::LocalFunctionInterface< LE, D, d, R, rL, rCL >& localFunctionEntity,
                   const Dune::Stuff::LocalFunctionInterface< LN, D, d, R, rL, rCL >& localFunctionNeighbor,
                   const BaseFunctionSetInterface< TE, D, d, R, rT, rCT >& testBaseEntity,
                   const BaseFunctionSetInterface< AE, D, d, R, rA, rCA >& ansatzBaseEntity,
                   const BaseFunctionSetInterface< TN, D, d, R, rT, rCT >& testBaseNeighbor,
                   const BaseFunctionSetInterface< AN, D, d, R, rA, rCA >& ansatzBaseNeighbor)
  {
    if (localFunctionEntity.order() < 0 || localFunctionNeighbor.order() < 0)
      return -1;
    else
      return std::max(localFunctionEntity.order(), localFunctionNeighbor.order())
          + std::max(testBaseEntity.order(), testBaseNeighbor.order())
          + std::max(ansatzBaseEntity.order(), ansatzBaseNeighbor.order());
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template< class LE, class LN, class TE, class AE, class TN, class AN, class IntersectionType,
            class D, int d, class R, int rT, int rCT, int rA, int rCA >
  void evaluate(const std::tuple< LE >& localFunctionsEntity,
                const std::tuple< LN >& localFunctionsNeighbor,
                const BaseFunctionSetInterface< TE, D, d, R, rT, rCT >& testBaseEntity,
                const BaseFunctionSetInterface< AE, D, d, R, rA, rCA >& ansatzBaseEntity,
                const BaseFunctionSetInterface< TN, D, d, R, rT, rCT >& testBaseNeighbor,
                const BaseFunctionSetInterface< AN, D, d, R, rA, rCA >& ansatzBaseNeighbor,
                const IntersectionType& intersection,
                const Dune::FieldVector< D, d - 1 >& localPoint,
                Dune::DynamicMatrix< R >& entityEntityRet,
                Dune::DynamicMatrix< R >& neighborNeighborRet,
                Dune::DynamicMatrix< R >& entityNeighborRet,
                Dune::DynamicMatrix< R >& neighborEntityRet) const
  {
    const auto& localFunctionEntity = std::get< 0 >(localFunctionsEntity);
    const auto& localFunctionNeighbor = std::get< 0 >(localFunctionsNeighbor);
    evaluate(localFunctionEntity, localFunctionNeighbor,
             testBaseEntity, ansatzBaseEntity,
             testBaseNeighbor, ansatzBaseNeighbor,
             intersection, localPoint,
             entityEntityRet,
             neighborNeighborRet,
             entityNeighborRet,
             neighborEntityRet);
  }

  template< class LE, class LN, class TE, class AE, class TN, class AN, class IntersectionType,
            class D, int d, class R, int rL, int rCL, int rT, int rCT, int rA, int rCA >
  void evaluate(const Dune::Stuff::LocalFunctionInterface< LE, D, d, R, rL, rCL >& /*localFunctionEntity*/,
                const Dune::Stuff::LocalFunctionInterface< LN, D, d, R, rL, rCL >& /*localFunctionNeighbor*/,
                const BaseFunctionSetInterface< TE, D, d, R, rT, rCT >& /*testBaseEntity*/,
                const BaseFunctionSetInterface< AE, D, d, R, rA, rCA >& /*ansatzBaseEntity*/,
                const BaseFunctionSetInterface< TN, D, d, R, rT, rCT >& /*testBaseNeighbor*/,
                const BaseFunctionSetInterface< AN, D, d, R, rA, rCA >& /*ansatzBaseNeighbor*/,
                const IntersectionType& /*intersection*/,
                const Dune::FieldVector< D, d - 1 >& /*localPoint*/,
                Dune::DynamicMatrix< R >& /*entityEntityRet*/,
                Dune::DynamicMatrix< R >& /*neighborNeighborRet*/,
                Dune::DynamicMatrix< R >& /*entityNeighborRet*/,
                Dune::DynamicMatrix< R >& /*neighborEntityRet*/) const
  {
    dune_static_assert(Dune::AlwaysFalse< R >::value,
                       "ERROR: not implemented for this combination of dimensions!");
  }

  /**
   *  \brief  Computes the ipdg fluxes in a primal setting.
   *  \tparam LE        Traits of the entity Dune::Stuff::LocalFunctionInterface implementation
   *  \tparam LN        Traits of the neighbor Dune::Stuff::LocalFunctionInterface implementation
   *  \tparam TE        Traits of the entity test BaseFunctionSetInterface implementation
   *  \tparam AE        Traits of the entity ansatz BaseFunctionSetInterface implementation
   *  \tparam TN        Traits of the neighbor test BaseFunctionSetInterface implementation
   *  \tparam AN        Traits of the neighbor ansatz BaseFunctionSetInterface implementation
   *  \tparam IntersectionType Type of the codim 1 Intersection
   *  \tparam D         DomainFieldType
   *  \tparam d         dimDomain
   *  \tparam R         RangeFieldType
   */
  template< class LE, class LN, class TE, class AE, class TN, class AN, class IntersectionType, class D, int d, class R >
  void evaluate(const Dune::Stuff::LocalFunctionInterface< LE, D, d, R, 1, 1 >& localFunctionEntity,
                const Dune::Stuff::LocalFunctionInterface< LN, D, d, R, 1, 1 >& localFunctionNeighbor,
                const BaseFunctionSetInterface< TE, D, d, R, 1, 1 >& testBaseEntity,
                const BaseFunctionSetInterface< AE, D, d, R, 1, 1 >& ansatzBaseEntity,
                const BaseFunctionSetInterface< TN, D, d, R, 1, 1 >& testBaseNeighbor,
                const BaseFunctionSetInterface< AN, D, d, R, 1, 1 >& ansatzBaseNeighbor,
                const IntersectionType& intersection,
                const Dune::FieldVector< D, d - 1 >& localPoint,
                Dune::DynamicMatrix< R >& entityEntityRet,
                Dune::DynamicMatrix< R >& neighborNeighborRet,
                Dune::DynamicMatrix< R >& entityNeighborRet,
                Dune::DynamicMatrix< R >& neighborEntityRet) const
  {
    typedef Dune::FieldVector< D, d >                                                 DomainType;
    typedef typename BaseFunctionSetInterface< TE, D, d, R, 1, 1 >::RangeType         RangeType;
    typedef typename BaseFunctionSetInterface< TE, D, d, R, 1, 1 >::JacobianRangeType JacobianRangeType;
    // convert local point (which is in intersection coordinates) to entity/neighbor coordinates
    const DomainType localPointEn = intersection.geometryInInside().global(localPoint);
    const DomainType localPointNe = intersection.geometryInOutside().global(localPoint);
    const DomainType unitOuterNormal = intersection.unitOuterNormal(localPoint);
    // evaluate local function
    const RangeType functionValueEn = localFunctionEntity.evaluate(localPointEn);
    const RangeType functionValueNe = localFunctionNeighbor.evaluate(localPointNe);
    // evaluate bases
    // * entity
    // ** test
    const size_t rowsEn = testBaseEntity.size();
    std::vector< RangeType > testValuesEn(rowsEn, RangeType(0));
    std::vector< JacobianRangeType > testGradientsEn(rowsEn, JacobianRangeType(0));
    testBaseEntity.evaluate(localPointEn, testValuesEn);
    testBaseEntity.jacobian(localPointEn, testGradientsEn);
    // ** ansatz
    const size_t colsEn = ansatzBaseEntity.size();
    std::vector< RangeType > ansatzValuesEn(colsEn, RangeType(0));
    std::vector< JacobianRangeType > ansatzGradientsEn(colsEn, JacobianRangeType(0));
    ansatzBaseEntity.evaluate(localPointEn, ansatzValuesEn);
    ansatzBaseEntity.jacobian(localPointEn, ansatzGradientsEn);
    // * neighbor
    // ** test
    const size_t rowsNe = testBaseNeighbor.size();
    std::vector< RangeType > testValuesNe(rowsNe, RangeType(0));
    std::vector< JacobianRangeType > testGradientsNe(rowsNe, JacobianRangeType(0));
    testBaseNeighbor.evaluate(localPointNe, testValuesNe);
    testBaseNeighbor.jacobian(localPointNe, testGradientsNe);
    // ** ansatz
    const size_t colsNe = ansatzBaseNeighbor.size();
    std::vector< RangeType > ansatzValuesNe(colsNe, RangeType(0));
    std::vector< JacobianRangeType > ansatzGradientsNe(colsNe, JacobianRangeType(0));
    ansatzBaseNeighbor.evaluate(localPointNe, ansatzValuesNe);
    ansatzBaseNeighbor.jacobian(localPointNe, ansatzGradientsNe);
    // evaluate penalty parameter
    const R penaltyParameter = penaltyFactor_ / std::pow(intersection.geometry().volume(), 1.0);
    // compute products
    assert(entityEntityRet.rows() >= rowsEn);
    assert(entityEntityRet.cols() >= colsEn);
    assert(entityNeighborRet.rows() >= rowsEn);
    assert(entityNeighborRet.cols() >= colsNe);
    assert(neighborEntityRet.rows() >= rowsNe);
    assert(neighborEntityRet.cols() >= colsEn);
    assert(neighborNeighborRet.rows() >= rowsNe);
    assert(neighborNeighborRet.cols() >= colsNe);
    R gradientTimesNormal(0);
    R normalTimesGradient(0);
    R evaluationTimesEvaluation(0);
    RangeType rangeTmp(0);
    // loop over all entity test basis functions
    for (size_t ii = 0; ii < rowsEn; ++ii) {
      auto& entityEntityRetRow = entityEntityRet[ii];
      auto& entityNeighborRetRow = entityNeighborRet[ii];
      gradientTimesNormal = testGradientsEn[ii][0] * unitOuterNormal;
      // loop over all entity ansatz basis functions
      for (size_t jj = 0; jj < colsEn; ++jj) {
        rangeTmp = ansatzValuesEn[jj] * gradientTimesNormal;
        entityEntityRetRow[jj] += -0.5 * functionValueEn * rangeTmp;
        normalTimesGradient = unitOuterNormal * ansatzGradientsEn[jj][0];
        rangeTmp = testValuesEn[ii] * normalTimesGradient;
        entityEntityRetRow[jj] += -0.5 * functionValueEn * rangeTmp;
        evaluationTimesEvaluation = ansatzValuesEn[jj] * testValuesEn[ii];
        entityEntityRetRow[jj] += penaltyParameter * evaluationTimesEvaluation;
      } // loop over all entity ansatz basis functions
      // loop over all neighbor ansatz basis functions
      for (size_t jj = 0; jj < colsNe; ++jj) {
        rangeTmp = ansatzValuesNe[jj] * gradientTimesNormal;
        entityNeighborRetRow[jj] += 0.5 * functionValueEn * rangeTmp;
        normalTimesGradient = unitOuterNormal * ansatzGradientsNe[jj][0];
        rangeTmp = testValuesEn[ii] * normalTimesGradient;
        entityNeighborRetRow[jj] += -0.5 * functionValueNe * rangeTmp;
        evaluationTimesEvaluation = testValuesEn[ii] * ansatzValuesNe[jj];
        entityNeighborRetRow[jj] += -1.0 * penaltyParameter * evaluationTimesEvaluation;
      } // loop over all neighbor ansatz basis functions
    } // loop over all entity test basis functions
    // loop over all neighbor test basis functions
    for (size_t ii = 0; ii < rowsNe; ++ii) {
      auto& neighborEntityRetRow = neighborEntityRet[ii];
      auto& neighborNeighborRetRow = neighborNeighborRet[ii];
      gradientTimesNormal = testGradientsNe[ii][0] * unitOuterNormal;
      // loop over all entity ansatz basis functions
      for (size_t jj = 0; jj < colsEn; ++jj) {
        rangeTmp = ansatzValuesEn[jj] * gradientTimesNormal;
        neighborEntityRetRow[jj] += -0.5 * functionValueNe * rangeTmp;
        normalTimesGradient = unitOuterNormal * ansatzGradientsEn[jj][0];
        rangeTmp = testValuesNe[ii] * normalTimesGradient;
        neighborEntityRetRow[jj] += 0.5 * functionValueNe * rangeTmp;
        rangeTmp = testValuesNe[ii] * ansatzValuesEn[jj];
        neighborEntityRetRow[jj] += -1.0 * penaltyParameter * rangeTmp;
      } // loop over all entity ansatz basis functions
      // loop over all neighbor ansatz basis functions
      for (size_t jj = 0; jj < colsNe; ++jj) {
        rangeTmp = ansatzValuesNe[jj] * gradientTimesNormal;
        neighborNeighborRetRow[jj] += 0.5 * functionValueNe * rangeTmp;
        normalTimesGradient = unitOuterNormal * ansatzGradientsNe[jj][0];
        rangeTmp = testValuesNe[ii] * normalTimesGradient;
        neighborNeighborRetRow[jj] += 0.5 * functionValueNe * rangeTmp;
        evaluationTimesEvaluation = testValuesNe[ii] * ansatzValuesNe[jj];
        neighborNeighborRetRow[jj] += penaltyParameter * evaluationTimesEvaluation;
      } // loop over all neighbor ansatz basis functions
    } // loop over all neighbor test basis functions
  } // void evaluate< ..., 1, 1 >(...) const

private:
  const LocalizableFunctionType& inducingFunction_;
  const double penaltyFactor_;
}; // CouplingPrimal


// forward, to be used in the traits
template< class LocalizableFunctionImp >
class BoundaryDirichletLHS;


/**
 * \brief Traits for the IPDGFluxes::BoundaryDirichletLHS evaluation.
 */
template< class LocalizableFunctionImp >
class BoundaryDirichletLHSTraits
{
public:
  typedef BoundaryDirichletLHS< LocalizableFunctionImp >  derived_type;
  typedef LocalizableFunctionImp                          LocalizableFunctionType;
  dune_static_assert((Dune::IsBaseOf< Dune::Stuff::LocalizableFunction, LocalizableFunctionImp >::value),
                     "ERROR: LocalizableFunctionImp is not a Dune::Stuff::LocalizableFunction.");
};


template< class LocalizableFunctionImp >
class BoundaryDirichletLHS
  : public LocalEvaluation::Codim1Interface< BoundaryDirichletLHSTraits< LocalizableFunctionImp >, 2 >
{
public:
  typedef BoundaryDirichletLHSTraits< LocalizableFunctionImp > Traits;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;

  BoundaryDirichletLHS(const LocalizableFunctionType& inducingFunction, const double penaltyFactor)
    : inducingFunction_(inducingFunction)
    , penaltyFactor_(penaltyFactor)
  {}

  template< class EntityType >
  std::tuple< typename LocalizableFunctionType::template LocalFunction< EntityType >::Type > localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(inducingFunction_.localFunction(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template< class L, class T, class A, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  static int order(const std::tuple< L >& localFunctions,
                   const BaseFunctionSetInterface< T, D, d, R, rT, rCT >& testBase,
                   const BaseFunctionSetInterface< A, D, d, R, rA, rCA >& ansatzBase)
  {
    const auto& localFunction = std::get< 0 >(localFunctions);
    return order(localFunction, testBase, ansatzBase);
  }

  template< class L, class T, class A, class D, int d, class R, int rL, int rCL, int rT, int rCT, int rA, int rCA >
  static int order(const Dune::Stuff::LocalFunctionInterface< L, D, d, R, rL, rCL >& localFunction,
                   const BaseFunctionSetInterface< T, D, d, R, rT, rCT >& testBase,
                   const BaseFunctionSetInterface< A, D, d, R, rA, rCA >& ansatzBase)
  {
    if (localFunction.order() < 0)
      return -1;
    else
      return std::max(int(localFunction.order() + testBase.order() + ansatzBase.order()), 0);
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template< class L, class T, class A, class IntersectionType,
            class D, int d, class R, int rT, int rCT, int rA, int rCA >
  void evaluate(const std::tuple< L >& localFuncs,
                const BaseFunctionSetInterface< T, D, d, R, rT, rCT >& testBase,
                const BaseFunctionSetInterface< A, D, d, R, rA, rCA >& ansatzBase,
                const IntersectionType& intersection,
                const Dune::FieldVector< D, d - 1 >& localPoint,
                Dune::DynamicMatrix< R >& ret) const
  {
    const auto& localFunction = std::get< 0 >(localFuncs);
    evaluate(localFunction, testBase, ansatzBase, intersection, localPoint, ret);
  }

  template< class L, class T, class A, class IntersectionType,
            class D, int d, class R, int rL, int rCL, int rT, int rCT, int rA, int rCA >
  void evaluate(const Dune::Stuff::LocalFunctionInterface< L, D, d, R, rL, rCL >& /*localFunction*/,
                const BaseFunctionSetInterface< T, D, d, R, rT, rCT >& /*testBase*/,
                const BaseFunctionSetInterface< A, D, d, R, rA, rCA >& /*ansatzBase*/,
                const IntersectionType& /*intersection*/,
                const Dune::FieldVector< D, d - 1 >& /*localPoint*/,
                Dune::DynamicMatrix< R >& /*ret*/) const
  {
    dune_static_assert(Dune::AlwaysFalse< R >::value,
                       "ERROR: not implemented for this combination of dimensions!");
  }

  template< class L, class T, class A, class IntersectionType, class D, int d, class R >
  void evaluate(const Dune::Stuff::LocalFunctionInterface< L, D, d, R, 1, 1 >& localFunction,
                const BaseFunctionSetInterface< T, D, d, R, 1, 1 >& testBase,
                const BaseFunctionSetInterface< A, D, d, R, 1, 1 >& ansatzBase,
                const IntersectionType& intersection,
                const Dune::FieldVector< D, d - 1 >& localPoint,
                Dune::DynamicMatrix< R >& ret) const
  {
    typedef Dune::FieldVector< D, d >                                                 DomainType;
    typedef typename BaseFunctionSetInterface< T, D, d, R, 1, 1 >::RangeType          RangeType;
    typedef typename BaseFunctionSetInterface< T, D, d, R, 1, 1 >::JacobianRangeType  JacobianRangeType;
    // get local point (which is in intersection coordinates) in entity coordinates
    const DomainType localPointEntity = intersection.geometryInInside().global(localPoint);
    const DomainType unitOuterNormal = intersection.unitOuterNormal(localPoint);
    // evaluate local function
    const RangeType functionValue = localFunction.evaluate(localPointEntity);
    // evaluate bases
    // * test
    const size_t rows = testBase.size();
    std::vector< RangeType > testValues(rows, RangeType(0));
    std::vector< JacobianRangeType > testGradients(rows, JacobianRangeType(0));
    testBase.evaluate(localPointEntity, testValues);
    testBase.jacobian(localPointEntity, testGradients);
    // * ansatz
    const size_t cols = ansatzBase.size();
    std::vector< RangeType > ansatzValues(cols, RangeType(0));
    std::vector< JacobianRangeType > ansatzGradients(cols, JacobianRangeType(0));
    ansatzBase.evaluate(localPointEntity, ansatzValues);
    ansatzBase.jacobian(localPointEntity, ansatzGradients);
    // evaluate penalty parameter
    const R penaltyParameter = penaltyFactor_ / std::pow(intersection.geometry().volume(), 1.0);
    // compute products
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    R gradientTimesNormal(0);
    R normalTimesGradient(0);
    R evaluationTimesEvaluation(0);
    RangeType rangeTmp(0);
    // loop over all test basis functions
    for (size_t ii = 0; ii < rows; ++ii) {
      auto& retRow = ret[ii];
      gradientTimesNormal = testGradients[ii][0] * unitOuterNormal;
      // loop over all ansatz basis functions
      for (size_t jj = 0; jj < cols; ++jj) {
        rangeTmp = ansatzValues[ii] * gradientTimesNormal;
        retRow[jj] += -1.0 * functionValue * rangeTmp;
        normalTimesGradient = unitOuterNormal * ansatzGradients[jj][0];
        rangeTmp = testValues[ii] * normalTimesGradient;
        retRow[jj] += -1.0 * functionValue * rangeTmp;
        evaluationTimesEvaluation = ansatzValues[jj] * testValues[ii];
        retRow[jj] += penaltyParameter * evaluationTimesEvaluation;
      } // loop over all ansatz basis functions
    } // loop over all test basis functions
  } // void evaluate(...) const

private:
  const LocalizableFunctionType& inducingFunction_;
  const double penaltyFactor_;
}; // class BoundaryDirichletLHS


// forward, to be used in the traits
template< class LocalizableDiffusionFunctionImp, class LocalizableDirichletFunctionImp >
class BoundaryDirichletRHS;


/**
 * \brief Traits for the IPDGFluxes::BoundaryDirichletRHS evaluation.
 */
template< class LocalizableDiffusionFunctionImp, class LocalizableDirichletFunctionImp >
class BoundaryDirichletRHSTraits
{
public:
  typedef BoundaryDirichletRHS< LocalizableDiffusionFunctionImp, LocalizableDirichletFunctionImp > derived_type;
  typedef LocalizableDiffusionFunctionImp LocalizableDiffusionFunctionType;
  typedef LocalizableDirichletFunctionImp LocalizableDirichletFunctionType;
  dune_static_assert((Dune::IsBaseOf< Dune::Stuff::LocalizableFunction, LocalizableDiffusionFunctionImp >::value),
                     "ERROR: LocalizableDiffusionFunctionImp is not a Dune::Stuff::LocalizableFunction.");
  dune_static_assert((Dune::IsBaseOf< Dune::Stuff::LocalizableFunction, LocalizableDirichletFunctionImp >::value),
                     "ERROR: LocalizableDirichletFunctionImp is not a Dune::Stuff::LocalizableFunction.");
};


template< class LocalizableDiffusionFunctionImp, class LocalizableDirichletFunctionImp >
class BoundaryDirichletRHS
  : public LocalEvaluation::Codim1Interface< BoundaryDirichletRHSTraits< LocalizableDiffusionFunctionImp, LocalizableDirichletFunctionImp >, 1 >
{
public:
  typedef BoundaryDirichletRHSTraits< LocalizableDiffusionFunctionImp, LocalizableDirichletFunctionImp > Traits;
  typedef typename Traits::LocalizableDiffusionFunctionType LocalizableDiffusionFunctionType;
  typedef typename Traits::LocalizableDirichletFunctionType LocalizableDirichletFunctionType;

  BoundaryDirichletRHS(const LocalizableDiffusionFunctionType& diffusion,
                       const LocalizableDirichletFunctionType& dirichlet,
                       const double penaltyFactor)
    : diffusion_(diffusion)
    , dirichlet_(dirichlet)
    , penaltyFactor_(penaltyFactor)
  {}

  template< class EntityType >
  std::tuple< typename LocalizableDiffusionFunctionType::template LocalFunction< EntityType >::Type
            , typename LocalizableDirichletFunctionType::template LocalFunction< EntityType >::Type >
  localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_.localFunction(entity), dirichlet_.localFunction(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template< class LDF, class LDR, class T, class D, int d, class R, int r, int rC >
  static int order(const std::tuple< LDF, LDR >& localFunctions,
                   const BaseFunctionSetInterface< T, D, d, R, r, rC >& testBase)
  {
    const auto& localDiffusion = std::get< 0 >(localFunctions);
    const auto& localDirichlet = std::get< 1 >(localFunctions);
    return order(localDiffusion, localDirichlet, testBase);
  }

  template< class LF, class LR, class T,
            class D, int d, class R, int rLF, int rCLF, int rLR, int rCLR, int rT, int rCT >
  static int order(const Dune::Stuff::LocalFunctionInterface< LF, D, d, R, rLF, rCLF >& localDiffusion,
                   const Dune::Stuff::LocalFunctionInterface< LR, D, d, R, rLR, rCLR >& localDirichlet,
                   const BaseFunctionSetInterface< T, D, d, R, rT, rCT >& testBase)
  {
    if (localDiffusion.order() < 0 || localDirichlet.order() < 0)
      return -1;
    else {
      // all order()s are at least 0
      const size_t testOrder = testBase.order();
      const size_t testGradientOrder = std::max(int(testOrder - 1), 0);
      const size_t diffusionOrder = localDiffusion.order();
      const size_t dirichletOrder = localDirichlet.order();
      return std::max(testOrder + dirichletOrder,
                      diffusionOrder + testGradientOrder + dirichletOrder);
    }
  } // static int order(...)

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template< class LDF, class LDR, class T, class IntersectionType, class D, int d, class R, int r, int rC >
  void evaluate(const std::tuple< LDF, LDR >& localFuncs,
                const BaseFunctionSetInterface< T, D, d, R, r, rC >& testBase,
                const IntersectionType& intersection,
                const Dune::FieldVector< D, d - 1 >& localPoint,
                Dune::DynamicVector< R >& ret) const
  {
    const auto& localDiffusion = std::get< 0 >(localFuncs);
    const auto& localDirichlet = std::get< 1 >(localFuncs);
    evaluate(localDiffusion, localDirichlet, testBase, intersection, localPoint, ret);
  }

  template< class LDF, class LDR, class T, class IntersectionType,
            class D, int d, class R, int rLDF, int rCLDF, int rLDR, int rCLDR, int rT, int rCT >
  void evaluate(const Dune::Stuff::LocalFunctionInterface< LDF, D, d, R, rLDF, rCLDF >& /*localDiffusion*/,
                const Dune::Stuff::LocalFunctionInterface< LDR, D, d, R, rLDR, rCLDR >& /*localDirichlet*/,
                const BaseFunctionSetInterface< T, D, d, R, rT, rCT >& /*testBase*/,
                const IntersectionType& /*intersection*/,
                const Dune::FieldVector< D, d - 1 >& /*localPoint*/,
                Dune::DynamicVector< R >& /*ret*/) const
  {
    dune_static_assert(Dune::AlwaysFalse< R >::value,
                       "ERROR: not implemented for this combination of dimensions!");
  }

  template< class LDF, class LDR, class T, class IntersectionType, class D, int d, class R >
  void evaluate(const Dune::Stuff::LocalFunctionInterface< LDF, D, d, R, 1, 1 >& localDiffusion,
                const Dune::Stuff::LocalFunctionInterface< LDR, D, d, R, 1, 1 >& localDirichlet,
                const BaseFunctionSetInterface< T, D, d, R, 1, 1 >& testBase,
                const IntersectionType& intersection,
                const Dune::FieldVector< D, d - 1 >& localPoint,
                Dune::DynamicVector< R >& ret) const
  {
    typedef Dune::FieldVector< D, d >                                                 DomainType;
    typedef typename BaseFunctionSetInterface< T, D, d, R, 1, 1 >::RangeType          RangeType;
    typedef typename BaseFunctionSetInterface< T, D, d, R, 1, 1 >::JacobianRangeType  JacobianRangeType;
    // get local point (which is in intersection coordinates) in entity coordinates
    const DomainType localPointEntity = intersection.geometryInInside().global(localPoint);
    const DomainType unitOuterNormal = intersection.unitOuterNormal(localPoint);
    // evaluate local functions
    const RangeType diffusionValue = localDiffusion.evaluate(localPointEntity);
    const RangeType dirichletValue = localDirichlet.evaluate(localPointEntity);
    // evaluate basis
    const size_t size = testBase.size();
    std::vector< RangeType > testValues(size, RangeType(0));
    std::vector< JacobianRangeType > testGradients(size, JacobianRangeType(0));
    testBase.evaluate(localPointEntity, testValues);
    testBase.jacobian(localPointEntity, testGradients);
    // evaluate penalty parameter
    const R penaltyParameter = penaltyFactor_ / std::pow(intersection.geometry().volume(), 1.0);
    // compute
    assert(ret.size() >= size);
    R gradientTimesNormal(0);
    // loop over all test basis functions
    for (size_t ii = 0; ii < size; ++ii) {
      ret[ii] += penaltyParameter * (testValues[ii] * dirichletValue);
      gradientTimesNormal = testGradients[ii][0] * unitOuterNormal;
      ret[ii] += -1.0 * diffusionValue * gradientTimesNormal * dirichletValue;
    } // loop over all test basis functions
  } // void evaluate(...) const

private:
  const LocalizableDiffusionFunctionType& diffusion_;
  const LocalizableDirichletFunctionType& dirichlet_;
  const double penaltyFactor_;
}; // class BoundaryDirichletRHS


} // namespace IPDGFluxes
} // namespace LocalEvaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_EVALUATION_IPDG_FLUXES_HH
