// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LOCALEVALUATION_SIPDG_HH
#define DUNE_GDT_LOCALEVALUATION_SIPDG_HH

#include <tuple>
#include <type_traits>

#include <dune/stuff/common/disable_warnings.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/common/color.hh>
#include <dune/stuff/functions/interfaces.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace LocalEvaluation {

/**
 *  \brief      Contains local evaluations for the symmetric interior penalty discontinuous Galerkin (SIPDG)
 *              discretization.
 *
 *              For the choice of penalization and the role of the user input see Epshteyn, Riviere (2007):
 *              "Estimation of penalty parameters for symmetric interior penalty Galerkin methods"
 */
namespace SIPDG {


template <class LocalizableFunctionImp>
class Inner;


template <class LocalizableFunctionImp>
class InnerTraits
{
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, LocalizableFunctionImp>::value,
                "LocalizableFunctionImp has to be tagged as Stuff::IsLocalizableFunction!");

public:
  typedef Inner<LocalizableFunctionImp> derived_type;
  typedef LocalizableFunctionImp LocalizableFunctionType;
  typedef typename LocalizableFunctionType::EntityType EntityType;
  typedef typename LocalizableFunctionType::DomainFieldType DomainFieldType;
  typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;
  typedef std::tuple<std::shared_ptr<LocalfunctionType>> LocalfunctionTupleType;
  static const unsigned int dimDomain = LocalizableFunctionType::dimDomain;
};


/**
 *  see Epshteyn, Riviere, 2007 for the meaning of beta
 */
template <class LocalizableFunctionImp>
class Inner : public LocalEvaluation::Codim1Interface<InnerTraits<LocalizableFunctionImp>, 4>
{
public:
  typedef InnerTraits<LocalizableFunctionImp> Traits;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = Traits::dimDomain;

  Inner(const LocalizableFunctionType& inducingFunction, const double beta = 1.0 / (dimDomain - 1.0))
    : inducingFunction_(inducingFunction)
    , beta_(beta)
  {
  }

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(inducingFunction_.local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class R, int rT, int rCT, int rA, int rCA>
  static size_t
  order(const LocalfunctionTupleType& localFunctionsEntity, const LocalfunctionTupleType& localFunctionsNeighbor,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBaseEntity,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBaseEntity,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBaseNeighbor,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBaseNeighbor)
  {
    const auto localFunctionEntity   = std::get<0>(localFunctionsEntity);
    const auto localFunctionNeighbor = std::get<0>(localFunctionsNeighbor);
    return order(*localFunctionEntity,
                 *localFunctionNeighbor,
                 testBaseEntity,
                 ansatzBaseEntity,
                 testBaseNeighbor,
                 ansatzBaseNeighbor);
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class IntersectionType, class R, int rT, int rCT, int rA, int rCA>
  void evaluate(
      const LocalfunctionTupleType& localFunctionsEntity, const LocalfunctionTupleType& localFunctionsNeighbor,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBaseEntity,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBaseEntity,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBaseNeighbor,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBaseNeighbor,
      const IntersectionType& intersection, const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint,
      Dune::DynamicMatrix<R>& entityEntityRet, Dune::DynamicMatrix<R>& neighborNeighborRet,
      Dune::DynamicMatrix<R>& entityNeighborRet, Dune::DynamicMatrix<R>& neighborEntityRet) const
  {
    const auto localFunctionEntity   = std::get<0>(localFunctionsEntity);
    const auto localFunctionNeighbor = std::get<0>(localFunctionsNeighbor);
    evaluate(*localFunctionEntity,
             *localFunctionNeighbor,
             testBaseEntity,
             ansatzBaseEntity,
             testBaseNeighbor,
             ansatzBaseNeighbor,
             intersection,
             localPoint,
             entityEntityRet,
             neighborNeighborRet,
             entityNeighborRet,
             neighborEntityRet);
  }

private:
  template <class R, int rL, int rCL, int rT, int rCT, int rA, int rCA>
  static size_t
  order(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& localFunctionEntity,
        const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& localFunctionNeighbor,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBaseEntity,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBaseEntity,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBaseNeighbor,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBaseNeighbor)
  {
    return std::max(localFunctionEntity.order(), localFunctionNeighbor.order())
           + std::max(testBaseEntity.order(), testBaseNeighbor.order())
           + std::max(ansatzBaseEntity.order(), ansatzBaseNeighbor.order());
  }

  template <class IntersectionType, class R, int rL, int rCL, int rT, int rCT, int rA, int rCA>
  void evaluate(
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& /*localFunctionEntity*/,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL,
                                          rCL>& /*localFunctionNeighbor*/,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& /*testBaseEntity*/,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& /*ansatzBaseEntity*/,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& /*testBaseNeighbor*/,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA,
                                             rCA>& /*ansatzBaseNeighbor*/,
      const IntersectionType& /*intersection*/, const Dune::FieldVector<DomainFieldType, dimDomain - 1>& /*localPoint*/,
      Dune::DynamicMatrix<R>& /*entityEntityRet*/, Dune::DynamicMatrix<R>& /*neighborNeighborRet*/,
      Dune::DynamicMatrix<R>& /*entityNeighborRet*/, Dune::DynamicMatrix<R>& /*neighborEntityRet*/) const
  {
    static_assert(Dune::AlwaysFalse<R>::value, "Not implemented for these dimensions!");
  }

  /**
   *  \brief  Computes the ipdg fluxes in a primal setting.
   *  \tparam IntersectionType Type of the codim 1 Intersection
   *  \tparam R         RangeFieldType
   */
  template <class IntersectionType, class R>
  void evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, 2, R, 1, 1>& localFunctionEntity,
                const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, 2, R, 1, 1>& localFunctionNeighbor,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, 2, R, 1, 1>& testBaseEntity,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, 2, R, 1, 1>& ansatzBaseEntity,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, 2, R, 1, 1>& testBaseNeighbor,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, 2, R, 1, 1>& ansatzBaseNeighbor,
                const IntersectionType& intersection, const Dune::FieldVector<DomainFieldType, 1>& localPoint,
                Dune::DynamicMatrix<R>& entityEntityRet, Dune::DynamicMatrix<R>& neighborNeighborRet,
                Dune::DynamicMatrix<R>& entityNeighborRet, Dune::DynamicMatrix<R>& neighborEntityRet) const
  {
    // clear ret
    entityEntityRet *= 0.0;
    neighborNeighborRet *= 0.0;
    entityNeighborRet *= 0.0;
    neighborEntityRet *= 0.0;
    typedef typename Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, 2, R, 1, 1>::DomainType DomainType;
    typedef typename Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, 2, R, 1, 1>::RangeType RangeType;
    typedef typename Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, 2, R, 1, 1>::JacobianRangeType
        JacobianRangeType;
    // convert local point (which is in intersection coordinates) to entity/neighbor coordinates
    const DomainType localPointEn    = intersection.geometryInInside().global(localPoint);
    const DomainType localPointNe    = intersection.geometryInOutside().global(localPoint);
    const DomainType unitOuterNormal = intersection.unitOuterNormal(localPoint);
    // evaluate local function
    const RangeType functionValueEn = localFunctionEntity.evaluate(localPointEn);
    const RangeType functionValueNe = localFunctionNeighbor.evaluate(localPointNe);
    // compute penalty (see Epshteyn, Riviere, 2007)
    const size_t max_polorder =
        std::max(testBaseEntity.order(),
                 std::max(ansatzBaseEntity.order(), std::max(testBaseNeighbor.order(), ansatzBaseNeighbor.order())));
    R sigma = 1.0;
    if (max_polorder <= 1)
      sigma *= 8.0;
    else if (max_polorder <= 2)
      sigma *= 20.0;
    else if (max_polorder <= 3)
      sigma *= 38.0;
    else {
#ifndef NDEBUG
#ifndef DUNE_GDT_LOCALEVALUATION_SIPDG_DISABLE_WARNINGS
      std::cout << "\n" << Dune::Stuff::Common::colorString("WARNING(dune.gdt.localevaluation.sipdg.inner):")
                << " a polynomial order of " << max_polorder << " is untested!" << std::endl
                << "  (#define DUNE_GDT_LOCALEVALUATION_SIPDG_DISABLE_WARNINGS to disable this warning)!" << std::endl;
#endif
#endif
      sigma *= 50.0;
    }
    const R penalty = sigma / std::pow(intersection.geometry().volume(), beta_);
    // evaluate bases
    // * entity
    //   * test
    const size_t rowsEn = testBaseEntity.size();
    std::vector<RangeType> testValuesEn(rowsEn, RangeType(0));
    std::vector<JacobianRangeType> testGradientsEn(rowsEn, JacobianRangeType(0));
    testBaseEntity.evaluate(localPointEn, testValuesEn);
    testBaseEntity.jacobian(localPointEn, testGradientsEn);
    //   * ansatz
    const size_t colsEn = ansatzBaseEntity.size();
    std::vector<RangeType> ansatzValuesEn(colsEn, RangeType(0));
    std::vector<JacobianRangeType> ansatzGradientsEn(colsEn, JacobianRangeType(0));
    ansatzBaseEntity.evaluate(localPointEn, ansatzValuesEn);
    ansatzBaseEntity.jacobian(localPointEn, ansatzGradientsEn);
    // * neighbor
    //   * test
    const size_t rowsNe = testBaseNeighbor.size();
    std::vector<RangeType> testValuesNe(rowsNe, RangeType(0));
    std::vector<JacobianRangeType> testGradientsNe(rowsNe, JacobianRangeType(0));
    testBaseNeighbor.evaluate(localPointNe, testValuesNe);
    testBaseNeighbor.jacobian(localPointNe, testGradientsNe);
    //   * ansatz
    const size_t colsNe = ansatzBaseNeighbor.size();
    std::vector<RangeType> ansatzValuesNe(colsNe, RangeType(0));
    std::vector<JacobianRangeType> ansatzGradientsNe(colsNe, JacobianRangeType(0));
    ansatzBaseNeighbor.evaluate(localPointNe, ansatzValuesNe);
    ansatzBaseNeighbor.jacobian(localPointNe, ansatzGradientsNe);
    // compute the evaluations
    assert(entityEntityRet.rows() >= rowsEn);
    assert(entityEntityRet.cols() >= colsEn);
    assert(entityNeighborRet.rows() >= rowsEn);
    assert(entityNeighborRet.cols() >= colsNe);
    assert(neighborEntityRet.rows() >= rowsNe);
    assert(neighborEntityRet.cols() >= colsEn);
    assert(neighborNeighborRet.rows() >= rowsNe);
    assert(neighborNeighborRet.cols() >= colsNe);
    // loop over all entity test basis functions
    for (size_t ii = 0; ii < rowsEn; ++ii) {
      auto& entityEntityRetRow   = entityEntityRet[ii];
      auto& entityNeighborRetRow = entityNeighborRet[ii];
      // loop over all entity ansatz basis functions
      for (size_t jj = 0; jj < colsEn; ++jj) {
        // consistency term
        entityEntityRetRow[jj] +=
            -0.5 * functionValueEn * (ansatzGradientsEn[jj][0] * unitOuterNormal) * testValuesEn[ii];
        // symmetry term
        entityEntityRetRow[jj] +=
            -0.5 * ansatzValuesEn[jj] * functionValueEn * (testGradientsEn[ii][0] * unitOuterNormal);
        // penalty term
        entityEntityRetRow[jj] += penalty * ansatzValuesEn[jj] * testValuesEn[ii];
      } // loop over all entity ansatz basis functions
      // loop over all neighbor ansatz basis functions
      for (size_t jj = 0; jj < colsNe; ++jj) {
        // consistency term
        entityNeighborRetRow[jj] +=
            -0.5 * functionValueNe * (ansatzGradientsNe[jj][0] * unitOuterNormal) * testValuesEn[ii];
        // symmetry term
        entityNeighborRetRow[jj] +=
            0.5 * ansatzValuesNe[jj] * functionValueEn * (testGradientsEn[ii][0] * unitOuterNormal);
        // penalty term
        entityNeighborRetRow[jj] += -1.0 * penalty * ansatzValuesNe[jj] * testValuesEn[ii];
      } // loop over all neighbor ansatz basis functions
    } // loop over all entity test basis functions
    // loop over all neighbor test basis functions
    for (size_t ii = 0; ii < rowsNe; ++ii) {
      auto& neighborEntityRetRow   = neighborEntityRet[ii];
      auto& neighborNeighborRetRow = neighborNeighborRet[ii];
      // loop over all entity ansatz basis functions
      for (size_t jj = 0; jj < colsEn; ++jj) {
        // consistency term
        neighborEntityRetRow[jj] +=
            0.5 * functionValueEn * (ansatzGradientsEn[jj][0] * unitOuterNormal) * testValuesNe[ii];
        // symmetry term
        neighborEntityRetRow[jj] +=
            -0.5 * ansatzValuesEn[jj] * functionValueNe * (testGradientsNe[ii][0] * unitOuterNormal);
        // penalty term
        neighborEntityRetRow[jj] += -1.0 * penalty * ansatzValuesEn[jj] * testValuesNe[ii];
      } // loop over all entity ansatz basis functions
      // loop over all neighbor ansatz basis functions
      for (size_t jj = 0; jj < colsNe; ++jj) {
        // consistency term
        neighborNeighborRetRow[jj] +=
            0.5 * functionValueNe * (ansatzGradientsNe[jj][0] * unitOuterNormal) * testValuesNe[ii];
        // symmetry term
        neighborNeighborRetRow[jj] +=
            0.5 * ansatzValuesNe[jj] * functionValueNe * (testGradientsNe[ii][0] * unitOuterNormal);
        // penalty term
        neighborNeighborRetRow[jj] += penalty * ansatzValuesNe[jj] * testValuesNe[ii];
      } // loop over all neighbor ansatz basis functions
    } // loop over all neighbor test basis functions
  } // void evaluate< ..., 1, 1 >(...) const
  const LocalizableFunctionType& inducingFunction_;
  const double beta_;
}; // CouplingPrimal


// forward
template <class LocalizableFunctionImp>
class BoundaryLHS;


template <class LocalizableFunctionImp>
class BoundaryLHSTraits
{
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, LocalizableFunctionImp>::value,
                "LocalizableFunctionImp has to be tagged as Stuff::IsLocalizableFunction!");

public:
  typedef LocalizableFunctionImp LocalizableFunctionType;
  typedef BoundaryLHS<LocalizableFunctionImp> derived_type;
  typedef typename LocalizableFunctionType::EntityType EntityType;
  typedef typename LocalizableFunctionType::DomainFieldType DomainFieldType;
  typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;
  typedef std::tuple<std::shared_ptr<LocalfunctionType>> LocalfunctionTupleType;
  static const unsigned int dimDomain = LocalizableFunctionType::dimDomain;
};


template <class LocalizableFunctionImp>
class BoundaryLHS : public LocalEvaluation::Codim1Interface<BoundaryLHSTraits<LocalizableFunctionImp>, 2>
{
public:
  typedef BoundaryLHSTraits<LocalizableFunctionImp> Traits;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = Traits::dimDomain;

  BoundaryLHS(const LocalizableFunctionType& inducingFunction, const double beta = 1.0 / (dimDomain - 1.0))
    : inducingFunction_(inducingFunction)
    , beta_(beta)
  {
  }

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(inducingFunction_.local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class R, int rT, int rCT, int rA, int rCA>
  static size_t
  order(const LocalfunctionTupleType localFuncs,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase)
  {
    const auto localFunction = std::get<0>(localFuncs);
    return redirect_order(*localFunction, testBase, ansatzBase);
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class IntersectionType, class R, int rT, int rCT, int rA, int rCA>
  void evaluate(const LocalfunctionTupleType localFuncs,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase,
                const IntersectionType& intersection,
                const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    const auto localFunction = std::get<0>(localFuncs);
    redirect_evaluate(*localFunction, testBase, ansatzBase, intersection, localPoint, ret);
  }

private:
  /**
   * \return localFunction.order() + testBase.order() + ansatzBase.order();
   */
  template <class R, int rL, int rCL, int rT, int rCT, int rA, int rCA>
  static size_t
  redirect_order(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& localFunction,
                 const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
                 const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase)
  {
    return localFunction.order() + testBase.order() + ansatzBase.order();
  }

  template <class IntersectionType, class R, int rL, int rCL, int rT, int rCT, int rA, int rCA>
  void redirect_evaluate(
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& /*localFunction*/,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& /*testBase*/,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& /*ansatzBase*/,
      const IntersectionType& /*intersection*/, const Dune::FieldVector<DomainFieldType, dimDomain - 1>& /*localPoint*/,
      Dune::DynamicMatrix<R>& /*ret*/) const
  {
    static_assert(Dune::AlwaysFalse<R>::value, "Not implemented for these dimensions!");
  }

  template <class IntersectionType, class R>
  void redirect_evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, 2, R, 1, 1>& localFunction,
                         const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, 2, R, 1, 1>& testBase,
                         const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, 2, R, 1, 1>& ansatzBase,
                         const IntersectionType& intersection, const Dune::FieldVector<DomainFieldType, 1>& localPoint,
                         Dune::DynamicMatrix<R>& ret) const
  {
    // clear ret
    ret *= 0.0;
    typedef typename Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, 2, R, 1, 1>::DomainType DomainType;
    typedef typename Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, 2, R, 1, 1>::RangeType RangeType;
    typedef typename Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, 2, R, 1, 1>::JacobianRangeType
        JacobianRangeType;
    // get local point (which is in intersection coordinates) in entity coordinates
    const DomainType localPointEntity = intersection.geometryInInside().global(localPoint);
    const DomainType unitOuterNormal  = intersection.unitOuterNormal(localPoint);
    // evaluate local function
    const RangeType functionValue = localFunction.evaluate(localPointEntity);
    // compute penalty (see Epshteyn, Riviere, 2007)
    const size_t max_polorder = std::max(testBase.order(), ansatzBase.order());
    R sigma = 1.0;
    if (max_polorder <= 1)
      sigma *= 14.0;
    else if (max_polorder <= 2)
      sigma *= 38.0;
    else if (max_polorder <= 3)
      sigma *= 74.0;
    else {
#ifndef NDEBUG
#ifndef DUNE_GDT_LOCALEVALUATION_SIPDG_DISABLE_WARNINGS
      std::cout << "\n" << Dune::Stuff::Common::colorString("WARNING(dune.gdt.localevaluation.sipdg.boundarylhs):")
                << " a polynomial order of " << max_polorder << " is untested!" << std::endl
                << "  (#define DUNE_GDT_LOCALEVALUATION_SIPDG_DISABLE_WARNINGS to disable this warning)!" << std::endl;
#endif
#endif
      sigma *= 100.0;
    }
    const R penalty = sigma / std::pow(intersection.geometry().volume(), beta_);
    // evaluate bases
    // * test
    const size_t rows = testBase.size();
    std::vector<RangeType> testValues(rows, RangeType(0));
    std::vector<JacobianRangeType> testGradients(rows, JacobianRangeType(0));
    testBase.evaluate(localPointEntity, testValues);
    testBase.jacobian(localPointEntity, testGradients);
    // * ansatz
    const size_t cols = ansatzBase.size();
    std::vector<RangeType> ansatzValues(cols, RangeType(0));
    std::vector<JacobianRangeType> ansatzGradients(cols, JacobianRangeType(0));
    ansatzBase.evaluate(localPointEntity, ansatzValues);
    ansatzBase.jacobian(localPointEntity, ansatzGradients);
    // compute products
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    // loop over all test basis functions
    for (size_t ii = 0; ii < rows; ++ii) {
      auto& retRow = ret[ii];
      // loop over all ansatz basis functions
      for (size_t jj = 0; jj < cols; ++jj) {
        // consistency term
        retRow[jj] += -1.0 * functionValue * (ansatzGradients[jj][0] * unitOuterNormal) * testValues[ii];
        // symmetry term
        retRow[jj] += -1.0 * ansatzValues[jj] * functionValue * (testGradients[ii][0] * unitOuterNormal);
        // penalty term
        retRow[jj] += penalty * ansatzValues[jj] * testValues[ii];
      } // loop over all ansatz basis functions
    } // loop over all test basis functions
  } // void evaluate(...) const

  const LocalizableFunctionType& inducingFunction_;
  const double beta_;
}; // class BoundaryDirichletLHS


// forward
template <class LocalizableDiffusionFunctionImp, class LocalizableDirichletFunctionImp>
class BoundaryRHS;


template <class LocalizableDiffusionFunctionImp, class LocalizableDirichletFunctionImp>
class BoundaryRHSTraits
{
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, LocalizableDiffusionFunctionImp>::value,
                "LocalizableDiffusionFunctionImp has to be tagged as Stuff::IsLocalizableFunction!");
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, LocalizableDirichletFunctionImp>::value,
                "LocalizableDirichletFunctionImp has to be tagged as Stuff::IsLocalizableFunction!");
  static_assert(std::is_same<typename LocalizableDiffusionFunctionImp::EntityType,
                             typename LocalizableDirichletFunctionImp::EntityType>::value,
                "EntityTypes have to agree!");
  static_assert(std::is_same<typename LocalizableDiffusionFunctionImp::DomainFieldType,
                             typename LocalizableDirichletFunctionImp::DomainFieldType>::value,
                "DomainFieldTypes have to agree!");
  static_assert(LocalizableDiffusionFunctionImp::dimDomain == LocalizableDirichletFunctionImp::dimDomain,
                "Dimensions of domains have to agree");

public:
  typedef BoundaryRHS<LocalizableDiffusionFunctionImp, LocalizableDirichletFunctionImp> derived_type;
  typedef LocalizableDiffusionFunctionImp LocalizableDiffusionFunctionType;
  typedef LocalizableDirichletFunctionImp LocalizableDirichletFunctionType;
  typedef typename LocalizableDiffusionFunctionType::LocalfunctionType LocalDiffusionFunctionType;
  typedef typename LocalizableDirichletFunctionType::LocalfunctionType LocalDirichletFunctionType;
  typedef std::tuple<std::shared_ptr<LocalDiffusionFunctionType>, std::shared_ptr<LocalDirichletFunctionType>>
      LocalfunctionTupleType;
  typedef typename LocalizableDiffusionFunctionType::EntityType EntityType;
  typedef typename LocalizableDiffusionFunctionType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = LocalizableDiffusionFunctionType::dimDomain;
};


template <class LocalizableDiffusionFunctionImp, class LocalizableDirichletFunctionImp>
class BoundaryRHS
    : public LocalEvaluation::
          Codim1Interface<BoundaryRHSTraits<LocalizableDiffusionFunctionImp, LocalizableDirichletFunctionImp>, 1>
{
public:
  typedef BoundaryRHSTraits<LocalizableDiffusionFunctionImp, LocalizableDirichletFunctionImp> Traits;
  typedef typename Traits::LocalizableDiffusionFunctionType LocalizableDiffusionFunctionType;
  typedef typename Traits::LocalizableDirichletFunctionType LocalizableDirichletFunctionType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = Traits::dimDomain;

  BoundaryRHS(const LocalizableDiffusionFunctionType& diffusion, const LocalizableDirichletFunctionType& dirichlet,
              const double beta = 1.0 / (dimDomain - 1.0))
    : diffusion_(diffusion)
    , dirichlet_(dirichlet)
    , beta_(beta)
  {
  }

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_.local_function(entity), dirichlet_.local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class R, int r, int rC>
  static size_t
  order(const LocalfunctionTupleType localFuncs,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, rC>& testBase)
  {
    const auto localDiffusion = std::get<0>(localFuncs);
    const auto localDirichlet = std::get<1>(localFuncs);
    return redirect_order(*localDiffusion, *localDirichlet, testBase);
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class IntersectionType, class R, int r, int rC>
  void evaluate(const LocalfunctionTupleType localFuncs,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, rC>& testBase,
                const IntersectionType& intersection,
                const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint, Dune::DynamicVector<R>& ret) const
  {
    const auto localDiffusion = std::get<0>(localFuncs);
    const auto localDirichlet = std::get<1>(localFuncs);
    redirect_evaluate(*localDiffusion, *localDirichlet, testBase, intersection, localPoint, ret);
  }

private:
  /**
   *  \return std::max(testOrder + dirichletOrder, diffusionOrder + testGradientOrder + dirichletOrder);
   */
  template <class R, int rLF, int rCLF, int rLR, int rCLR, int rT, int rCT>
  static size_t redirect_order(
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rLF, rCLF>& localDiffusion,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rLR, rCLR>& localDirichlet,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase)
  {
    const size_t testOrder         = testBase.order();
    const size_t testGradientOrder = std::max(ssize_t(testOrder) - 1, ssize_t(0));
    const size_t diffusionOrder    = localDiffusion.order();
    const size_t dirichletOrder = localDirichlet.order();
    return std::max(testOrder + dirichletOrder, diffusionOrder + testGradientOrder + dirichletOrder);
  } // static size_t redirect_order(...)

  template <class IntersectionType, class R, int rLDF, int rCLDF, int rLDR, int rCLDR, int rT, int rCT>
  void redirect_evaluate(
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rLDF, rCLDF>& /*localDiffusion*/,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rLDR, rCLDR>& /*localDirichlet*/,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& /*testBase*/,
      const IntersectionType& /*intersection*/, const Dune::FieldVector<DomainFieldType, dimDomain - 1>& /*localPoint*/,
      Dune::DynamicVector<R>& /*ret*/) const
  {
    static_assert(Dune::AlwaysFalse<R>::value, "Not implemented for these dimensions!");
  } // void redirect_evaluate(...) const

  template <class IntersectionType, class R>
  void redirect_evaluate(
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localDiffusion,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localDirichlet,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& testBase,
      const IntersectionType& intersection, const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint,
      Dune::DynamicVector<R>& ret) const
  {
    // clear ret
    ret *= 0.0;
    typedef typename Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>::DomainType
        DomainType;
    typedef
        typename Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>::RangeType RangeType;
    typedef
        typename Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>::JacobianRangeType
            JacobianRangeType;
    // get local point (which is in intersection coordinates) in entity coordinates
    const DomainType localPointEntity = intersection.geometryInInside().global(localPoint);
    const DomainType unitOuterNormal  = intersection.unitOuterNormal(localPoint);
    // evaluate local functions
    const RangeType diffusionValue = localDiffusion.evaluate(localPointEntity);
    const RangeType dirichletValue = localDirichlet.evaluate(localPointEntity);
    // compute penalty (see Epshteyn, Riviere, 2007)
    const size_t polorder = testBase.order();
    R sigma = 1.0;
    if (polorder <= 1)
      sigma *= 14.0;
    else if (polorder <= 2)
      sigma *= 38.0;
    else if (polorder <= 3)
      sigma *= 74.0;
    else {
#ifndef NDEBUG
#ifndef DUNE_GDT_LOCALEVALUATION_SIPDG_DISABLE_WARNINGS
      std::cout << "\n" << Dune::Stuff::Common::colorString("WARNING(dune.gdt.localevaluation.sipdg.boundaryrhs):")
                << " a polynomial order of " << polorder << " is untested!" << std::endl
                << "  (#define DUNE_GDT_LOCALEVALUATION_SIPDG_DISABLE_WARNINGS to disable this warning)!" << std::endl;
#endif
#endif
      sigma *= 100.0;
    }
    const R penalty = sigma / std::pow(intersection.geometry().volume(), beta_);
    // evaluate basis
    const size_t size = testBase.size();
    std::vector<RangeType> testValues(size, RangeType(0));
    std::vector<JacobianRangeType> testGradients(size, JacobianRangeType(0));
    testBase.evaluate(localPointEntity, testValues);
    testBase.jacobian(localPointEntity, testGradients);
    // compute
    assert(ret.size() >= size);
    // loop over all test basis functions
    for (size_t ii = 0; ii < size; ++ii) {
      // symmetry term
      ret[ii] += -1.0 * dirichletValue * diffusionValue * (testGradients[ii][0] * unitOuterNormal);
      // penalty term
      ret[ii] += penalty * dirichletValue * testValues[ii];
    } // loop over all test basis functions
  } // void redirect_evaluate(...) const

  const LocalizableDiffusionFunctionType& diffusion_;
  const LocalizableDirichletFunctionType& dirichlet_;
  const double beta_;
}; // class BoundaryDirichletRHS


} // namespace SIPDG
} // namespace LocalEvaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALEVALUATION_SIPDG_HH
