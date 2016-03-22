// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LOCALEVALUATION_SWIPDG_HH
#define DUNE_GDT_LOCALEVALUATION_SWIPDG_HH

#include <tuple>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/densematrix.hh>

#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/functions/interfaces.hh>

#include "interface.hh"
#include "sipdg.hh"

namespace Dune {
namespace GDT {
namespace LocalEvaluation {

/**
 *  \brief      Contains local evaluations for the symmetric weighted interior penalty discontinuous Galerkin (SWIPDG)
 *              discretization.
 *
 *              For the choice of penalization and the role of the user input see Epshteyn, Riviere (2007):
 *              "Estimation of penalty parameters for symmetric interior penalty Galerkin methods"
 *              For the coice of the weighting see Ern, Stephansen, Zunino (2007): "A discontinuous Galerkin method with
 *              weighted averages for advection-diffusion equations with locally small and anisotropic diffusivity"
 */
namespace SWIPDG {


// forwards
template <class DiffusionFactorImp, class DiffusionTensorImp = void>
class Inner;


template <class DiffusionFactorImp, class DiffusionTensorImp = void>
class BoundaryLHS;


template <class DiffusionFactorImp, class DiffusionTensorImp>
class InnerPenalty;


template <class DiffusionFactorImp, class DiffusionTensorImp>
class BoundaryLHSPenalty;


template <class DiffusionFactorImp, class DirichletImp, class DiffusionTensorImp = void>
class BoundaryRHS;


namespace internal {


template <class DiffusionFactorImp, class DiffusionTensorImp>
class InnerTraits
{
  static_assert(Stuff::is_localizable_function<DiffusionFactorImp>::value,
                "DiffusionFactorImp has to be a localizable function!");
  static_assert(Stuff::is_localizable_function<DiffusionTensorImp>::value,
                "DiffusionTensorImp has to be a localizable function!");
  static_assert(std::is_same<typename DiffusionFactorImp::EntityType, typename DiffusionTensorImp::EntityType>::value,
                "EntityTypes have to agree!");
  static_assert(
      std::is_same<typename DiffusionFactorImp::DomainFieldType, typename DiffusionTensorImp::DomainFieldType>::value,
      "DomainFieldTypes have to agree!");
  static_assert(DiffusionFactorImp::dimDomain == DiffusionTensorImp::dimDomain, "Dimensions have to agree!");

public:
  typedef Inner<DiffusionFactorImp, DiffusionTensorImp> derived_type;
  typedef DiffusionFactorImp LocalizableDiffusionFactorFunctionType;
  typedef DiffusionTensorImp LocalizableDiffusionTensorFunctionType;
  typedef typename LocalizableDiffusionFactorFunctionType::LocalfunctionType LocalDiffusionFactorFunctionType;
  typedef typename LocalizableDiffusionTensorFunctionType::LocalfunctionType LocalDiffusionTensorFunctionType;
  typedef std::tuple<std::shared_ptr<LocalDiffusionFactorFunctionType>,
                     std::shared_ptr<LocalDiffusionTensorFunctionType>> LocalfunctionTupleType;
  typedef typename LocalizableDiffusionFactorFunctionType::EntityType EntityType;
  typedef typename LocalizableDiffusionFactorFunctionType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = LocalizableDiffusionFactorFunctionType::dimDomain;
}; // class InnerTraits


template <class LocalizableFunctionImp>
class InnerTraits<LocalizableFunctionImp, void>
{
  static_assert(Stuff::is_localizable_function<LocalizableFunctionImp>::value,
                "LocalizableFunctionImp has to be a localizable function!");

public:
  typedef Inner<LocalizableFunctionImp, void> derived_type;
  typedef LocalizableFunctionImp LocalizableFunctionType;
  typedef typename LocalizableFunctionType::EntityType EntityType;
  typedef typename LocalizableFunctionType::DomainFieldType DomainFieldType;
  typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;
  typedef std::tuple<std::shared_ptr<LocalfunctionType>> LocalfunctionTupleType;
  static const size_t dimDomain = LocalizableFunctionType::dimDomain;
}; // class InnerTraits< ..., void >


template <class DiffusionFactorImp, class DiffusionTensorImp>
class InnerPenaltyTraits
{
  static_assert(Stuff::is_localizable_function<DiffusionFactorImp>::value,
                "DiffusionFactorImp has to be a localizable function!");
  static_assert(Stuff::is_localizable_function<DiffusionTensorImp>::value,
                "DiffusionTensorImp has to be a localizable function!");
  static_assert(std::is_same<typename DiffusionFactorImp::EntityType, typename DiffusionTensorImp::EntityType>::value,
                "EntityTypes have to agree!");
  static_assert(
      std::is_same<typename DiffusionFactorImp::DomainFieldType, typename DiffusionTensorImp::DomainFieldType>::value,
      "DomainFieldTypes have to agree!");
  static_assert(DiffusionFactorImp::dimDomain == DiffusionTensorImp::dimDomain, "Dimensions have to agree!");

public:
  typedef InnerPenalty<DiffusionFactorImp, DiffusionTensorImp> derived_type;
  typedef DiffusionFactorImp LocalizableDiffusionFactorFunctionType;
  typedef DiffusionTensorImp LocalizableDiffusionTensorFunctionType;
  typedef typename LocalizableDiffusionFactorFunctionType::LocalfunctionType LocalDiffusionFactorFunctionType;
  typedef typename LocalizableDiffusionTensorFunctionType::LocalfunctionType LocalDiffusionTensorFunctionType;
  typedef std::tuple<std::shared_ptr<LocalDiffusionFactorFunctionType>,
                     std::shared_ptr<LocalDiffusionTensorFunctionType>> LocalfunctionTupleType;
  typedef typename LocalizableDiffusionFactorFunctionType::EntityType EntityType;
  typedef typename LocalizableDiffusionFactorFunctionType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = LocalizableDiffusionFactorFunctionType::dimDomain;
}; // class InnerPenaltyTraits


template <class DiffusionFactorImp, class DiffusionTensorImp>
class BoundaryLHSTraits
{
  static_assert(Stuff::is_localizable_function<DiffusionFactorImp>::value,
                "DiffusionFactorImp has to be a localizable function!");
  static_assert(Stuff::is_localizable_function<DiffusionTensorImp>::value,
                "DiffusionTensorImp has to be a localizable function!");
  static_assert(std::is_same<typename DiffusionFactorImp::EntityType, typename DiffusionTensorImp::EntityType>::value,
                "EntityTypes have to agree!");
  static_assert(
      std::is_same<typename DiffusionFactorImp::DomainFieldType, typename DiffusionTensorImp::DomainFieldType>::value,
      "DomainFieldTypes have to agree!");
  static_assert(DiffusionFactorImp::dimDomain == DiffusionTensorImp::dimDomain, "Dimensions have to agree!");

public:
  typedef BoundaryLHS<DiffusionFactorImp, DiffusionTensorImp> derived_type;
  typedef DiffusionFactorImp LocalizableDiffusionFactorFunctionType;
  typedef DiffusionTensorImp LocalizableDiffusionTensorFunctionType;
  typedef typename LocalizableDiffusionFactorFunctionType::LocalfunctionType LocalDiffusionFactorFunctionType;
  typedef typename LocalizableDiffusionTensorFunctionType::LocalfunctionType LocalDiffusionTensorFunctionType;
  typedef std::tuple<std::shared_ptr<LocalDiffusionFactorFunctionType>,
                     std::shared_ptr<LocalDiffusionTensorFunctionType>> LocalfunctionTupleType;
  typedef typename LocalizableDiffusionFactorFunctionType::EntityType EntityType;
  typedef typename LocalizableDiffusionFactorFunctionType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = LocalizableDiffusionFactorFunctionType::dimDomain;
}; // class BoundaryLHSTraits


template <class DiffusionFactorImp, class DiffusionTensorImp>
class BoundaryLHSPenaltyTraits
{
  static_assert(Stuff::is_localizable_function<DiffusionFactorImp>::value,
                "DiffusionFactorImp has to be a localizable function!");
  static_assert(Stuff::is_localizable_function<DiffusionTensorImp>::value,
                "DiffusionTensorImp has to be a localizable function!");
  static_assert(std::is_same<typename DiffusionFactorImp::EntityType, typename DiffusionTensorImp::EntityType>::value,
                "EntityTypes have to agree!");
  static_assert(
      std::is_same<typename DiffusionFactorImp::DomainFieldType, typename DiffusionTensorImp::DomainFieldType>::value,
      "DomainFieldTypes have to agree!");
  static_assert(DiffusionFactorImp::dimDomain == DiffusionTensorImp::dimDomain, "Dimensions have to agree!");

public:
  typedef BoundaryLHSPenalty<DiffusionFactorImp, DiffusionTensorImp> derived_type;
  typedef DiffusionFactorImp LocalizableDiffusionFactorFunctionType;
  typedef DiffusionTensorImp LocalizableDiffusionTensorFunctionType;
  typedef typename LocalizableDiffusionFactorFunctionType::LocalfunctionType LocalDiffusionFactorFunctionType;
  typedef typename LocalizableDiffusionTensorFunctionType::LocalfunctionType LocalDiffusionTensorFunctionType;
  typedef std::tuple<std::shared_ptr<LocalDiffusionFactorFunctionType>,
                     std::shared_ptr<LocalDiffusionTensorFunctionType>> LocalfunctionTupleType;
  typedef typename LocalizableDiffusionFactorFunctionType::EntityType EntityType;
  typedef typename LocalizableDiffusionFactorFunctionType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = LocalizableDiffusionFactorFunctionType::dimDomain;
}; // class BoundaryLHSPenaltyTraits


template <class LocalizableFunctionImp>
class BoundaryLHSTraits<LocalizableFunctionImp, void>
{
  static_assert(Stuff::is_localizable_function<LocalizableFunctionImp>::value,
                "LocalizableFunctionImp has to be a localizable function!");

public:
  typedef LocalizableFunctionImp LocalizableFunctionType;
  typedef BoundaryLHS<LocalizableFunctionType> derived_type;
  typedef typename LocalizableFunctionType::EntityType EntityType;
  typedef typename LocalizableFunctionType::DomainFieldType DomainFieldType;
  typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;
  typedef std::tuple<std::shared_ptr<LocalfunctionType>> LocalfunctionTupleType;
  static const size_t dimDomain = LocalizableFunctionType::dimDomain;
}; // class BoundaryLHSTraits< ..., void >


template <class DiffusionFactorImp, class DirichletImp, class DiffusionTensorImp>
class BoundaryRHSTraits
{
  static_assert(Stuff::is_localizable_function<DiffusionFactorImp>::value,
                "DiffusionFactorImp has to be a localizable function!");
  static_assert(Stuff::is_localizable_function<DirichletImp>::value, "DirichletImp has to be a localizable function!");
  static_assert(Stuff::is_localizable_function<DiffusionTensorImp>::value,
                "DiffusionTensorImp has to be a localizable function!");
  static_assert(std::is_same<typename DiffusionFactorImp::EntityType, typename DiffusionTensorImp::EntityType>::value
                    && std::is_same<typename DiffusionFactorImp::EntityType, typename DirichletImp::EntityType>::value,
                "EntityTypes have to agree!");
  static_assert(
      std::is_same<typename DiffusionFactorImp::DomainFieldType, typename DiffusionTensorImp::DomainFieldType>::value
          && std::is_same<typename DiffusionFactorImp::DomainFieldType, typename DirichletImp::DomainFieldType>::value,
      "DomainFieldTypes have to agree!");
  static_assert(DiffusionFactorImp::dimDomain == DiffusionTensorImp::dimDomain
                    && DiffusionFactorImp::dimDomain == DirichletImp::dimDomain,
                "Dimensions have to agree!");

public:
  typedef BoundaryRHS<DiffusionFactorImp, DirichletImp, DiffusionTensorImp> derived_type;
  typedef DiffusionFactorImp LocalizableDiffusionFactorFunctionType;
  typedef DirichletImp LocalizableDirichletFunctionType;
  typedef DiffusionTensorImp LocalizableDiffusionTensorFunctionType;
  typedef typename LocalizableDiffusionFactorFunctionType::LocalfunctionType LocalDiffusionFactorFunctionType;
  typedef typename LocalizableDirichletFunctionType::LocalfunctionType LocalDirichletFunctionType;
  typedef typename LocalizableDiffusionTensorFunctionType::LocalfunctionType LocalDiffusionTensorFunctionType;
  typedef std::tuple<std::shared_ptr<LocalDiffusionFactorFunctionType>,
                     std::shared_ptr<LocalDiffusionTensorFunctionType>,
                     std::shared_ptr<LocalDirichletFunctionType>> LocalfunctionTupleType;
  typedef typename LocalizableDiffusionFactorFunctionType::EntityType EntityType;
  typedef typename LocalizableDiffusionFactorFunctionType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = LocalizableDiffusionFactorFunctionType::dimDomain;
}; // class BoundaryRHSTraits


template <class LocalizableDiffusionFunctionImp, class LocalizableDirichletFunctionImp>
class BoundaryRHSTraits<LocalizableDiffusionFunctionImp, LocalizableDirichletFunctionImp, void>
{
  static_assert(Stuff::is_localizable_function<LocalizableDiffusionFunctionImp>::value,
                "LocalizableDiffusionFunctionImp has to be a localizable function!");
  static_assert(Stuff::is_localizable_function<LocalizableDirichletFunctionImp>::value,
                "LocalizableDirichletFunctionImp has to be a localizable function!");
  static_assert(std::is_same<typename LocalizableDiffusionFunctionImp::EntityType,
                             typename LocalizableDirichletFunctionImp::EntityType>::value,
                "EntityTypes have to agree!");
  static_assert(std::is_same<typename LocalizableDiffusionFunctionImp::DomainFieldType,
                             typename LocalizableDirichletFunctionImp::DomainFieldType>::value,
                "DomainFieldTypes have to agree!");
  static_assert(LocalizableDiffusionFunctionImp::dimDomain == LocalizableDirichletFunctionImp::dimDomain,
                "Dimensions have to agree!");

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
  static const size_t dimDomain = LocalizableDiffusionFunctionType::dimDomain;
}; // class BoundaryRHSTraits< ..., void >


} // namespace internal


/**
 *  see Epshteyn, Riviere, 2007 for the meaning of beta
 */
template <class LocalizableFunctionImp>
class Inner<LocalizableFunctionImp, void>
    : public LocalEvaluation::Codim1Interface<internal::InnerTraits<LocalizableFunctionImp, void>, 4>
{
public:
  typedef internal::InnerTraits<LocalizableFunctionImp, void> Traits;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const size_t dimDomain = Traits::dimDomain;

  Inner(const LocalizableFunctionType& inducingFunction,
        const double beta = SIPDG::internal::default_beta(LocalizableFunctionType::dimDomain))
    : inducingFunction_(inducingFunction)
    , beta_(beta)
  {
  }

  /// \name Required by LocalEvaluation::Codim1Interface< ..., 4 >
  /// \{

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(inducingFunction_.local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t
  order(const LocalfunctionTupleType& localFunctionsEntity, const LocalfunctionTupleType& localFunctionsNeighbor,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBaseEntity,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBaseEntity,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBaseNeighbor,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBaseNeighbor)
      const
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
  template <class IntersectionType, class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
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

  /// \}
  /// \name Actual implementation of order
  /// \{

  template <class R, size_t rL, size_t rCL, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t
  order(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& localFunctionEntity,
        const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& localFunctionNeighbor,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBaseEntity,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBaseEntity,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBaseNeighbor,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBaseNeighbor)
      const
  {
    return std::max(localFunctionEntity.order(), localFunctionNeighbor.order())
           + std::max(testBaseEntity.order(), testBaseNeighbor.order())
           + std::max(ansatzBaseEntity.order(), ansatzBaseNeighbor.order());
  }

  /// \}
  /// \name Actual implementation of evaluate
  /// \{

  /**
   *  \brief  Computes the swipdg fluxes in a primal setting.
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
    // convert local point (which is in intersection coordinates) to entity/neighbor coordinates
    const auto localPointEn    = intersection.geometryInInside().global(localPoint);
    const auto localPointNe    = intersection.geometryInOutside().global(localPoint);
    const auto unitOuterNormal = intersection.unitOuterNormal(localPoint);
    // evaluate local function
    const auto functionValueEn = localFunctionEntity.evaluate(localPointEn);
    const auto functionValueNe = localFunctionNeighbor.evaluate(localPointNe);
    // compute penalty factor (see Epshteyn, Riviere, 2007)
    const size_t max_polorder =
        std::max(testBaseEntity.order(),
                 std::max(ansatzBaseEntity.order(), std::max(testBaseNeighbor.order(), ansatzBaseNeighbor.order())));
    const R sigma = SIPDG::internal::inner_sigma(max_polorder);
    // compute weighting (see Ern, Stephansen, Zunino 2007)
    const R delta_plus   = /*unitOuterNormal * (*/ functionValueNe /** unitOuterNormal)*/;
    const R delta_minus  = /*unitOuterNormal * (*/ functionValueEn /** unitOuterNormal)*/;
    const R gamma        = (delta_plus * delta_minus) / (delta_plus + delta_minus);
    const R penalty      = (sigma * gamma) / std::pow(intersection.geometry().volume(), beta_);
    const R weight_plus  = delta_minus / (delta_plus + delta_minus);
    const R weight_minus = delta_plus / (delta_plus + delta_minus);
    // evaluate bases
    // * entity
    //   * test
    const size_t rowsEn        = testBaseEntity.size();
    const auto testValuesEn    = testBaseEntity.evaluate(localPointEn);
    const auto testGradientsEn = testBaseEntity.jacobian(localPointEn);
    //   * ansatz
    const size_t colsEn          = ansatzBaseEntity.size();
    const auto ansatzValuesEn    = ansatzBaseEntity.evaluate(localPointEn);
    const auto ansatzGradientsEn = ansatzBaseEntity.jacobian(localPointEn);
    // * neighbor
    //   * test
    const size_t rowsNe        = testBaseNeighbor.size();
    const auto testValuesNe    = testBaseNeighbor.evaluate(localPointNe);
    const auto testGradientsNe = testBaseNeighbor.jacobian(localPointNe);
    //   * ansatz
    const size_t colsNe          = ansatzBaseNeighbor.size();
    const auto ansatzValuesNe    = ansatzBaseNeighbor.evaluate(localPointNe);
    const auto ansatzGradientsNe = ansatzBaseNeighbor.jacobian(localPointNe);
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
            -weight_minus * functionValueEn * (ansatzGradientsEn[jj][0] * unitOuterNormal) * testValuesEn[ii];
        // symmetry term
        entityEntityRetRow[jj] +=
            -weight_minus * ansatzValuesEn[jj] * functionValueEn * (testGradientsEn[ii][0] * unitOuterNormal);
        // penalty term
        entityEntityRetRow[jj] += penalty * ansatzValuesEn[jj] * testValuesEn[ii];
      } // loop over all entity ansatz basis functions
      // loop over all neighbor ansatz basis functions
      for (size_t jj = 0; jj < colsNe; ++jj) {
        // consistency term
        entityNeighborRetRow[jj] +=
            -weight_plus * functionValueNe * (ansatzGradientsNe[jj][0] * unitOuterNormal) * testValuesEn[ii];
        // symmetry term
        entityNeighborRetRow[jj] +=
            weight_minus * ansatzValuesNe[jj] * functionValueEn * (testGradientsEn[ii][0] * unitOuterNormal);
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
            weight_minus * functionValueEn * (ansatzGradientsEn[jj][0] * unitOuterNormal) * testValuesNe[ii];
        // symmetry term
        neighborEntityRetRow[jj] +=
            -weight_plus * ansatzValuesEn[jj] * functionValueNe * (testGradientsNe[ii][0] * unitOuterNormal);
        // penalty term
        neighborEntityRetRow[jj] += -1.0 * penalty * ansatzValuesEn[jj] * testValuesNe[ii];
      } // loop over all entity ansatz basis functions
      // loop over all neighbor ansatz basis functions
      for (size_t jj = 0; jj < colsNe; ++jj) {
        // consistency term
        neighborNeighborRetRow[jj] +=
            weight_plus * functionValueNe * (ansatzGradientsNe[jj][0] * unitOuterNormal) * testValuesNe[ii];
        // symmetry term
        neighborNeighborRetRow[jj] +=
            weight_plus * ansatzValuesNe[jj] * functionValueNe * (testGradientsNe[ii][0] * unitOuterNormal);
        // penalty term
        neighborNeighborRetRow[jj] += penalty * ansatzValuesNe[jj] * testValuesNe[ii];
      } // loop over all neighbor ansatz basis functions
    } // loop over all neighbor test basis functions
  } // ... evaluate(...)

  /// \}

private:
  const LocalizableFunctionType& inducingFunction_;
  const double beta_;
}; // class Inner< ..., void >


template <class LocalizableFunctionImp>
class BoundaryLHS<LocalizableFunctionImp, void>
    : public LocalEvaluation::Codim1Interface<internal::BoundaryLHSTraits<LocalizableFunctionImp, void>, 2>
{
public:
  typedef internal::BoundaryLHSTraits<LocalizableFunctionImp, void> Traits;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const size_t dimDomain = Traits::dimDomain;

  BoundaryLHS(const LocalizableFunctionType& inducingFunction,
              const double beta = SIPDG::internal::default_beta(dimDomain))
    : inducingFunction_(inducingFunction)
    , beta_(beta)
  {
  }

  /// \name Required by LocalEvaluation::Codim1Interface< ..., 2 >
  /// \{

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(inducingFunction_.local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t
  order(const LocalfunctionTupleType localFuncs,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase) const
  {
    return order(*std::get<0>(localFuncs), testBase, ansatzBase);
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class IntersectionType, class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  void evaluate(const LocalfunctionTupleType localFuncs,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase,
                const IntersectionType& intersection,
                const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    evaluate(*std::get<0>(localFuncs), testBase, ansatzBase, intersection, localPoint, ret);
  }

  /// \}
  /// \name Actual implementation of order
  /// \{

  /**
   *  \return localFunction.order() + testBase.order() + ansatzBase.order()
   */
  template <class R, size_t rL, size_t rCL, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t
  order(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& localFunction,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase) const
  {
    return localFunction.order() + testBase.order() + ansatzBase.order();
  }

  /// \}
  /// \name Actual implementation of evaluate
  /// \{

  template <class IntersectionType, class R>
  void evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, 2, R, 1, 1>& localFunction,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, 2, R, 1, 1>& testBase,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, 2, R, 1, 1>& ansatzBase,
                const IntersectionType& intersection, const Dune::FieldVector<DomainFieldType, 1>& localPoint,
                Dune::DynamicMatrix<R>& ret) const
  {
    // clear ret
    ret *= 0.0;
    // get local point (which is in intersection coordinates) in entity coordinates
    const auto localPointEntity = intersection.geometryInInside().global(localPoint);
    const auto unitOuterNormal  = intersection.unitOuterNormal(localPoint);
    // evaluate local function
    const auto functionValue = localFunction.evaluate(localPointEntity);
    // compute penalty (see Epshteyn, Riviere, 2007)
    const size_t max_polorder = std::max(testBase.order(), ansatzBase.order());
    const R sigma             = SIPDG::internal::boundary_sigma(max_polorder);
    // compute weighting (see Ern, Stephansen, Zunino 2007)
    const R gamma   = /*unitOuterNormal * (*/ functionValue /** unitOuterNormal)*/;
    const R penalty = (sigma * gamma) / std::pow(intersection.geometry().volume(), beta_);
    // evaluate bases
    // * test
    const size_t rows        = testBase.size();
    const auto testValues    = testBase.evaluate(localPointEntity);
    const auto testGradients = testBase.jacobian(localPointEntity);
    // * ansatz
    const size_t cols          = ansatzBase.size();
    const auto ansatzValues    = ansatzBase.evaluate(localPointEntity);
    const auto ansatzGradients = ansatzBase.jacobian(localPointEntity);
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
  } // ... evaluate(...)

  /// \}

private:
  const LocalizableFunctionType& inducingFunction_;
  const double beta_;
}; // class BoundaryLHS< ..., void >


template <class LocalizableDiffusionFunctionImp, class LocalizableDirichletFunctionImp>
class BoundaryRHS<LocalizableDiffusionFunctionImp, LocalizableDirichletFunctionImp, void>
    : public LocalEvaluation::Codim1Interface<internal::BoundaryRHSTraits<LocalizableDiffusionFunctionImp,
                                                                          LocalizableDirichletFunctionImp, void>,
                                              1>
{
public:
  typedef internal::BoundaryRHSTraits<LocalizableDiffusionFunctionImp, LocalizableDirichletFunctionImp, void> Traits;
  typedef typename Traits::LocalizableDiffusionFunctionType LocalizableDiffusionFunctionType;
  typedef typename Traits::LocalizableDirichletFunctionType LocalizableDirichletFunctionType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const size_t dimDomain = Traits::dimDomain;

  BoundaryRHS(const LocalizableDiffusionFunctionType& diffusion, const LocalizableDirichletFunctionType& dirichlet,
              const double beta = SIPDG::internal::default_beta(dimDomain))
    : diffusion_(diffusion)
    , dirichlet_(dirichlet)
    , beta_(beta)
  {
  }

  /// \name Required by LocalEvaluation::Codim1Interface< ..., 1 >
  /// \{

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_.local_function(entity), dirichlet_.local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class R, size_t r, size_t rC>
  size_t order(const LocalfunctionTupleType localFuncs,
               const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, rC>& testBase) const
  {
    const auto localDiffusion = std::get<0>(localFuncs);
    const auto localDirichlet = std::get<1>(localFuncs);
    return order(*localDiffusion, *localDirichlet, testBase);
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class IntersectionType, class R, size_t r, size_t rC>
  void evaluate(const LocalfunctionTupleType localFuncs,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, rC>& testBase,
                const IntersectionType& intersection,
                const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint, Dune::DynamicVector<R>& ret) const
  {
    const auto localDiffusion = std::get<0>(localFuncs);
    const auto localDirichlet = std::get<1>(localFuncs);
    evaluate(*localDiffusion, *localDirichlet, testBase, intersection, localPoint, ret);
  }

  /// \}
  /// \name Actual implementation of order
  /// \{

  /**
   *  \return std::max(testOrder + dirichletOrder, diffusionOrder + testGradientOrder + dirichletOrder);
   */
  template <class R, size_t rLF, size_t rCLF, size_t rLR, size_t rCLR, size_t rT, size_t rCT>
  size_t
  order(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rLF, rCLF>& localDiffusion,
        const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rLR, rCLR>& localDirichlet,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase) const
  {
    const size_t testOrder         = testBase.order();
    const size_t testGradientOrder = boost::numeric_cast<size_t>(std::max(ssize_t(testOrder) - 1, ssize_t(0)));
    const size_t diffusionOrder    = localDiffusion.order();
    const size_t dirichletOrder = localDirichlet.order();
    return std::max(testOrder + dirichletOrder, diffusionOrder + testGradientOrder + dirichletOrder);
  } // ... order(...)

  /// \}
  /// \name Actual implementation of evaluate
  /// \{

  template <class IntersectionType, class R>
  void evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localDiffusion,
                const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localDirichlet,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& testBase,
                const IntersectionType& intersection,
                const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint, Dune::DynamicVector<R>& ret) const
  {
    // clear ret
    ret *= 0.0;
    // get local point (which is in intersection coordinates) in entity coordinates
    const auto localPointEntity = intersection.geometryInInside().global(localPoint);
    const auto unitOuterNormal  = intersection.unitOuterNormal(localPoint);
    // evaluate local functions
    const auto diffusionValue = localDiffusion.evaluate(localPointEntity);
    const auto dirichletValue = localDirichlet.evaluate(localPointEntity);
    // compute penalty (see Epshteyn, Riviere, 2007)
    const auto polorder = testBase.order();
    const R sigma       = SIPDG::internal::boundary_sigma(polorder);
    // compute weighting (see Ern, Stephansen, Zunino 2007)
    const R gamma   = /*unitOuterNormal * (*/ diffusionValue /** unitOuterNormal)*/;
    const R penalty = (sigma * gamma) / std::pow(intersection.geometry().volume(), beta_);
    // evaluate basis
    const auto size          = testBase.size();
    const auto testValues    = testBase.evaluate(localPointEntity);
    const auto testGradients = testBase.jacobian(localPointEntity);
    // compute
    assert(ret.size() >= size);
    // loop over all test basis functions
    for (size_t ii = 0; ii < size; ++ii) {
      // symmetry term
      ret[ii] += -1.0 * dirichletValue * diffusionValue * (testGradients[ii][0] * unitOuterNormal);
      // penalty term
      ret[ii] += penalty * dirichletValue * testValues[ii];
    } // loop over all test basis functions
  } // ... evaluate(...)

private:
  const LocalizableDiffusionFunctionType& diffusion_;
  const LocalizableDirichletFunctionType& dirichlet_;
  const double beta_;
}; // class BoundaryRHS< ..., void >


} // namespace SWIPDG
} // namespace LocalEvaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALEVALUATION_SWIPDG_HH
