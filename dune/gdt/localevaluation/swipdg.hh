// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LOCALEVALUATION_SWIPDG_HH
#define DUNE_GDT_LOCALEVALUATION_SWIPDG_HH

#include <tuple>
#include <type_traits>

#include <dune/common/densematrix.hh>

#ifndef NDEBUG
#ifndef DUNE_GDT_LOCALEVALUATION_SWIPDG_DISABLE_WARNINGS
#include <dune/stuff/common/logging.hh>
#endif
#endif
#include <dune/stuff/functions/interfaces.hh>

#include "interface.hh"

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


/**
 * \note see Epshteyn, Riviere, 2007
 */
static inline double default_beta(const size_t dimDomain)
{
  return 1.0 / (dimDomain - 1.0);
}


/**
 * \note see Epshteyn, Riviere, 2007
 */
static inline double inner_sigma(const size_t pol_order)
{
  double sigma = 1.0;
  if (pol_order <= 1)
    sigma *= 8.0;
  else if (pol_order <= 2)
    sigma *= 20.0;
  else if (pol_order <= 3)
    sigma *= 38.0;
  else {
#ifndef NDEBUG
#ifndef DUNE_GDT_LOCALEVALUATION_SWIPDG_DISABLE_WARNINGS
    DSC::TimedLogger().get("gdt.localevaluation.swipdg.inner").warn()
        << "a polynomial order of " << pol_order << " is untested!\n"
        << "  #define DUNE_GDT_LOCALEVALUATION_SWIPDG_DISABLE_WARNINGS to statically disable this warning\n"
        << "  or dynamically disable warnings of the TimedLogger() instance!" << std::endl;
#endif
#endif
    sigma *= 50.0;
  }
  return sigma;
} // ... inner_sigma(...)


/**
 * \note see Epshteyn, Riviere, 2007
 */
static inline double boundary_sigma(const size_t pol_order)
{
  double sigma = 1.0;
  if (pol_order <= 1)
    sigma *= 14.0;
  else if (pol_order <= 2)
    sigma *= 38.0;
  else if (pol_order <= 3)
    sigma *= 74.0;
  else {
#ifndef NDEBUG
#ifndef DUNE_GDT_LOCALEVALUATION_SWIPDG_DISABLE_WARNINGS
    DSC::TimedLogger().get("gdt.localevaluation.swipdg.inner").warn()
        << "a polynomial order of " << pol_order << " is untested!\n"
        << "  #define DUNE_GDT_LOCALEVALUATION_SWIPDG_DISABLE_WARNINGS to statically disable this warning\n"
        << "  or dynamically disable warnings of the TimedLogger() instance!" << std::endl;
#endif
#endif
    sigma *= 100.0;
  }
  return sigma;
} // ... boundary_sigma(...)


template <class DiffusionFactorImp, class DiffusionTensorImp>
class InnerTraits
{
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, DiffusionFactorImp>::value,
                "DiffusionFactorImp has to be tagged as Stuff::IsLocalizableFunction!");
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, DiffusionTensorImp>::value,
                "DiffusionTensorImp has to be tagged as Stuff::IsLocalizableFunction!");
  static_assert(std::is_same<typename DiffusionFactorImp::EntityType, typename DiffusionTensorImp::EntityType>::value,
                "EntityTypes have to agree!");
  static_assert(
      std::is_same<typename DiffusionFactorImp::DomainFieldType, typename DiffusionTensorImp::DomainFieldType>::value,
      "DomainFieldTypes have to agree!");
  static_assert(DiffusionFactorImp::dimDomain == DiffusionTensorImp::dimDomain, "Dimensions of domains have to agree");

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
  static const unsigned int dimDomain = LocalizableDiffusionFactorFunctionType::dimDomain;
};


template <class LocalizableFunctionImp>
class InnerTraits<LocalizableFunctionImp, void>
{
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, LocalizableFunctionImp>::value,
                "LocalizableFunctionImp has to be tagged as Stuff::IsLocalizableFunction!");

public:
  typedef Inner<LocalizableFunctionImp, void> derived_type;
  typedef LocalizableFunctionImp LocalizableFunctionType;
  typedef typename LocalizableFunctionType::EntityType EntityType;
  typedef typename LocalizableFunctionType::DomainFieldType DomainFieldType;
  typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;
  typedef std::tuple<std::shared_ptr<LocalfunctionType>> LocalfunctionTupleType;
  static const unsigned int dimDomain = LocalizableFunctionType::dimDomain;
};


template <class DiffusionFactorImp, class DiffusionTensorImp>
class InnerPenaltyTraits
{
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, DiffusionFactorImp>::value,
                "DiffusionFactorImp has to be tagged as Stuff::IsLocalizableFunction!");
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, DiffusionTensorImp>::value,
                "DiffusionTensorImp has to be tagged as Stuff::IsLocalizableFunction!");
  static_assert(std::is_same<typename DiffusionFactorImp::EntityType, typename DiffusionTensorImp::EntityType>::value,
                "EntityTypes have to agree!");
  static_assert(
      std::is_same<typename DiffusionFactorImp::DomainFieldType, typename DiffusionTensorImp::DomainFieldType>::value,
      "DomainFieldTypes have to agree!");
  static_assert(DiffusionFactorImp::dimDomain == DiffusionTensorImp::dimDomain, "Dimensions of domains have to agree");

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
  static const unsigned int dimDomain = LocalizableDiffusionFactorFunctionType::dimDomain;
};


template <class DiffusionFactorImp, class DiffusionTensorImp>
class BoundaryLHSTraits
{
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, DiffusionFactorImp>::value,
                "DiffusionFactorImp has to be tagged as Stuff::IsLocalizableFunction!");
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, DiffusionTensorImp>::value,
                "DiffusionTensorImp has to be tagged as Stuff::IsLocalizableFunction!");
  static_assert(std::is_same<typename DiffusionFactorImp::EntityType, typename DiffusionTensorImp::EntityType>::value,
                "EntityTypes have to agree!");
  static_assert(
      std::is_same<typename DiffusionFactorImp::DomainFieldType, typename DiffusionTensorImp::DomainFieldType>::value,
      "DomainFieldTypes have to agree!");
  static_assert(DiffusionFactorImp::dimDomain == DiffusionTensorImp::dimDomain, "Dimensions of domains have to agree");

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
  static const unsigned int dimDomain = LocalizableDiffusionFactorFunctionType::dimDomain;
};


template <class DiffusionFactorImp, class DiffusionTensorImp>
class BoundaryLHSPenaltyTraits
{
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, DiffusionFactorImp>::value,
                "DiffusionFactorImp has to be tagged as Stuff::IsLocalizableFunction!");
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, DiffusionTensorImp>::value,
                "DiffusionTensorImp has to be tagged as Stuff::IsLocalizableFunction!");
  static_assert(std::is_same<typename DiffusionFactorImp::EntityType, typename DiffusionTensorImp::EntityType>::value,
                "EntityTypes have to agree!");
  static_assert(
      std::is_same<typename DiffusionFactorImp::DomainFieldType, typename DiffusionTensorImp::DomainFieldType>::value,
      "DomainFieldTypes have to agree!");
  static_assert(DiffusionFactorImp::dimDomain == DiffusionTensorImp::dimDomain, "Dimensions of domains have to agree");

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
  static const unsigned int dimDomain = LocalizableDiffusionFactorFunctionType::dimDomain;
};


template <class LocalizableFunctionImp>
class BoundaryLHSTraits<LocalizableFunctionImp, void>
{
  static_assert(std::is_base_of<Dune::Stuff::IsLocalizableFunction, LocalizableFunctionImp>::value,
                "LocalizableFunctionImp has to be derived from Stuff::IsLocalizableFunction.");

public:
  typedef LocalizableFunctionImp LocalizableFunctionType;
  typedef BoundaryLHS<LocalizableFunctionType> derived_type;
  typedef typename LocalizableFunctionType::EntityType EntityType;
  typedef typename LocalizableFunctionType::DomainFieldType DomainFieldType;
  typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;
  typedef std::tuple<std::shared_ptr<LocalfunctionType>> LocalfunctionTupleType;
  static const unsigned int dimDomain = LocalizableFunctionType::dimDomain;
};


template <class DiffusionFactorImp, class DirichletImp, class DiffusionTensorImp>
class BoundaryRHSTraits
{
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, DiffusionFactorImp>::value,
                "DiffusionFactorImp has to be tagged as Stuff::IsLocalizableFunction!");
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, DirichletImp>::value,
                "DirichletImp has to be tagged as Stuff::IsLocalizableFunction!");
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, DiffusionTensorImp>::value,
                "DiffusionTensorImp has to be tagged as Stuff::IsLocalizableFunction!");
  static_assert(std::is_same<typename DiffusionFactorImp::EntityType, typename DiffusionTensorImp::EntityType>::value
                    && std::is_same<typename DiffusionFactorImp::EntityType, typename DirichletImp::EntityType>::value,
                "EntityTypes have to agree!");
  static_assert(
      std::is_same<typename DiffusionFactorImp::DomainFieldType, typename DiffusionTensorImp::DomainFieldType>::value
          && std::is_same<typename DiffusionFactorImp::DomainFieldType, typename DirichletImp::DomainFieldType>::value,
      "DomainFieldTypes have to agree!");
  static_assert(DiffusionFactorImp::dimDomain == DiffusionTensorImp::dimDomain
                    && DiffusionFactorImp::dimDomain == DirichletImp::dimDomain,
                "Dimensions of domains have to agree");

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
  static const unsigned int dimDomain = LocalizableDiffusionFactorFunctionType::dimDomain;
};


template <class LocalizableDiffusionFunctionImp, class LocalizableDirichletFunctionImp>
class BoundaryRHSTraits<LocalizableDiffusionFunctionImp, LocalizableDirichletFunctionImp, void>
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
  static const unsigned int dimDomain = Traits::dimDomain;

  Inner(const LocalizableFunctionType& inducingFunction,
        const double beta = internal::default_beta(LocalizableFunctionType::dimDomain))
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

  template <class R, int rL, int rCL, int rT, int rCT, int rA, int rCA>
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
    // compute penalty factor (see Epshteyn, Riviere, 2007)
    const size_t max_polorder =
        std::max(testBaseEntity.order(),
                 std::max(ansatzBaseEntity.order(), std::max(testBaseNeighbor.order(), ansatzBaseNeighbor.order())));
    const R sigma = internal::inner_sigma(max_polorder);
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
  } // void evaluate< ..., 1, 1 >(...) const

private:
  const LocalizableFunctionType& inducingFunction_;
  const double beta_;
}; // CouplingPrimal


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
  static const unsigned int dimDomain = Traits::dimDomain;

  BoundaryLHS(const LocalizableFunctionType& inducingFunction, const double beta = internal::default_beta(dimDomain))
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
  size_t
  order(const LocalfunctionTupleType localFuncs,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase) const
  {
    const auto localFunction = std::get<0>(localFuncs);
    return order(*localFunction, testBase, ansatzBase);
  }

  /**
   *  \return localFunction.order() + testBase.order() + ansatzBase.order()
   */
  template <class R, int rL, int rCL, int rT, int rCT, int rA, int rCA>
  size_t
  order(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& localFunction,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase) const
  {
    return localFunction.order() + testBase.order() + ansatzBase.order();
  } // size_t order(...)

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
    evaluate(*localFunction, testBase, ansatzBase, intersection, localPoint, ret);
  }

  template <class IntersectionType, class R, int rL, int rCL, int rT, int rCT, int rA, int rCA>
  void
  evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& /*localFunction*/,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& /*testBase*/,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& /*ansatzBase*/,
           const IntersectionType& /*intersection*/,
           const Dune::FieldVector<DomainFieldType, dimDomain - 1>& /*localPoint*/,
           Dune::DynamicMatrix<R>& /*ret*/) const
  {
    static_assert(Dune::AlwaysFalse<R>::value, "Not implemented for these dimensions!");
  } // void evaluate(...) const

  template <class IntersectionType, class R>
  void evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, 2, R, 1, 1>& localFunction,
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
    const R sigma             = internal::boundary_sigma(max_polorder);
    // compute weighting (see Ern, Stephansen, Zunino 2007)
    const R gamma   = /*unitOuterNormal * (*/ functionValue /** unitOuterNormal)*/;
    const R penalty = (sigma * gamma) / std::pow(intersection.geometry().volume(), beta_);
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

private:
  const LocalizableFunctionType& inducingFunction_;
  const double beta_;
}; // class BoundaryLHS


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
  static const unsigned int dimDomain = Traits::dimDomain;

  BoundaryRHS(const LocalizableDiffusionFunctionType& diffusion, const LocalizableDirichletFunctionType& dirichlet,
              const double beta = internal::default_beta(dimDomain))
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
  size_t order(const LocalfunctionTupleType localFuncs,
               const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, rC>& testBase) const
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
  size_t redirect_order(
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rLF, rCLF>& localDiffusion,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rLR, rCLR>& localDirichlet,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase) const
  {
    const size_t testOrder         = testBase.order();
    const size_t testGradientOrder = std::max(ssize_t(testOrder) - 1, ssize_t(0));
    const size_t diffusionOrder    = localDiffusion.order();
    const size_t dirichletOrder = localDirichlet.order();
    return std::max(testOrder + dirichletOrder, diffusionOrder + testGradientOrder + dirichletOrder);
  } // ... redirect_order(...)

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
    const R sigma         = internal::boundary_sigma(polorder);
    // compute weighting (see Ern, Stephansen, Zunino 2007)
    const R gamma   = /*unitOuterNormal * (*/ diffusionValue /** unitOuterNormal)*/;
    const R penalty = (sigma * gamma) / std::pow(intersection.geometry().volume(), beta_);
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
}; // class BoundaryRHS


} // namespace SWIPDG
} // namespace LocalEvaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALEVALUATION_SWIPDG_HH
