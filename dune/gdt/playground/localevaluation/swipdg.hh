// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PLAYGROUND_LOCALEVALUATION_SWIPDG_HH
#define DUNE_GDT_PLAYGROUND_LOCALEVALUATION_SWIPDG_HH

#include <dune/stuff/common/print.hh>
#include <dune/stuff/common/fmatrix.hh>
#include <dune/stuff/common/timedlogging.hh>
#include <dune/stuff/common/type_utils.hh>

#include <dune/gdt/localevaluation/swipdg.hh>

namespace Dune {
namespace GDT {
namespace LocalEvaluation {
namespace SWIPDG {


//! \todo Make use of InnerPenalty
template <class DiffusionFactorType, class DiffusionTensorType>
class Inner
    : public LocalEvaluation::Codim1Interface<internal::InnerTraits<DiffusionFactorType, DiffusionTensorType>, 4>
{
public:
  typedef internal::InnerTraits<DiffusionFactorType, DiffusionTensorType> Traits;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const size_t dimDomain = Traits::dimDomain;

  Inner(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& inducingFunction,
        const double beta = SIPDG::internal::default_beta(dimDomain))
    : diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(inducingFunction)
    , beta_(beta)
  {
  }

  /// \name Required by LocalEvaluation::Codim1Interface< ..., 4 >.
  /// \{

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_factor_.local_function(entity), diffusion_tensor_.local_function(entity));
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
    const auto local_diffusion_factor_entity   = std::get<0>(localFunctionsEntity);
    const auto local_diffusion_tensor_entity   = std::get<1>(localFunctionsEntity);
    const auto local_diffusion_factor_neighbor = std::get<0>(localFunctionsNeighbor);
    const auto local_diffusion_tensor_neighbor = std::get<1>(localFunctionsNeighbor);
    return order(*local_diffusion_factor_entity,
                 *local_diffusion_tensor_entity,
                 *local_diffusion_factor_neighbor,
                 *local_diffusion_tensor_neighbor,
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
    const auto local_diffusion_factor_entity   = std::get<0>(localFunctionsEntity);
    const auto local_diffusion_tensor_entity   = std::get<1>(localFunctionsEntity);
    const auto local_diffusion_factor_neighbor = std::get<0>(localFunctionsNeighbor);
    const auto local_diffusion_tensor_neighbor = std::get<1>(localFunctionsNeighbor);
    evaluate(*local_diffusion_factor_entity,
             *local_diffusion_tensor_entity,
             *local_diffusion_factor_neighbor,
             *local_diffusion_tensor_neighbor,
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
  /// \name Actual Typelementation of order.
  /// \{

  template <class R, size_t rDF, size_t rCDF, size_t rDT, size_t rCDT, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t
  order(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDF, rCDF>&
            localDiffusionFactorEntity,
        const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDT, rCDT>&
            localDiffusionTensorEntity,
        const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDF, rCDF>&
            localDiffusionFactorNeighbor,
        const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDT, rCDT>&
            localDiffusionTensorNeighbor,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBaseEntity,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBaseEntity,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBaseNeighbor,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBaseNeighbor)
      const
  {
    return std::max(localDiffusionFactorEntity.order(), localDiffusionFactorNeighbor.order())
           + std::max(localDiffusionTensorEntity.order(), localDiffusionTensorNeighbor.order())
           + std::max(testBaseEntity.order(), testBaseNeighbor.order())
           + std::max(ansatzBaseEntity.order(), ansatzBaseNeighbor.order());
  } // size_t order(...)

  /// \}
  /// \name Actual Typelementation of evaluate.
  /// \{

  template <class R, class IntersectionType>
  void evaluate(
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localDiffusionFactorEntity,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain, dimDomain>&
          localDiffusionTensorEntity,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>&
          localDiffusionFactorNeighbor,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain, dimDomain>&
          localDiffusionTensorNeighbor,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& testBaseEntity,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& ansatzBaseEntity,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& testBaseNeighbor,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& ansatzBaseNeighbor,
      const IntersectionType& intersection, const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint,
      Dune::DynamicMatrix<R>& entityEntityRet, Dune::DynamicMatrix<R>& neighborNeighborRet,
      Dune::DynamicMatrix<R>& entityNeighborRet, Dune::DynamicMatrix<R>& neighborEntityRet) const
  {
    typedef Stuff::Common::FieldMatrix<R, dimDomain, dimDomain> TensorType;
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
    const auto local_diffusion_factor_en       = localDiffusionFactorEntity.evaluate(localPointEn);
    const TensorType local_diffusion_tensor_en = localDiffusionTensorEntity.evaluate(localPointEn);
    const auto local_diffusion_factor_ne       = localDiffusionFactorNeighbor.evaluate(localPointNe);
    const TensorType local_diffusion_tensor_ne = localDiffusionTensorNeighbor.evaluate(localPointNe);
    // compute penalty factor (see Epshteyn, Riviere, 2007)
    const size_t max_polorder =
        std::max(testBaseEntity.order(),
                 std::max(ansatzBaseEntity.order(), std::max(testBaseNeighbor.order(), ansatzBaseNeighbor.order())));
    const R sigma = SIPDG::internal::inner_sigma(max_polorder);
    // compute weighting (see Ern, Stephansen, Zunino 2007)
    // this evaluation has to be linear wrt the diffusion factor, so no other averaging method is allowed here!
    const auto local_diffusion_factor = (local_diffusion_factor_en + local_diffusion_factor_ne) * 0.5;
    const R delta_plus                = unitOuterNormal * (local_diffusion_tensor_ne * unitOuterNormal);
    const R delta_minus               = unitOuterNormal * (local_diffusion_tensor_en * unitOuterNormal);
    const R gamma                     = (delta_plus * delta_minus) / (delta_plus + delta_minus);
    const R penalty                   = (local_diffusion_factor * sigma * gamma) / std::pow(intersection.geometry().volume(), beta_);
    const R weight_plus               = delta_minus / (delta_plus + delta_minus);
    const R weight_minus              = delta_plus / (delta_plus + delta_minus);
    // compute diffusion value (should be factor * tensor, but this is the same)
    const auto diffusion_value_en = local_diffusion_tensor_en * local_diffusion_factor_en;
    const auto diffusion_value_ne = local_diffusion_tensor_ne * local_diffusion_factor_ne;
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
    assert(entityEntityRet.rows() >= rowsEn && entityEntityRet.cols() >= colsEn);
    assert(entityNeighborRet.rows() >= rowsEn && entityNeighborRet.cols() >= colsNe);
    assert(neighborEntityRet.rows() >= rowsNe && neighborEntityRet.cols() >= colsEn);
    assert(neighborNeighborRet.rows() >= rowsNe && neighborNeighborRet.cols() >= colsNe);
    // loop over all entity test basis functions
    for (size_t ii = 0; ii < rowsEn; ++ii) {
      auto& entityEntityRetRow   = entityEntityRet[ii];
      auto& entityNeighborRetRow = entityNeighborRet[ii];
      // loop over all entity ansatz basis functions
      for (size_t jj = 0; jj < colsEn; ++jj) {
        // consistency term
        entityEntityRetRow[jj] +=
            -weight_minus * ((diffusion_value_en * ansatzGradientsEn[jj][0]) * unitOuterNormal) * testValuesEn[ii];
        // symmetry term
        entityEntityRetRow[jj] +=
            -weight_minus * ansatzValuesEn[jj] * ((diffusion_value_en * testGradientsEn[ii][0]) * unitOuterNormal);
        // penalty term
        entityEntityRetRow[jj] += penalty * ansatzValuesEn[jj] * testValuesEn[ii];
      } // loop over all entity ansatz basis functions
      // loop over all neighbor ansatz basis functions
      for (size_t jj = 0; jj < colsNe; ++jj) {
        // consistency term
        entityNeighborRetRow[jj] +=
            -weight_plus * ((diffusion_value_ne * ansatzGradientsNe[jj][0]) * unitOuterNormal) * testValuesEn[ii];
        // symmetry term
        entityNeighborRetRow[jj] +=
            weight_minus * ansatzValuesNe[jj] * ((diffusion_value_en * testGradientsEn[ii][0]) * unitOuterNormal);
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
            weight_minus * ((diffusion_value_en * ansatzGradientsEn[jj][0]) * unitOuterNormal) * testValuesNe[ii];
        // symmetry term
        neighborEntityRetRow[jj] +=
            -weight_plus * ansatzValuesEn[jj] * ((diffusion_value_ne * testGradientsNe[ii][0]) * unitOuterNormal);
        // penalty term
        neighborEntityRetRow[jj] += -1.0 * penalty * ansatzValuesEn[jj] * testValuesNe[ii];
      } // loop over all entity ansatz basis functions
      // loop over all neighbor ansatz basis functions
      for (size_t jj = 0; jj < colsNe; ++jj) {
        // consistency term
        neighborNeighborRetRow[jj] +=
            weight_plus * ((diffusion_value_ne * ansatzGradientsNe[jj][0]) * unitOuterNormal) * testValuesNe[ii];
        // symmetry term
        neighborNeighborRetRow[jj] +=
            weight_plus * ansatzValuesNe[jj] * ((diffusion_value_ne * testGradientsNe[ii][0]) * unitOuterNormal);
        // penalty term
        neighborNeighborRetRow[jj] += penalty * ansatzValuesNe[jj] * testValuesNe[ii];
      } // loop over all neighbor ansatz basis functions
    } // loop over all neighbor test basis functions
  } // ... evaluate(...)

  /// \}

private:
  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
  const double beta_;
}; // Inner


template <class DiffusionFactorType, class DiffusionTensorType>
class InnerPenalty
    : public LocalEvaluation::Codim1Interface<internal::InnerPenaltyTraits<DiffusionFactorType, DiffusionTensorType>, 4>
{
public:
  typedef internal::InnerPenaltyTraits<DiffusionFactorType, DiffusionTensorType> Traits;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const size_t dimDomain = Traits::dimDomain;

  InnerPenalty(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& inducingFunction,
               const double beta = SIPDG::internal::default_beta(dimDomain))
    : diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(inducingFunction)
    , beta_(beta)
  {
  }

  /// \name Required by LocalEvaluation::Codim1Interface< ..., 4 >.
  /// \{

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_factor_.local_function(entity), diffusion_tensor_.local_function(entity));
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
    const auto local_diffusion_factor_entity   = std::get<0>(localFunctionsEntity);
    const auto local_diffusion_tensor_entity   = std::get<1>(localFunctionsEntity);
    const auto local_diffusion_factor_neighbor = std::get<0>(localFunctionsNeighbor);
    const auto local_diffusion_tensor_neighbor = std::get<1>(localFunctionsNeighbor);
    return order(*local_diffusion_factor_entity,
                 *local_diffusion_tensor_entity,
                 *local_diffusion_factor_neighbor,
                 *local_diffusion_tensor_neighbor,
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
    const auto local_diffusion_factor_entity   = std::get<0>(localFunctionsEntity);
    const auto local_diffusion_tensor_entity   = std::get<1>(localFunctionsEntity);
    const auto local_diffusion_factor_neighbor = std::get<0>(localFunctionsNeighbor);
    const auto local_diffusion_tensor_neighbor = std::get<1>(localFunctionsNeighbor);
    evaluate(*local_diffusion_factor_entity,
             *local_diffusion_tensor_entity,
             *local_diffusion_factor_neighbor,
             *local_diffusion_tensor_neighbor,
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
  /// \name Actual implementation of order.
  /// \{

  template <class R, size_t rDF, size_t rCDF, size_t rDT, size_t rCDT, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t
  order(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDF, rCDF>&
            localDiffusionFactorEntity,
        const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDT, rCDT>&
            localDiffusionTensorEntity,
        const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDF, rCDF>&
            localDiffusionFactorNeighbor,
        const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDT, rCDT>&
            localDiffusionTensorNeighbor,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBaseEntity,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBaseEntity,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBaseNeighbor,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBaseNeighbor)
      const
  {
    return std::max(localDiffusionFactorEntity.order(), localDiffusionFactorNeighbor.order())
           + std::max(localDiffusionTensorEntity.order(), localDiffusionTensorNeighbor.order())
           + std::max(testBaseEntity.order(), testBaseNeighbor.order())
           + std::max(ansatzBaseEntity.order(), ansatzBaseNeighbor.order());
  } // size_t order(...)

  /// \}
  /// \name Actual implementation of evaluation.
  /// \{

  template <class R, class IntersectionType>
  void evaluate(
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localDiffusionFactorEntity,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain, dimDomain>&
          localDiffusionTensorEntity,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>&
          localDiffusionFactorNeighbor,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain, dimDomain>&
          localDiffusionTensorNeighbor,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& testBaseEntity,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& ansatzBaseEntity,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& testBaseNeighbor,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& ansatzBaseNeighbor,
      const IntersectionType& intersection, const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint,
      Dune::DynamicMatrix<R>& entityEntityRet, Dune::DynamicMatrix<R>& neighborNeighborRet,
      Dune::DynamicMatrix<R>& entityNeighborRet, Dune::DynamicMatrix<R>& neighborEntityRet) const
  {
    typedef Stuff::Common::FieldMatrix<R, dimDomain, dimDomain> TensorType;
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
    const auto local_diffusion_factor_en       = localDiffusionFactorEntity.evaluate(localPointEn);
    const TensorType local_diffusion_tensor_en = localDiffusionTensorEntity.evaluate(localPointEn);
    const auto local_diffusion_factor_ne       = localDiffusionFactorNeighbor.evaluate(localPointNe);
    const TensorType local_diffusion_tensor_ne = localDiffusionTensorNeighbor.evaluate(localPointNe);
    // compute penalty factor (see Epshteyn, Riviere, 2007)
    const size_t max_polorder =
        std::max(testBaseEntity.order(),
                 std::max(ansatzBaseEntity.order(), std::max(testBaseNeighbor.order(), ansatzBaseNeighbor.order())));
    const R sigma = SIPDG::internal::inner_sigma(max_polorder);
    // compute weighting (see Ern, Stephansen, Zunino 2007)
    // this evaluation has to be linear wrt the diffusion factor, so no other averaging method is allowed here!
    const auto local_diffusion_factor = (local_diffusion_factor_en + local_diffusion_factor_ne) * 0.5;
    const R delta_plus                = unitOuterNormal * (local_diffusion_tensor_ne * unitOuterNormal);
    const R delta_minus               = unitOuterNormal * (local_diffusion_tensor_en * unitOuterNormal);
    const R gamma                     = (delta_plus * delta_minus) / (delta_plus + delta_minus);
    const R penalty                   = (local_diffusion_factor * sigma * gamma) / std::pow(intersection.geometry().volume(), beta_);
    // evaluate bases
    // * entity
    //   * test
    const size_t rowsEn     = testBaseEntity.size();
    const auto testValuesEn = testBaseEntity.evaluate(localPointEn);
    //   * ansatz
    const size_t colsEn       = ansatzBaseEntity.size();
    const auto ansatzValuesEn = ansatzBaseEntity.evaluate(localPointEn);
    // * neighbor
    //   * test
    const size_t rowsNe     = testBaseNeighbor.size();
    const auto testValuesNe = testBaseNeighbor.evaluate(localPointNe);
    //   * ansatz
    const size_t colsNe       = ansatzBaseNeighbor.size();
    const auto ansatzValuesNe = ansatzBaseNeighbor.evaluate(localPointNe);
    // compute the evaluations
    assert(entityEntityRet.rows() >= rowsEn && entityEntityRet.cols() >= colsEn);
    assert(entityNeighborRet.rows() >= rowsEn && entityNeighborRet.cols() >= colsNe);
    assert(neighborEntityRet.rows() >= rowsNe && neighborEntityRet.cols() >= colsEn);
    assert(neighborNeighborRet.rows() >= rowsNe && neighborNeighborRet.cols() >= colsNe);
    // loop over all entity test basis functions
    for (size_t ii = 0; ii < rowsEn; ++ii) {
      auto& entityEntityRetRow   = entityEntityRet[ii];
      auto& entityNeighborRetRow = entityNeighborRet[ii];
      // loop over all entity ansatz basis functions
      for (size_t jj = 0; jj < colsEn; ++jj) {
        // penalty term
        entityEntityRetRow[jj] += penalty * ansatzValuesEn[jj] * testValuesEn[ii];
      } // loop over all entity ansatz basis functions
      // loop over all neighbor ansatz basis functions
      for (size_t jj = 0; jj < colsNe; ++jj) {
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
        // penalty term
        neighborEntityRetRow[jj] += -1.0 * penalty * ansatzValuesEn[jj] * testValuesNe[ii];
      } // loop over all entity ansatz basis functions
      // loop over all neighbor ansatz basis functions
      for (size_t jj = 0; jj < colsNe; ++jj) {
        // penalty term
        neighborNeighborRetRow[jj] += penalty * ansatzValuesNe[jj] * testValuesNe[ii];
      } // loop over all neighbor ansatz basis functions
    } // loop over all neighbor test basis functions
  } // ... evaluate(...)

  /// \}

private:
  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
  const double beta_;
}; // InnerPenalty


template <class DiffusionFactorType, class DiffusionTensorType>
class BoundaryLHS
    : public LocalEvaluation::Codim1Interface<internal::BoundaryLHSTraits<DiffusionFactorType, DiffusionTensorType>, 2>
{
public:
  typedef internal::BoundaryLHSTraits<DiffusionFactorType, DiffusionTensorType> Traits;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const size_t dimDomain = Traits::dimDomain;

  BoundaryLHS(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
              const double beta = SIPDG::internal::default_beta(dimDomain))
    : diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , beta_(beta)
  {
  }

  /// \name Required by LocalEvaluation::Codim1Interface< ..., 2 >

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_factor_.local_function(entity), diffusion_tensor_.local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t
  order(const LocalfunctionTupleType& localFuncs,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase) const
  {
    const auto local_diffusion_factor = std::get<0>(localFuncs);
    const auto local_diffusion_tensor = std::get<1>(localFuncs);
    return order(*local_diffusion_factor, *local_diffusion_tensor, testBase, ansatzBase);
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class IntersectionType, class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  void evaluate(const LocalfunctionTupleType& localFuncs,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase,
                const IntersectionType& intersection,
                const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    const auto local_diffusion_factor = std::get<0>(localFuncs);
    const auto local_diffusion_tensor = std::get<1>(localFuncs);
    evaluate(*local_diffusion_factor, *local_diffusion_tensor, testBase, ansatzBase, intersection, localPoint, ret);
  }

  /// \}
  /// \name Actual implementation of order.
  /// \{

  template <class R, size_t rDF, size_t rCDF, size_t rDT, size_t rCDT, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t
  order(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDF, rCDF>& localDiffusionFactor,
        const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDT, rCDT>& localDiffusionTensor,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase) const
  {
    return localDiffusionFactor.order() + localDiffusionTensor.order() + testBase.order() + ansatzBase.order();
  }

  /// \}
  /// \name Actual implementation of evaluate.
  /// \{

  template <class R, class IntersectionType>
  void
  evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localDiffusionFactor,
           const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain, dimDomain>&
               localDiffusionTensor,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& testBase,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& ansatzBase,
           const IntersectionType& intersection, const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint,
           Dune::DynamicMatrix<R>& ret) const
  {
    typedef Stuff::Common::FieldMatrix<R, dimDomain, dimDomain> TensorType;
    // clear ret
    ret *= 0.0;
    // get local point (which is in intersection coordinates) in entity coordinates
    const auto localPointEntity = intersection.geometryInInside().global(localPoint);
    const auto unitOuterNormal  = intersection.unitOuterNormal(localPoint);
    // evaluate local function
    const auto diffusion_factor_value       = localDiffusionFactor.evaluate(localPointEntity);
    const TensorType diffusion_tensor_value = localDiffusionTensor.evaluate(localPointEntity);
    // compute penalty (see Epshteyn, Riviere, 2007)
    const size_t max_polorder = std::max(testBase.order(), ansatzBase.order());
    const R sigma             = SIPDG::internal::boundary_sigma(max_polorder);
    // compute weighting (see Ern, Stephansen, Zunino 2007)
    const R gamma   = unitOuterNormal * (diffusion_tensor_value * unitOuterNormal);
    const R penalty = (diffusion_factor_value * sigma * gamma) / std::pow(intersection.geometry().volume(), beta_);
    // compute diffusion value (should be factor * tensor, but this is the same)
    const auto diffusion_value = diffusion_tensor_value * diffusion_factor_value;
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
    assert(ret.rows() >= rows && ret.cols() >= cols);
    // loop over all test basis functions
    for (size_t ii = 0; ii < rows; ++ii) {
      auto& retRow = ret[ii];
      // loop over all ansatz basis functions
      for (size_t jj = 0; jj < cols; ++jj) {
        // consistency term
        retRow[jj] += -1.0 * ((diffusion_value * ansatzGradients[jj][0]) * unitOuterNormal) * testValues[ii];
        // symmetry term
        retRow[jj] += -1.0 * ansatzValues[jj] * ((diffusion_value * testGradients[ii][0]) * unitOuterNormal);
        // penalty term
        retRow[jj] += penalty * ansatzValues[jj] * testValues[ii];
      } // loop over all ansatz basis functions
    } // loop over all test basis functions
  } // void evaluate(...)

  /// \}

private:
  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
  const double beta_;
}; // class BoundaryLHS


template <class DiffusionFactorType, class DiffusionTensorType>
class BoundaryLHSPenalty
    : public LocalEvaluation::
          Codim1Interface<internal::BoundaryLHSPenaltyTraits<DiffusionFactorType, DiffusionTensorType>, 2>
{
public:
  typedef internal::BoundaryLHSPenaltyTraits<DiffusionFactorType, DiffusionTensorType> Traits;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const size_t dimDomain = Traits::dimDomain;

  BoundaryLHSPenalty(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
                     const double beta = SIPDG::internal::default_beta(dimDomain))
    : diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , beta_(beta)
  {
  }

  /// \name Required by LocalEvaluation::Codim1Interface< ..., 2 >
  /// \{

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_factor_.local_function(entity), diffusion_tensor_.local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t
  order(const LocalfunctionTupleType& localFuncs,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase) const
  {
    const auto local_diffusion_factor = std::get<0>(localFuncs);
    const auto local_diffusion_tensor = std::get<1>(localFuncs);
    return order(*local_diffusion_factor, *local_diffusion_tensor, testBase, ansatzBase);
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class IntersectionType, class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  void evaluate(const LocalfunctionTupleType& localFuncs,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase,
                const IntersectionType& intersection,
                const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    const auto local_diffusion_factor = std::get<0>(localFuncs);
    const auto local_diffusion_tensor = std::get<1>(localFuncs);
    evaluate(*local_diffusion_factor, *local_diffusion_tensor, testBase, ansatzBase, intersection, localPoint, ret);
  }

  /// \}
  /// \name Actual implementation of order.
  /// \{

  template <class R, size_t rDF, size_t rCDF, size_t rDT, size_t rCDT, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t
  order(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDF, rCDF>& localDiffusionFactor,
        const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDT, rCDT>& localDiffusionTensor,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase) const
  {
    return localDiffusionFactor.order() + localDiffusionTensor.order() + testBase.order() + ansatzBase.order();
  }

  /// \}
  /// \name Actual implementation of evaluate.
  /// \{

  template <class R, class IntersectionType>
  void
  evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localDiffusionFactor,
           const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain, dimDomain>&
               localDiffusionTensor,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& testBase,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& ansatzBase,
           const IntersectionType& intersection, const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint,
           Dune::DynamicMatrix<R>& ret) const
  {
    typedef Stuff::Common::FieldMatrix<R, dimDomain, dimDomain> TensorType;
#ifndef NDEBUG
    auto logger = DSC::TimedLogger().get("gdt.localevaluation.swipdg.boundarylhspenalty");
#endif
    // clear ret
    ret *= 0.0;
    // get local point (which is in intersection coordinates) in entity coordinates
    const auto localPointEntity = intersection.geometryInInside().global(localPoint);
    const auto unitOuterNormal  = intersection.unitOuterNormal(localPoint);
#ifndef NDEBUG
    DSC::print(intersection.geometry().global(localPoint), "global(localPoint)", logger.debug());
    DSC::print(localPointEntity, "localPointEntity", logger.debug(), "  ");
    DSC::print(unitOuterNormal, "unitOuterNormal", logger.debug(), "  ");
#endif // NDEBUG
    // evaluate local function
    const auto diffusion_factor_value       = localDiffusionFactor.evaluate(localPointEntity);
    const TensorType diffusion_tensor_value = localDiffusionTensor.evaluate(localPointEntity);
#ifndef NDEBUG
    DSC::print(diffusion_factor_value, "diffusion_factor_value", logger.debug(), "  ");
    DSC::print(diffusion_tensor_value, "diffusion_tensor_value", logger.debug(), "  ");
#endif // NDEBUG
    // compute penalty (see Epshteyn, Riviere, 2007)
    const size_t max_polorder = std::max(testBase.order(), ansatzBase.order());
    const R sigma             = SIPDG::internal::boundary_sigma(max_polorder);
    // compute weighting (see Ern, Stephansen, Zunino 2007)
    const R gamma   = unitOuterNormal * (diffusion_tensor_value * unitOuterNormal);
    const R penalty = (diffusion_factor_value * sigma * gamma) / std::pow(intersection.geometry().volume(), beta_);
#ifndef NDEBUG
    logger.debug() << "  max_polorder = " << max_polorder << std::endl;
    logger.debug() << "  sigma = " << sigma << std::endl;
    logger.debug() << "  gamma = " << gamma << std::endl;
    logger.debug() << "  intersection.geometry().volume() = " << intersection.geometry().volume() << std::endl;
    logger.debug() << "  beta_ = " << beta_ << std::endl;
    logger.debug() << "  penalty = " << penalty << std::endl;
#endif // NDEBUG
    // evaluate bases
    // * test
    const size_t rows     = testBase.size();
    const auto testValues = testBase.evaluate(localPointEntity);
    // * ansatz
    const size_t cols       = ansatzBase.size();
    const auto ansatzValues = ansatzBase.evaluate(localPointEntity);
    // compute products
    assert(ret.rows() >= rows && ret.cols() >= cols);
    // loop over all test basis functions
    for (size_t ii = 0; ii < rows; ++ii) {
      auto& retRow = ret[ii];
      // loop over all ansatz basis functions
      for (size_t jj = 0; jj < cols; ++jj) {
        // penalty term
        retRow[jj] += penalty * ansatzValues[jj] * testValues[ii];
      } // loop over all ansatz basis functions
    } // loop over all test basis functions
#ifndef NDEBUG
    DSC::print(testValues, "testValues", logger.debug(), "  ");
    DSC::print(ansatzValues, "ansatzValues", logger.debug(), "  ");
    DSC::print(ret, "ret", logger.debug(), "  ");
#endif // NDEBUG
  } // void evaluate(...) const

  /// \}

private:
  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
  const double beta_;
}; // class BoundaryLHSPenalty


template <class DiffusionFactorType, class DirichletType, class DiffusionTensorType>
class BoundaryRHS
    : public LocalEvaluation::
          Codim1Interface<internal::BoundaryRHSTraits<DiffusionFactorType, DirichletType, DiffusionTensorType>, 1>
{
public:
  typedef internal::BoundaryRHSTraits<DiffusionFactorType, DirichletType, DiffusionTensorType> Traits;
  typedef typename Traits::LocalizableDirichletFunctionType LocalizableDirichletFunctionType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const size_t dimDomain = Traits::dimDomain;

  BoundaryRHS(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
              const LocalizableDirichletFunctionType& dirichlet,
              const double beta = SIPDG::internal::default_beta(dimDomain))
    : diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , dirichlet_(dirichlet)
    , beta_(beta)
  {
  }

  /// \name Required by LocalEvaluation::Codim1Interface< ..., 1 >.
  /// \{

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_factor_.local_function(entity),
                           diffusion_tensor_.local_function(entity),
                           dirichlet_.local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class R, size_t r, size_t rC>
  size_t order(const LocalfunctionTupleType& localFuncs,
               const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, rC>& testBase) const
  {
    const auto localDiffusionFactor = std::get<0>(localFuncs);
    const auto localDiffusionTensor = std::get<1>(localFuncs);
    const auto localDirichlet = std::get<2>(localFuncs);
    return order(*localDiffusionFactor, *localDiffusionTensor, *localDirichlet, testBase);
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class IntersectionType, class R, size_t r, size_t rC>
  void evaluate(const LocalfunctionTupleType& localFuncs,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, rC>& testBase,
                const IntersectionType& intersection,
                const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint, Dune::DynamicVector<R>& ret) const
  {
    const auto localDiffusionFactor = std::get<0>(localFuncs);
    const auto localDiffusionTensor = std::get<1>(localFuncs);
    const auto localDirichlet = std::get<2>(localFuncs);
    evaluate(*localDiffusionFactor, *localDiffusionTensor, *localDirichlet, testBase, intersection, localPoint, ret);
  }

  /// \}
  /// \name Actual implementation of order.
  /// \{

  template <class R, size_t rDF, size_t rCDF, size_t rDT, size_t rCDT, size_t rLR, size_t rCLR, size_t rT, size_t rCT>
  size_t
  order(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDF, rCDF>& localDiffusionFactor,
        const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDT, rCDT>& localDiffusionTensor,
        const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rLR, rCLR>& localDirichlet,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase) const
  {
    const size_t testOrder         = testBase.order();
    const size_t testGradientOrder = std::max(ssize_t(testOrder) - 1, ssize_t(0));
    const size_t diffusionOrder    = localDiffusionFactor.order() + localDiffusionTensor.order();
    const size_t dirichletOrder = localDirichlet.order();
    return std::max(testOrder + dirichletOrder, diffusionOrder + testGradientOrder + dirichletOrder);
  } // ... order(...)

  /// \}
  /// \name Actual implementation of evaluate.
  /// \{

  template <class R, class IntersectionType>
  void
  evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localDiffusionFactor,
           const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain, dimDomain>&
               localDiffusionTensor,
           const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localDirichlet,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& testBase,
           const IntersectionType& intersection, const Dune::FieldVector<DomainFieldType, dimDomain - 1>& localPoint,
           Dune::DynamicVector<R>& ret) const
  {
    typedef Stuff::Common::FieldMatrix<R, dimDomain, dimDomain> TensorType;
    // clear ret
    ret *= 0.0;
    // get local point (which is in intersection coordinates) in entity coordinates
    const auto localPointEntity = intersection.geometryInInside().global(localPoint);
    const auto unitOuterNormal  = intersection.unitOuterNormal(localPoint);
    // evaluate local functions
    const auto diffusionFactorValue       = localDiffusionFactor.evaluate(localPointEntity);
    const TensorType diffusionTensorValue = localDiffusionTensor.evaluate(localPointEntity);
    const auto dirichletValue             = localDirichlet.evaluate(localPointEntity);
    // compute penalty (see Epshteyn, Riviere, 2007)
    const size_t polorder = testBase.order();
    const R sigma         = SIPDG::internal::boundary_sigma(polorder);
    // compute weighting (see Ern, Stephansen, Zunino 2007)
    const R gamma   = unitOuterNormal * (diffusionTensorValue * unitOuterNormal);
    const R penalty = (diffusionFactorValue * sigma * gamma) / std::pow(intersection.geometry().volume(), beta_);
    // compute diffusion value (should be factor * tensor, but this is the same)
    const auto diffusionValue = diffusionTensorValue * diffusionFactorValue;
    // evaluate basis
    const size_t size        = testBase.size();
    const auto testValues    = testBase.evaluate(localPointEntity);
    const auto testGradients = testBase.jacobian(localPointEntity);
    // compute
    assert(ret.size() >= size);
    // loop over all test basis functions
    for (size_t ii = 0; ii < size; ++ii) {
      // symmetry term
      ret[ii] += -1.0 * dirichletValue * ((diffusionValue * testGradients[ii][0]) * unitOuterNormal);
      // penalty term
      ret[ii] += penalty * dirichletValue * testValues[ii];
    } // loop over all test basis functions
  } // ... evaluate(...)

  /// \}

  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
  const LocalizableDirichletFunctionType& dirichlet_;
  const double beta_;
}; // class BoundaryRHS


} // namespace SWIPDG
} // namespace LocalEvaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_LOCALEVALUATION_SWIPDG_HH
