// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Kirsten Weber

#ifndef DUNE_GDT_EVALUATION_ELLIPTIC_HH
#define DUNE_GDT_EVALUATION_ELLIPTIC_HH

#include <tuple>

#include <dune/stuff/common/disable_warnings.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/typetraits.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/functions/interfaces.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace LocalEvaluation {


// forwards
template <class DiffusionFactorImp, class DiffusionTensorImp = void>
class Elliptic;


namespace internal {


/**
 *  \brief  Traits for the Elliptic evaluation.
 */
template <class DiffusionFactorImp, class DiffusionTensorImp>
class EllipticTraits
{
  static_assert(std::is_base_of<Dune::Stuff::IsLocalizableFunction, DiffusionFactorImp>::value,
                "DiffusionFactorImp has to be derived from Stuff::IsLocalizableFunction.");
  static_assert(std::is_base_of<Dune::Stuff::IsLocalizableFunction, DiffusionTensorImp>::value,
                "DiffusionTensorImp has to be derived from Stuff::IsLocalizableFunction.");
  static_assert(std::is_same<typename DiffusionFactorImp::EntityType, typename DiffusionTensorImp::EntityType>::value,
                "EntityTypes have to agree!");
  static_assert(
      std::is_same<typename DiffusionFactorImp::DomainFieldType, typename DiffusionTensorImp::DomainFieldType>::value,
      "DomainFieldTypes have to agree!");
  static_assert(DiffusionFactorImp::dimDomain == DiffusionTensorImp::dimDomain, "Dimensions of domains have to agree");

public:
  typedef Elliptic<DiffusionFactorImp, DiffusionTensorImp> derived_type;
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
class EllipticTraits<LocalizableFunctionImp, void>
{
  static_assert(std::is_base_of<Dune::Stuff::IsLocalizableFunction, LocalizableFunctionImp>::value,
                "LocalizableFunctionImp has to be derived from Stuff::IsLocalizableFunction.");

public:
  typedef LocalizableFunctionImp LocalizableFunctionType;
  typedef Elliptic<LocalizableFunctionType, void> derived_type;
  typedef typename LocalizableFunctionType::EntityType EntityType;
  typedef typename LocalizableFunctionType::DomainFieldType DomainFieldType;
  typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;
  typedef std::tuple<std::shared_ptr<LocalfunctionType>> LocalfunctionTupleType;
  static const unsigned int dimDomain = LocalizableFunctionType::dimDomain;
};
} // namespace internal


/**
 *  \brief  Computes an elliptic evaluation.
 */
template <class LocalizableFunctionImp>
class Elliptic<LocalizableFunctionImp, void>
    : public LocalEvaluation::Codim0Interface<internal::EllipticTraits<LocalizableFunctionImp, void>, 2>
{
public:
  typedef internal::EllipticTraits<LocalizableFunctionImp, void> Traits;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = Traits::dimDomain;

  explicit Elliptic(const LocalizableFunctionType& inducingFunction)
    : inducingFunction_(inducingFunction)
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
  order(const LocalfunctionTupleType& localFuncs,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase) const
  {
    const auto localFunction = std::get<0>(localFuncs);
    return redirect_order(*localFunction, testBase, ansatzBase);
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class R, int rT, int rCT, int rA, int rCA>
  void evaluate(const LocalfunctionTupleType& localFuncs,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase,
                const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    const auto localFunction = std::get<0>(localFuncs);
    evaluate(*localFunction, testBase, ansatzBase, localPoint, ret);
  }

  /**
   *  \brief  Computes an elliptic evaluation for a 2x2 matrix-valued local function and matrix-valued basefunctionsets.
   *  \tparam R RangeFieldType
   *  \note   Unfortunately we need this explicit specialization, otherwise the compiler will complain for 1d grids.
   */
  template <class R>
  void evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 2, 2>& localFunction,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& testBase,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& ansatzBase,
                const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    evaluate_matrix_valued_(localFunction, testBase, ansatzBase, localPoint, ret);
  }

  /**
   *  \brief  Computes an elliptic evaluation for a 3x3 matrix-valued local function and matrix-valued basefunctionsets.
   *  \tparam R RangeFieldType
   *  \note   Unfortunately we need this explicit specialization, otherwise the compiler will complain for 1d grids.
   */
  template <class R>
  void evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 3, 3>& localFunction,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& testBase,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& ansatzBase,
                const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    evaluate_matrix_valued_(localFunction, testBase, ansatzBase, localPoint, ret);
  }

private:
  /**
   *  \return localFunction.order() + (testBase.order() - 1) + (ansatzBase.order() - 1)
   */
  template <class R, int rL, int rCL, int rT, int rCT, int rA, int rCA>
  size_t redirect_order(
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& localFunction,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase) const
  {
    return localFunction.order() + std::max(ssize_t(testBase.order()) - 1, ssize_t(0))
           + std::max(ssize_t(ansatzBase.order()) - 1, ssize_t(0));
  } // size_t redirect_order( ... )

  //  template< class R, int rL, int rCL, int rT, int rCT, int rA, int rCA >
  //  void redirect_evaluate(const Stuff::LocalfunctionInterface
  //                             < EntityType, DomainFieldType, dimDomain, R, rL, rCL >& /*localFunction*/,
  //                         const Stuff::LocalfunctionSetInterface
  //                             < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& /*testBase*/,
  //                         const Stuff::LocalfunctionSetInterface
  //                             < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& /*ansatzBase*/,
  //                         const Dune::FieldVector< DomainFieldType, dimDomain >& /*localPoint*/,
  //                         Dune::DynamicMatrix< R >& /*ret*/) const
  //  {
  //    static_assert(Dune::AlwaysFalse< R >::value, "Not implemented for these dimensions!");
  //  } // ... redirect_evaluate< ... >(...)

  /**
   *  \brief  Computes an elliptic evaluation for a scalar local function and scalar or vector valued basefunctionsets.
   *  \tparam R RangeFieldType
   */
  template <class R, int r>
  void evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& localFunction,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& testBase,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& ansatzBase,
                const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    typedef
        typename Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>::JacobianRangeType
            JacobianRangeType;
    // evaluate local function
    const auto functionValue = localFunction.evaluate(localPoint);
    // evaluate test gradient
    const size_t rows = testBase.size();
    std::vector<JacobianRangeType> testGradients(rows, JacobianRangeType(0));
    testBase.jacobian(localPoint, testGradients);
    // evaluate ansatz gradient
    const size_t cols = ansatzBase.size();
    std::vector<JacobianRangeType> ansatzGradients(cols, JacobianRangeType(0));
    ansatzBase.jacobian(localPoint, ansatzGradients);
    // compute products
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    for (size_t ii = 0; ii < rows; ++ii) {
      auto& retRow = ret[ii];
      for (size_t jj = 0; jj < cols; ++jj) {
        retRow[jj] = functionValue * (ansatzGradients[jj][0] * testGradients[ii][0]);
      }
    }
  } // ... redirect_evaluate< ..., 1, 1 >(...)

  template <class R>
  void evaluate_matrix_valued_(
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain, dimDomain>&
          localFunction,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& testBase,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& ansatzBase,
      const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    typedef typename Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain, dimDomain>::
        RangeType DiffusionRangeType;
    typedef
        typename Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>::JacobianRangeType
            JacobianRangeType;
    // evaluate local function
    const DiffusionRangeType functionValue = localFunction.evaluate(localPoint);
    // evaluate test gradient
    const size_t rows = testBase.size();
    std::vector<JacobianRangeType> testGradients(rows, JacobianRangeType(0));
    testBase.jacobian(localPoint, testGradients);
    // evaluate ansatz gradient
    const size_t cols = ansatzBase.size();
    std::vector<JacobianRangeType> ansatzGradients(cols, JacobianRangeType(0));
    ansatzBase.jacobian(localPoint, ansatzGradients);
    // compute products
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    FieldVector<DomainFieldType, dimDomain> product(0.0);
    for (size_t ii = 0; ii < rows; ++ii) {
      auto& retRow = ret[ii];
      for (size_t jj = 0; jj < cols; ++jj) {
        functionValue.mv(ansatzGradients[jj][0], product);
        retRow[jj] = product * testGradients[ii][0];
      }
    }
  } // ... evaluate_matrix_valued_(...)

  const LocalizableFunctionType& inducingFunction_;
}; // class LocalElliptic


} // namespace Evaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_EVALUATION_ELLIPTIC_HH
