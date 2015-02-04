// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Kirsten Weber

#ifndef DUNE_GDT_EVALUATION_ELLIPTIC_HH
#define DUNE_GDT_EVALUATION_ELLIPTIC_HH

#include <tuple>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/dynmatrix.hh>
#include <dune/common/typetraits.hh>

#include <dune/stuff/functions/interfaces.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace LocalEvaluation {


// forward
template <class DiffusionFactorImp, class DiffusionTensorImp = void>
class Elliptic;


namespace internal {


/**
 *  \brief  Traits for the Elliptic evaluation.
 */
template <class DiffusionFactorType, class DiffusionTensorType>
class EllipticTraits
{
  static_assert(Stuff::is_localizable_function<DiffusionFactorType>::value,
                "DiffusionFactorType has to be a localizable function!");
  static_assert(Stuff::is_localizable_function<DiffusionTensorType>::value,
                "DiffusionTensorType has to be a localizable function!");
  static_assert(std::is_same<typename DiffusionFactorType::EntityType, typename DiffusionTensorType::EntityType>::value,
                "EntityTypes have to agree!");
  static_assert(
      std::is_same<typename DiffusionFactorType::DomainFieldType, typename DiffusionTensorType::DomainFieldType>::value,
      "DomainFieldTypes have to agree!");
  static_assert(DiffusionFactorType::dimDomain == DiffusionTensorType::dimDomain, "Dimensions have to agree!");

public:
  typedef Elliptic<DiffusionFactorType, DiffusionTensorType> derived_type;
  typedef std::tuple<std::shared_ptr<typename DiffusionFactorType::LocalfunctionType>,
                     std::shared_ptr<typename DiffusionTensorType::LocalfunctionType>> LocalfunctionTupleType;
  typedef typename DiffusionFactorType::EntityType EntityType;
  typedef typename DiffusionFactorType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = DiffusionFactorType::dimDomain;
}; // class EllipticTraits


template <class DiffusionType>
class EllipticTraits<DiffusionType, void>
{
  static_assert(Stuff::is_localizable_function<DiffusionType>::value,
                "DiffusionType has to be a localizable function!");

public:
  typedef Elliptic<DiffusionType, void> derived_type;
  typedef typename DiffusionType::EntityType EntityType;
  typedef typename DiffusionType::DomainFieldType DomainFieldType;
  typedef std::tuple<std::shared_ptr<typename DiffusionType::LocalfunctionType>> LocalfunctionTupleType;
  static const unsigned int dimDomain = DiffusionType::dimDomain;
}; // class EllipticTraits< ..., void >


} // namespace internal


/**
 *  \brief  Computes an elliptic evaluation.
 */
template <class DiffusionImp>
class Elliptic<DiffusionImp, void>
    : public LocalEvaluation::Codim0Interface<internal::EllipticTraits<DiffusionImp, void>, 2>
{
public:
  typedef DiffusionImp DiffusionType;
  typedef internal::EllipticTraits<DiffusionType, void> Traits;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = Traits::dimDomain;

  explicit Elliptic(const DiffusionType& inducingFunction)
    : diffusion_(inducingFunction)
  {
  }

  /// \name Required by LocalEvaluation::Codim0Interface< ..., 2 >
  /// \{

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_.local_function(entity));
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
    return order(*std::get<0>(localFuncs), testBase, ansatzBase);
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
    evaluate(*std::get<0>(localFuncs), testBase, ansatzBase, localPoint, ret);
  }

  /// \}
  /// \name Actual implementations of order
  /// \{

  /**
   *  \return localFunction.order() + (testBase.order() - 1) + (ansatzBase.order() - 1)
   */
  template <class R, int rL, int rCL, int rT, int rCT, int rA, int rCA>
  size_t
  order(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rL, rCL>& localFunction,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase) const
  {
    return localFunction.order() + boost::numeric_cast<size_t>(std::max(ssize_t(testBase.order()) - 1, ssize_t(0)))
           + boost::numeric_cast<size_t>(std::max(ssize_t(ansatzBase.order()) - 1, ssize_t(0)));
  } // ... order( ... )

  /// \}
  /// \name Actual implementations of evaluate
  /// \{

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
  } // ... redirect_evaluate< ..., 1, ... >(...)

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

  /// \}

private:
  template <class R>
  void evaluate_matrix_valued_(
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain, dimDomain>&
          localFunction,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& testBase,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& ansatzBase,
      const Dune::FieldVector<DomainFieldType, dimDomain>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    // evaluate local function
    const auto functionValue = localFunction.evaluate(localPoint);
    // evaluate test gradient
    const size_t rows        = testBase.size();
    const auto testGradients = testBase.jacobian(localPoint);
    // evaluate ansatz gradient
    const size_t cols          = ansatzBase.size();
    const auto ansatzGradients = ansatzBase.jacobian(localPoint);
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

  const DiffusionType& diffusion_;
}; // class Elliptic< ...., void >


} // namespace Evaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_EVALUATION_ELLIPTIC_HH
