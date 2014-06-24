// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Kirsten Weber

#ifndef DUNE_GDT_EVALUATION_ELLIPTIC_HH
#define DUNE_GDT_EVALUATION_ELLIPTIC_HH

#include <tuple>
#include <type_traits>

#include <dune/stuff/common/disable_warnings.hh>
# include <dune/common/dynmatrix.hh>
#include <dune/stuff/common/reenable_warnings.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/static_assert.hh>

#include <dune/stuff/functions/interfaces.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace LocalEvaluation {


// forwards
template< class DiffusionFactorType, class DiffusionTensorType = void >
class Elliptic;


namespace internal {


/**
 *  \brief  Traits for the Elliptic evaluation.
 */
template< class DiffusionFactorType, class DiffusionTensorType >
class EllipticTraits
{
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, DiffusionFactorType >::value,
                "DiffusionFactorType has to be derived from Stuff::IsLocalizableFunction.");
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, DiffusionTensorType >::value,
                "DiffusionTensorType has to be derived from Stuff::IsLocalizableFunction.");
public:
  typedef Elliptic< DiffusionFactorType, DiffusionTensorType > derived_type;
};


template< class LocalizableFunctionType >
class EllipticTraits< LocalizableFunctionType, void >
{
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, LocalizableFunctionType >::value,
                "LocalizableFunctionType has to be derived from Stuff::IsLocalizableFunction.");
public:
  typedef Elliptic< LocalizableFunctionType, void > derived_type;
};
} // namespace internal


/**
 *  \brief  Computes an elliptic evaluation.
 */
template< class LocalizableFunctionType >
class Elliptic< LocalizableFunctionType, void >
  : public LocalEvaluation::Codim0Interface< internal::EllipticTraits< LocalizableFunctionType, void >, 2 >
{
  typedef typename LocalizableFunctionType::EntityType EntityType;
public:
  typedef internal::EllipticTraits< LocalizableFunctionType, void > Traits;

  Elliptic(const LocalizableFunctionType& inducingFunction)
    : inducingFunction_(inducingFunction)
  {}


  class LocalfunctionTuple
  {
    typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;
  public:
    typedef std::tuple< std::shared_ptr< LocalfunctionType > > Type;
  };


  typename LocalfunctionTuple::Type localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(inducingFunction_.local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template< class E, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  size_t order(const typename LocalfunctionTuple::Type& localFuncs,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase) const
  {
    const auto localFunction = std::get< 0 >(localFuncs);
    return order(*localFunction, testBase, ansatzBase);
  }

  /**
   *  \todo add copydoc
   *  \return localFunction.order() + (testBase.order() - 1) + (ansatzBase.order() - 1)
   */
  template< class E, class D, int d, class R, int rL, int rCL, int rT, int rCT, int rA, int rCA >
  size_t order(const Stuff::LocalfunctionInterface< E, D, d, R, rL, rCL >& localFunction,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase) const
  {
    return localFunction.order()
        + std::max(ssize_t(testBase.order()) - 1, ssize_t(0))
        + std::max(ssize_t(ansatzBase.order()) - 1, ssize_t(0));
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template< class E, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  void evaluate(const typename LocalfunctionTuple::Type& localFuncs,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase,
                const Dune::FieldVector< D, d >& localPoint,
                Dune::DynamicMatrix< R >& ret) const
  {
    const auto localFunction = std::get< 0 >(localFuncs);
    evaluate(*localFunction, testBase, ansatzBase, localPoint, ret);
  }

  template< class E, class D, int d, class R, int rL, int rCL, int rT, int rCT, int rA, int rCA >
  void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, rL, rCL >& /*localFunction*/,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& /*testBase*/,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& /*ansatzBase*/,
                const Dune::FieldVector< D, d >& /*localPoint*/,
                Dune::DynamicMatrix< R >& /*ret*/) const
  {
    static_assert(Dune::AlwaysFalse< R >::value, "Not implemented for these dimensions!");
  }

  /**
   *  \brief  Computes an elliptic evaluation for a scalar local function and scalar or vector valued basefunctionsets.
   *  \tparam E EntityType
   *  \tparam D DomainFieldType
   *  \tparam d dimDomain
   *  \tparam R RangeFieldType
   */
  template< class E, class D, int d, class R, int r >
  void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, 1, 1 >& localFunction,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, r, 1 >& testBase,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, r, 1 >& ansatzBase,
                const Dune::FieldVector< D, d >& localPoint,
                Dune::DynamicMatrix< R >& ret) const
  {
    typedef typename Stuff::LocalfunctionSetInterface< E, D, d, R, r, 1 >::JacobianRangeType JacobianRangeType;
    // evaluate local function
    const auto functionValue = localFunction.evaluate(localPoint);
    // evaluate test gradient
    const size_t rows = testBase.size();
    std::vector< JacobianRangeType > testGradients(rows, JacobianRangeType(0));
    testBase.jacobian(localPoint, testGradients);
    // evaluate ansatz gradient
    const size_t cols = ansatzBase.size();
    std::vector< JacobianRangeType > ansatzGradients(cols, JacobianRangeType(0));
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
  } // ... evaluate< ..., 1, 1 >(...)

  /**
   *  \brief  Computes an elliptic evaluation for a 2x2 matrix-valued local function and matrix-valued basefunctionsets.
   *  \tparam E EntityType
   *  \tparam D DomainFieldType
   *  \tparam d dimDomain
   *  \tparam R RangeFieldType
   *  \note   Unfortunately we need this explicit specialization, otherwise the compiler will complain for 1d grids.
   */
  template< class E, class D, int d, class R >
  void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, 2, 2 >& localFunction,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >& testBase,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >& ansatzBase,
                const Dune::FieldVector< D, d >& localPoint,
                Dune::DynamicMatrix< R >& ret) const
  {
    evaluate_matrix_valued_(localFunction, testBase, ansatzBase, localPoint, ret);
  }

  /**
   *  \brief  Computes an elliptic evaluation for a 3x3 matrix-valued local function and matrix-valued basefunctionsets.
   *  \tparam E EntityType
   *  \tparam D DomainFieldType
   *  \tparam d dimDomain
   *  \tparam R RangeFieldType
   *  \note   Unfortunately we need this explicit specialization, otherwise the compiler will complain for 1d grids.
   */
  template< class E, class D, int d, class R >
  void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, 3, 3 >& localFunction,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >& testBase,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >& ansatzBase,
                const Dune::FieldVector< D, d >& localPoint,
                Dune::DynamicMatrix< R >& ret) const
  {
    evaluate_matrix_valued_(localFunction, testBase, ansatzBase, localPoint, ret);
  }

private:
  template< class E, class D, int d, class R >
  void evaluate_matrix_valued_(const Stuff::LocalfunctionInterface< E, D, d, R, d, d >& localFunction,
                               const Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >& testBase,
                               const Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >& ansatzBase,
                               const Dune::FieldVector< D, d >& localPoint,
                               Dune::DynamicMatrix< R >& ret) const
  {
    typedef typename Stuff::LocalfunctionInterface< E, D, d, R, d, d >::RangeType             DiffusionRangeType;
    typedef typename Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >::JacobianRangeType  JacobianRangeType;
    // evaluate local function
    const DiffusionRangeType functionValue = localFunction.evaluate(localPoint);
    // evaluate test gradient
    const size_t rows = testBase.size();
    std::vector< JacobianRangeType > testGradients(rows, JacobianRangeType(0));
    testBase.jacobian(localPoint, testGradients);
    // evaluate ansatz gradient
    const size_t cols = ansatzBase.size();
    std::vector< JacobianRangeType > ansatzGradients(cols, JacobianRangeType(0));
    ansatzBase.jacobian(localPoint, ansatzGradients);
    // compute products
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    FieldVector< D, d > product(0.0);
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
