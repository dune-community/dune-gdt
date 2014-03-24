// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Kirsten Weber

#ifndef DUNE_GDT_EVALUATION_ELLIPTIC_HH
#define DUNE_GDT_EVALUATION_ELLIPTIC_HH

#include <tuple>
#include <type_traits>

#include <dune/common/dynmatrix.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/static_assert.hh>

#include <dune/stuff/functions/interfaces.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace LocalEvaluation {


// forward, to be used in the traits
template< class LocalizableFunctionImp >
class Elliptic;


/**
 *  \brief  Traits for the Elliptic evaluation.
 */
template< class LocalizableFunctionImp >
class EllipticTraits
{
public:
  typedef Elliptic< LocalizableFunctionImp >  derived_type;
  typedef LocalizableFunctionImp              LocalizableFunctionType;
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, LocalizableFunctionImp >::value,
                "LocalizableFunctionImp has to be derived from Stuff::IsLocalizableFunction.");
};


/**
 *  \brief  Computes an elliptic evaluation.
 */
template< class LocalizableFunctionImp >
class Elliptic
  : public LocalEvaluation::Codim0Interface< EllipticTraits< LocalizableFunctionImp >, 2 >
{
public:
  typedef EllipticTraits< LocalizableFunctionImp > Traits;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;

  Elliptic(const LocalizableFunctionType& inducingFunction)
    : inducingFunction_(inducingFunction)
  {}

  template< class EntityType >
  class LocalfunctionTuple
  {
    typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;
  public:
    typedef std::tuple< std::shared_ptr< LocalfunctionType > > Type;
  };

  template< class EntityType >
  typename LocalfunctionTuple< EntityType >::Type localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(inducingFunction_.local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template< class E, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  size_t order(const typename LocalfunctionTuple< E >::Type& localFuncs,
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
    return std::max(int(localFunction.order() + testBase.order() + ansatzBase.order() - 2), 0);
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template< class E, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  static void evaluate(const typename LocalfunctionTuple< E >::Type& localFuncs,
                       const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
                       const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase,
                       const Dune::FieldVector< D, d >& localPoint,
                       Dune::DynamicMatrix< R >& ret)
  {
    const auto localFunction = std::get< 0 >(localFuncs);
    evaluate(*localFunction, testBase, ansatzBase, localPoint, ret);
  }

  template< class E, class D, int d, class R, int rL, int rCL, int rT, int rCT, int rA, int rCA >
  static void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, rL, rCL >& /*localFunction*/,
                       const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& /*testBase*/,
                       const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& /*ansatzBase*/,
                       const Dune::FieldVector< D, d >& /*localPoint*/,
                       Dune::DynamicMatrix< R >& /*ret*/)
  {
    static_assert(Dune::AlwaysFalse< R >::value, "Not implemented for these dimensions!");
  }

  /**
   *  \brief  Computes an elliptic evaluation for a scalar local function and scalar basefunctionsets.
   *  \tparam E EntityType
   *  \tparam D DomainFieldType
   *  \tparam d dimDomain
   *  \tparam R RangeFieldType
   */
  template< class E, class D, int d, class R >
  static void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, 1, 1 >& localFunction,
                       const Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >& testBase,
                       const Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >& ansatzBase,
                       const Dune::FieldVector< D, d >& localPoint,
                       Dune::DynamicMatrix< R >& ret)
  {
    typedef typename Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >::RangeType          RangeType;
    typedef typename Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >::JacobianRangeType  JacobianRangeType;
    // evaluate local function
    const RangeType functionValue = localFunction.evaluate(localPoint);
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
        const R gradientProduct = ansatzGradients[jj][0] * testGradients[ii][0];
        retRow[jj] = functionValue * gradientProduct;
      }
    }
  } // ... evaluate< ..., 1, 1 >(...)

  /**
   *  \brief  Computes an elliptic evaluation for a matrix-valued local function and matrix-valued basefunctionsets.
   *  \tparam E EntityType
   *  \tparam D DomainFieldType
   *  \tparam d dimDomain
   *  \tparam R RangeFieldType
   */
  template< class E, class D, int d, class R>
  static void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, d, d >& localFunction,
                       const Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >& testBase,
                       const Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >& ansatzBase,
                       const Dune::FieldVector< D, d >& localPoint,
                       Dune::DynamicMatrix< R >& ret)
  {
    typedef typename Stuff::LocalfunctionInterface< E, D, d, R, d, d >::RangeType             DiffusionRangeType;
    typedef typename Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >::JacobianRangeType  JacobianRangeType;
    typedef typename Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >::RangeType          RangeType;

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
  } // ... evaluate< ..., d, d >(...)

private:
  const LocalizableFunctionType& inducingFunction_;
}; // class LocalElliptic


} // namespace Evaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_EVALUATION_ELLIPTIC_HH
