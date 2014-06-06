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

#include <dune/common/dynmatrix.hh>
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
  typedef typename LocalizableFunctionType::EntityType                                         EntityType;
  typedef typename LocalizableFunctionType::DomainFieldType                                    DomainFieldType;
  typedef std::tuple< std::shared_ptr< typename LocalizableFunctionType::LocalfunctionType > > LocalFunctionTupleType;
  static const unsigned int dimDomain = LocalizableFunctionType::dimDomain;
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
  typedef typename Traits::EntityType             EntityType;
  typedef typename Traits::DomainFieldType        DomainFieldType;
  typedef typename Traits::LocalFunctionTupleType LocalFunctionTupleType;
  static const unsigned int dimDomain = Traits::dimDomain;

  Elliptic(const LocalizableFunctionType& inducing_function)
    : inducing_function_(inducing_function)
  {}

  class DUNE_DEPRECATED_MSG("Please use LocalFunctionTupleType!") LocalfunctionTuple
  {
    typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;
  public:
    typedef std::tuple< std::shared_ptr< LocalfunctionType > > Type;
  };


  LocalFunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(inducing_function_.local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template< class R, int rT, int rCT, int rA, int rCA >
  size_t order(const typename LocalfunctionTuple::Type& local_funcs,
               const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, rT, rCT >& test_base,
               const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatz_base) const
  {
    const auto local_function = std::get< 0 >(local_funcs);
    return order(*local_function, test_base, ansatz_base);
  }

  /**
   *  \todo add copydoc
   *  \return local_function.order() + (test_base.order() - 1) + (ansatz_base.order() - 1)
   */
  template< class R, int rL, int rCL, int rT, int rCT, int rA, int rCA >
  size_t order(const Stuff::LocalfunctionInterface< EntityType, DomainFieldType, dimDomain, R, rL, rCL >& local_function,
               const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, rT, rCT >& test_base,
               const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatz_base) const
  {
    return local_function.order()
        + std::max(ssize_t(test_base.order()) - 1, ssize_t(0))
        + std::max(ssize_t(ansatz_base.order()) - 1, ssize_t(0));
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template< class R, int rT, int rCT, int rA, int rCA >
  void evaluate(const typename LocalfunctionTuple::Type& local_funcs,
                const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, rT, rCT >& test_base,
                const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatz_base,
                const Dune::FieldVector< DomainFieldType, dimDomain >& local_point,
                Dune::DynamicMatrix< R >& ret) const
  {
    const auto local_function = std::get< 0 >(local_funcs);
    evaluate(*local_function, test_base, ansatz_base, local_point, ret);
  }

  template< class R, int rL, int rCL, int rT, int rCT, int rA, int rCA >
  void evaluate(const Stuff::LocalfunctionInterface< EntityType, DomainFieldType, dimDomain, R, rL, rCL >& /*local_function*/,
                const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, rT, rCT >& /*test_base*/,
                const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, rA, rCA >& /*ansatz_base*/,
                const Dune::FieldVector< DomainFieldType, dimDomain >& /*local_point*/,
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
  template< class R, int r >
  void evaluate(const Stuff::LocalfunctionInterface< EntityType, DomainFieldType, dimDomain, R, 1, 1 >& local_function,
                const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, r, 1 >& test_base,
                const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, r, 1 >& ansatz_base,
                const Dune::FieldVector< DomainFieldType, dimDomain >& local_point,
                Dune::DynamicMatrix< R >& ret) const
  {
    typedef typename Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, r, 1 >::JacobianRangeType
        JacobianRangeType;
    // evaluate local function
    const auto functionValue = local_function.evaluate(local_point);
    // evaluate test gradient
    const size_t rows = test_base.size();
    std::vector< JacobianRangeType > testGradients(rows, JacobianRangeType(0));
    test_base.jacobian(local_point, testGradients);
    // evaluate ansatz gradient
    const size_t cols = ansatz_base.size();
    std::vector< JacobianRangeType > ansatzGradients(cols, JacobianRangeType(0));
    ansatz_base.jacobian(local_point, ansatzGradients);
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
  template< class R >
  void evaluate(const Stuff::LocalfunctionInterface< EntityType, DomainFieldType, dimDomain, R, 2, 2 >& local_function,
                const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, 1, 1 >& test_base,
                const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, 1, 1 >& ansatz_base,
                const Dune::FieldVector< DomainFieldType, dimDomain >& local_point,
                Dune::DynamicMatrix< R >& ret) const
  {
    evaluate_matrix_valued_(local_function, test_base, ansatz_base, local_point, ret);
  }

  /**
   *  \brief  Computes an elliptic evaluation for a 3x3 matrix-valued local function and matrix-valued basefunctionsets.
   *  \tparam E EntityType
   *  \tparam D DomainFieldType
   *  \tparam d dimDomain
   *  \tparam R RangeFieldType
   *  \note   Unfortunately we need this explicit specialization, otherwise the compiler will complain for 1d grids.
   */
  template< class R >
  void evaluate(const Stuff::LocalfunctionInterface< EntityType, DomainFieldType, dimDomain, R, 3, 3 >& local_function,
                const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, 1, 1 >& test_base,
                const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, 1, 1 >& ansatz_base,
                const Dune::FieldVector< DomainFieldType, dimDomain >& local_point,
                Dune::DynamicMatrix< R >& ret) const
  {
    evaluate_matrix_valued_(local_function, test_base, ansatz_base, local_point, ret);
  }

private:
  template< class R >
  void evaluate_matrix_valued_(const Stuff::LocalfunctionInterface< EntityType, DomainFieldType, dimDomain, R, d, d >& local_function,
                               const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, 1, 1 >& test_base,
                               const Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, 1, 1 >& ansatz_base,
                               const Dune::FieldVector< DomainFieldType, dimDomain >& local_point,
                               Dune::DynamicMatrix< R >& ret) const
  {
    typedef typename Stuff::LocalfunctionInterface< EntityType, DomainFieldType, dimDomain, R, dimDomain, dimDomain >
        ::RangeType DiffusionRangeType;
    typedef typename Stuff::LocalfunctionSetInterface< EntityType, DomainFieldType, dimDomain, R, 1, 1 >
        ::JacobianRangeType JacobianRangeType;
    // evaluate local function
    const DiffusionRangeType functionValue = local_function.evaluate(local_point);
    // evaluate test gradient
    const size_t rows = test_base.size();
    std::vector< JacobianRangeType > testGradients(rows, JacobianRangeType(0));
    test_base.jacobian(local_point, testGradients);
    // evaluate ansatz gradient
    const size_t cols = ansatz_base.size();
    std::vector< JacobianRangeType > ansatzGradients(cols, JacobianRangeType(0));
    ansatz_base.jacobian(local_point, ansatzGradients);
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

  const LocalizableFunctionType& inducing_function_;
}; // class LocalElliptic


} // namespace Evaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_EVALUATION_ELLIPTIC_HH
