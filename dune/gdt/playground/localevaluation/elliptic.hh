// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PLAYGROUND_EVALUATION_ELLIPTIC_HH
#define DUNE_GDT_PLAYGROUND_EVALUATION_ELLIPTIC_HH

#include <dune/gdt/localevaluation/elliptic.hh>

namespace Dune {
namespace GDT {
namespace LocalEvaluation {


template< class DiffusionFactorType, class DiffusionTensorType >
class Elliptic
  : public LocalEvaluation::Codim0Interface< internal::EllipticTraits< DiffusionFactorType, DiffusionTensorType >, 2 >
{
public:
  typedef internal::EllipticTraits< DiffusionFactorType, DiffusionTensorType > Traits;

  Elliptic(const DiffusionFactorType& diffusion_factor,
           const DiffusionTensorType& diffusion_tensor)
    : diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
  {}

  template< class EntityType >
  class LocalfunctionTuple
  {
    typedef typename DiffusionFactorType::LocalfunctionType LocalDiffusionFactorFunctionType;
    typedef typename DiffusionTensorType::LocalfunctionType LocalDiffusionTensorFunctionType;
  public:
    typedef std::tuple< std::shared_ptr< LocalDiffusionFactorFunctionType >
                      , std::shared_ptr< LocalDiffusionTensorFunctionType > > Type;
  }; // class LocalfunctionTuple

  template< class EntityType >
  typename LocalfunctionTuple< EntityType >::Type localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_factor_.local_function(entity),
                           diffusion_tensor_.local_function(entity));
  } // ... localFunctions(...)

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template< class E, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  size_t order(const typename LocalfunctionTuple< E >::Type& local_functions_tuple,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase) const
  {
    const auto local_diffusion_factor = std::get< 0 >(local_functions_tuple);
    const auto local_diffusion_tensor = std::get< 1 >(local_functions_tuple);
    return order(*local_diffusion_factor, *local_diffusion_tensor, testBase, ansatzBase);
  } // ... order(...)

  template< class E, class D, int d, class R, int rDF, int rCDF, int rDT, int rCDT, int rT, int rCT, int rA, int rCA >
  size_t order(const Stuff::LocalfunctionInterface< E, D, d, R, rDF, rCDF >& local_diffusion_factor,
               const Stuff::LocalfunctionInterface< E, D, d, R, rDT, rCDT >& local_diffusion_tensor,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase) const
  {
    return local_diffusion_factor.order()
        + local_diffusion_tensor.order()
        + std::max(int(testBase.order() - 1), 0)
        + std::max(int(ansatzBase.order() - 1), 0);
  } // ... order(...)

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template< class E, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  void evaluate(const typename LocalfunctionTuple< E >::Type& local_functions_tuple,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase,
                const Dune::FieldVector< D, d >& localPoint,
                Dune::DynamicMatrix< R >& ret) const
  {
    const auto local_diffusion_factor = std::get< 0 >(local_functions_tuple);
    const auto local_diffusion_tensor = std::get< 1 >(local_functions_tuple);
    evaluate(*local_diffusion_factor, *local_diffusion_tensor, testBase, ansatzBase, localPoint, ret);
  }

  template< class E, class D, int d, class R, int rDF, int rCDF, int rDT, int rCDT, int rT, int rCT, int rA, int rCA >
  void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, rDF, rCDF >& /*local_diffusion_factor*/,
                const Stuff::LocalfunctionInterface< E, D, d, R, rDT, rCDT >& /*local_diffusion_tensor*/,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& /*testBase*/,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& /*ansatzBase*/,
                const Dune::FieldVector< D, d >& /*localPoint*/,
                Dune::DynamicMatrix< R >& /*ret*/) const
  {
    static_assert(Dune::AlwaysFalse< R >::value, "Not implemented for these dimensions!");
  }

  template< class E, class D, int d, class R, int r >
  void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, 1, 1 >& local_diffusion_factor,
                const Stuff::LocalfunctionInterface< E, D, d, R, 1, 1 >& local_diffusion_tensor,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, r, 1 >& testBase,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, r, 1 >& ansatzBase,
                const Dune::FieldVector< D, d >& localPoint,
                Dune::DynamicMatrix< R >& ret)
  {
    typedef typename Stuff::LocalfunctionSetInterface< E, D, d, R, r, 1 >::JacobianRangeType JacobianRangeType;
    // evaluate local functions
    const auto local_diffusion_factor_value = local_diffusion_factor.evaluate(localPoint);
    const auto local_diffusion_tensor_value = local_diffusion_tensor.evaluate(localPoint);
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
        retRow[jj] = local_diffusion_factor_value * local_diffusion_tensor_value
                     * (ansatzGradients[jj][0] * testGradients[ii][0]);
      }
    }
  } // ... evaluate< ..., 1, 1 >(...)

  template< class E, class D, int d, class R >
  void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, 1, 1 >& local_diffusion_factor,
                const Stuff::LocalfunctionInterface< E, D, d, R, 2, 2 >& local_diffusion_tensor,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >& testBase,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >& ansatzBase,
                const Dune::FieldVector< D, d >& localPoint,
                Dune::DynamicMatrix< R >& ret)
  {
    evaluate_matrix_valued_(local_diffusion_factor, local_diffusion_tensor, testBase, ansatzBase, localPoint, ret);
  }

  template< class E, class D, int d, class R >
  void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, 1, 1 >& local_diffusion_factor,
                const Stuff::LocalfunctionInterface< E, D, d, R, 3, 3 >& local_diffusion_tensor,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >& testBase,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >& ansatzBase,
                const Dune::FieldVector< D, d >& localPoint,
                Dune::DynamicMatrix< R >& ret)
  {
    evaluate_matrix_valued_(local_diffusion_factor, local_diffusion_tensor, testBase, ansatzBase, localPoint, ret);
  }

private:
  template< class E, class D, int d, class R >
  void evaluate_matrix_valued_(const Stuff::LocalfunctionInterface< E, D, d, R, 1, 1 >& local_diffusion_factor,
                               const Stuff::LocalfunctionInterface< E, D, d, R, d, d >& local_diffusion_tensor,
                               const Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >& testBase,
                               const Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >& ansatzBase,
                               const Dune::FieldVector< D, d >& localPoint,
                               Dune::DynamicMatrix< R >& ret)
  {
    typedef typename Stuff::LocalfunctionSetInterface< E, D, d, R, 1, 1 >::JacobianRangeType JacobianRangeType;
    // evaluate local functions
    const auto local_diffusion_factor_value = local_diffusion_factor.evaluate(localPoint);
    auto local_diffusion_tensor_value = local_diffusion_tensor.evaluate(localPoint);
    local_diffusion_tensor_value *= local_diffusion_factor_value[0];
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
        local_diffusion_tensor_value.mv(ansatzGradients[jj][0], product);
        retRow[jj] = product * testGradients[ii][0];
      }
    }
  } // ... evaluate_matrix_valued_(...)

  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
}; // class LocalElliptic


} // namespace Evaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_EVALUATION_ELLIPTIC_HH
