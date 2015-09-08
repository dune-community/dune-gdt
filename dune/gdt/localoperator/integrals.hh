// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LOCALOPERATOR_INTEGRALS_HH
#define DUNE_GDT_LOCALOPERATOR_INTEGRALS_HH

#include <type_traits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/densematrix.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/common/matrix.hh>

#include "../localevaluation/interface.hh"
#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forwards, to be used in the traits
template< class BinaryEvaluationType >
class LocalVolumeIntegralOperator;

namespace internal {


template< class BinaryEvaluationType >
class LocalVolumeIntegralOperatorTraits
{
  static_assert(std::is_base_of< LocalEvaluation::Codim0Interface< typename BinaryEvaluationType::Traits, 2 >,
                                 BinaryEvaluationType >::value,
                "BinaryEvaluationType has to be derived from LocalEvaluation::Codim0Interface< ..., 2 >!");
public:
  typedef LocalVolumeIntegralOperator< BinaryEvaluationType > derived_type;
};


} // namespace internal


template< class BinaryEvaluationType >
class LocalVolumeIntegralOperator
  : public LocalVolumeTwoFormInterface< internal::LocalVolumeIntegralOperatorTraits< BinaryEvaluationType > >
{
  typedef LocalVolumeIntegralOperator< BinaryEvaluationType >                 ThisType;
public:
  typedef internal::LocalVolumeIntegralOperatorTraits< BinaryEvaluationType > Traits;

  template< class... Args >
  explicit LocalVolumeIntegralOperator(Args&& ...args)
    : integrand_(std::forward< Args >(args)...)
    , over_integrate_(0)
  {}

  template< class... Args >
  explicit LocalVolumeIntegralOperator(const int over_integrate, Args&& ...args)
    : integrand_(std::forward< Args >(args)...)
    , over_integrate_(boost::numeric_cast< size_t >(over_integrate))
  {}

  template< class... Args >
  explicit LocalVolumeIntegralOperator(const size_t over_integrate, Args&& ...args)
    : integrand_(std::forward< Args >(args)...)
    , over_integrate_(over_integrate)
  {}

  LocalVolumeIntegralOperator(const ThisType& other) = default;
  LocalVolumeIntegralOperator(ThisType&& source) = default;

  template< class E, class D, size_t d, class R, size_t rT, size_t rCT, size_t rA, size_t rCA >
  void apply2(const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& test_base,
              const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatz_base,
              Dune::DynamicMatrix< R >& ret) const
  {
    const auto& entity = ansatz_base.entity();
    const auto local_functions = integrand_.localFunctions(entity);
    // create quadrature
    const size_t integrand_order = integrand_.order(local_functions, ansatz_base, test_base) + over_integrate_;
    const auto& quadrature = QuadratureRules< D, d >::rule(entity.type(),
                                                           boost::numeric_cast< int >(integrand_order));
    // prepare storage
    const size_t rows = test_base.size();
    const size_t cols = ansatz_base.size();
    ret *= 0.0;
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    Dune::DynamicMatrix< R > evaluation_result(rows, cols, 0.);
    // loop over all quadrature points
    for (const auto& quadrature_point : quadrature) {
      const auto xx = quadrature_point.position();
      // integration factors
      const auto integration_factor = entity.geometry().integrationElement(xx);
      const auto quadrature_weight = quadrature_point.weight();
      // evaluate the integrand
      integrand_.evaluate(local_functions, test_base, ansatz_base, xx, evaluation_result);
      // compute integral
      for (size_t ii = 0; ii < rows; ++ii) {
        const auto& evaluation_result_row = evaluation_result[ii];
        auto& ret_row = ret[ii];
        for (size_t jj = 0; jj < cols; ++jj)
          ret_row[jj] += evaluation_result_row[jj] * integration_factor * quadrature_weight;
      } // compute integral
    } // loop over all quadrature points
  } // ... apply2(...)

private:
  const BinaryEvaluationType integrand_;
  const size_t over_integrate_;
}; // class LocalVolumeIntegralOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALOPERATOR_INTEGRALS_HH
