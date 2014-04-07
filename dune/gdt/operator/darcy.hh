// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATOR_DARCY_HH
#define DUNE_GDT_OPERATOR_DARCY_HH

#include <limits>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/gdt/space/continuouslagrange/fem.hh>
#include <dune/gdt/discretefunction/default.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace Operator {


// forward, to be used in the traits
template <class GridViewImp, class FunctionImp>
class DarcyReconstruction;


template <class GridViewImp, class FunctionImp>
class DarcyReconstructionTraits
{
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, FunctionImp>::value,
                "FunctionImp has to be derived from Stuff::IsLocalizableFunction!");
  static_assert(std::is_same<typename GridViewImp::ctype, typename FunctionImp::DomainFieldType>::value,
                "Types do not match!");
  static_assert(GridViewImp::dimension == FunctionImp::dimDomain, "Dimensions do not match!");

public:
  typedef DarcyReconstruction<GridViewImp, FunctionImp> derived_type;
  typedef GridViewImp GridViewType;
  typedef typename FunctionImp::RangeFieldType FieldType;
}; // class DarcyReconstructionTraits


/**
  * \note Only works for scalar valued function atm.
  **/
template <class GridViewImp, class FunctionImp>
class DarcyReconstruction : public OperatorInterface<DarcyReconstructionTraits<GridViewImp, FunctionImp>>
{
public:
  typedef DarcyReconstructionTraits<GridViewImp, FunctionImp> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::FieldType FieldType;

  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;

  DarcyReconstruction(const GridViewType& grid_view, const FunctionImp& function)
    : grid_view_(grid_view)
    , function_(function)
  {
  }

  template <class E, class D, int d, class R, int r, int rC, class T, class V>
  void apply(const Stuff::LocalizableFunctionInterface<E, D, d, R, r, rC>& /*source*/,
             DiscreteFunction<SpaceInterface<T>, V>& /*range*/) const
  {
    static_assert((Dune::AlwaysFalse<E>::value), "Not implemented for this combination of source and range!");
  }

  /**
   * \brief Does an L2 projection of '- function * \gradient source' onto range.
   */
  template <class GP, int p, class V>
  void apply(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType, 1, 1>& source,
             DiscreteFunction<ContinuousLagrangeSpace::FemWrapper<GP, p, FieldType, dimDomain, 1>, V>& range) const
  {
#if HAVE_EIGEN
    typedef Stuff::LA::EigenRowMajorSparseMatrix<FieldType> MatrixType;
    typedef Stuff::LA::EigenDenseVector<FieldType> VectorType;
#elif HAVE_DUNE_ISTL
    typedef Stuff::LA::IstlRowMajorSparseMatrix<FieldType> MatrixType;
    typedef Stuff::LA::IstlDenseVector<FieldType> VectorType;
#else
    typedef Stuff::LA::CommonDenseMatrix<FieldType> MatrixType;
    typedef Stuff::LA::CommonDenseVector<FieldType> VectorType;
#endif
    MatrixType lhs(
        range.space().mapper().size(), range.space().mapper().size(), range.space().compute_volume_pattern());
    VectorType rhs(range.space().mapper().size());

    // walk the grid
    const auto entity_it_end = grid_view_.template end<0>();
    for (auto entity_it = grid_view_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity        = *entity_it;
      const auto local_function = function_.local_function(entity);
      const auto local_source   = source.local_function(entity);
      const auto basis          = range.space().base_function_set(entity);
      // do a volume quadrature
      const size_t integrand_order =
          std::max(local_function->order() + size_t(local_source->order() - 1), basis.order()) + basis.order();
      assert(integrand_order < std::numeric_limits<int>::max());
      const auto& quadrature       = QuadratureRules<DomainFieldType, dimDomain>::rule(entity.type(), int(integrand_order));
      const auto quadrature_it_end = quadrature.end();
      for (auto quadrature_it = quadrature.begin(); quadrature_it != quadrature_it_end; ++quadrature_it) {
        const auto xx                             = quadrature_it->position();
        const FieldType quadrature_weight         = quadrature_it->weight();
        const DomainFieldType integration_element = entity.geometry().integrationElement(xx);
        const auto function_value                 = local_function->evaluate(xx);
        const auto source_gradient                = local_source->jacobian(xx);
        const auto basis_value = basis.evaluate(xx);
        for (size_t ii = 0; ii < basis.size(); ++ii) {
          const size_t global_ii = range.space().mapper().mapToGlobal(entity, ii);
          rhs.add_to_entry(global_ii,
                           integration_element * quadrature_weight * -1.0 * function_value
                               * (source_gradient[0] * basis_value[ii]));
          for (size_t jj = 0; jj < basis.size(); ++jj) {
            const size_t global_jj = range.space().mapper().mapToGlobal(entity, jj);
            lhs.add_to_entry(
                global_ii, global_jj, integration_element * quadrature_weight * (basis_value[ii] * basis_value[jj]));
          }
        }
      } // do a volume quadrature
    } // walk the grid

    // solve
    Stuff::LA::Solver<MatrixType>(lhs).apply(rhs, range.vector());
  } // ... apply(...)

private:
  const GridViewType& grid_view_;
  const FunctionImp& function_;
}; // class DarcyReconstruction


} // namespace Operator
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATOR_DARCY_HH
