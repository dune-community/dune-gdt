// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_RHS_HH
#define DUNE_GDT_OPERATORS_FV_RHS_HH

#include <type_traits>

#include <dune/xt/grid/walker/functors.hh>

#include <dune/xt/la/container/interfaces.hh>

#include <dune/gdt/functionals/base.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/integrands/fv.hh>
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/local/operators/integrals.hh>

namespace Dune {
namespace GDT {


// forward
template <class RhsEvaluationImp>
class AdvectionRhsOperator;


namespace internal {


template <class RhsEvaluationImp>
class AdvectionRhsOperatorTraits
{
public:
  typedef AdvectionRhsOperator<RhsEvaluationImp> derived_type;
  typedef RhsEvaluationImp RhsEvaluationType;
  typedef typename RhsEvaluationImp::DomainFieldType FieldType;
  typedef typename RhsEvaluationImp::PartialURangeType JacobianType;
}; // class AdvectionRhsOperatorTraits


} // namespace internal


template <class GridLayerType, class MatrixType, class DiscreteFunctionType>
class MatrixSolveFunctor : public XT::Grid::Functor::Codim0<GridLayerType>
{
  typedef typename XT::Grid::Functor::Codim0<GridLayerType> BaseType;

public:
  using typename BaseType::EntityType;

  MatrixSolveFunctor(const std::vector<MatrixType>& matrices,
                     const DiscreteFunctionType& rhs,
                     DiscreteFunctionType& solution)
    : matrices_(matrices)
    , rhs_(rhs)
    , solution_(solution)
  {}

  virtual void apply_local(const EntityType& entity)
  {
    // get mapper
    const auto& mapper = rhs_.space().mapper();

    // copy rhs to DynamicVector
    DynamicVector<typename MatrixType::value_type> local_solution(mapper.numDofs(entity), 0.);
    DynamicVector<typename MatrixType::value_type> local_rhs(local_solution.size(), 0.);
    const auto& rhs_vector = rhs_.vector();
    auto& solution_vector = solution_.vector();
    const auto global_indices = mapper.globalIndices(entity);
    for (size_t ii = 0; ii < local_rhs.size(); ++ii)
      local_rhs[ii] = rhs_vector.get_entry(global_indices[ii]);
    // solve
    matrices_[rhs_.space().grid_layer().indexSet().index(entity)].solve(local_solution, local_rhs);
    // write solution
    for (size_t ii = 0; ii < local_rhs.size(); ++ii)
      solution_vector.set_entry(global_indices[ii], local_solution[ii]);
  }

private:
  const std::vector<MatrixType>& matrices_;
  const DiscreteFunctionType& rhs_;
  DiscreteFunctionType& solution_;
};

template <class GridLayerType, class MatrixType, class DiscreteFunctionType>
class MatrixApplyFunctor : public XT::Grid::Functor::Codim0<GridLayerType>
{
  typedef typename XT::Grid::Functor::Codim0<GridLayerType> BaseType;

public:
  using typename BaseType::EntityType;

  MatrixApplyFunctor(const std::vector<MatrixType>& matrices,
                     const DiscreteFunctionType& vector,
                     DiscreteFunctionType& result)
    : matrices_(matrices)
    , vector_(vector)
    , result_(result)
  {}

  virtual void apply_local(const EntityType& entity)
  {
    // get mapper
    const auto& mapper = vector_.space().mapper();

    // copy rhs to DynamicVector
    DynamicVector<typename MatrixType::value_type> local_vector(mapper.numDofs(entity), 0.);
    DynamicVector<typename MatrixType::value_type> local_result(local_vector.size(), 0.);
    const auto& vector_vector = vector_.vector();
    auto& result_vector = result_.vector();
    const auto global_indices = mapper.globalIndices(entity);
    for (size_t ii = 0; ii < local_vector.size(); ++ii)
      local_vector[ii] = vector_vector.get_entry(global_indices[ii]);
    matrices_[vector_.space().grid_layer().indexSet().index(entity)].mv(local_vector, local_result);

    // write solution
    for (size_t ii = 0; ii < local_vector.size(); ++ii)
      result_vector.set_entry(global_indices[ii], local_result[ii]);
  }

private:
  const std::vector<MatrixType>& matrices_;
  const DiscreteFunctionType& vector_;
  DiscreteFunctionType& result_;
};

template <class RhsEvaluationImp>
class AdvectionRhsOperator : public Dune::GDT::OperatorInterface<internal::AdvectionRhsOperatorTraits<RhsEvaluationImp>>
{
  //  static_assert(is_rhs_evaluation<RhsEvaluationImp>::value, "RhsEvaluationImp has to be derived from
  //  RhsInterface!");

public:
  typedef internal::AdvectionRhsOperatorTraits<RhsEvaluationImp> Traits;
  typedef typename Traits::RhsEvaluationType RhsEvaluationType;
  using BaseType = typename Dune::GDT::OperatorInterface<internal::AdvectionRhsOperatorTraits<RhsEvaluationImp>>;

  AdvectionRhsOperator(const RhsEvaluationType& rhs_evaluation)
    : rhs_evaluation_(rhs_evaluation)
  {}

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range, const XT::Common::Parameter& param) const
  {
    std::fill(range.vector().begin(), range.vector().end(), 0);
    LocalVolumeIntegralFunctional<LocalFvRhsIntegrand<RhsEvaluationType, SourceType>,
                                  typename RangeType::SpaceType::BaseFunctionSetType>
        local_functional(rhs_evaluation_, source, param);
    VectorFunctionalBase<typename RangeType::VectorType,
                         typename RangeType::SpaceType,
                         typename RangeType::SpaceType::GridLayerType,
                         typename RangeType::DomainFieldType>
        functional_assembler(range.vector(), range.space());
    functional_assembler.append(local_functional);
    functional_assembler.assemble(true);
  }

  const RhsEvaluationType& evaluation() const
  {
    return rhs_evaluation_;
  }

  // assembles jacobian (jacobian is assumed to be zero initially)
  template <class SourceType, class MatrixTraits>
  void assemble_jacobian(XT::LA::MatrixInterface<MatrixTraits, typename SourceType::RangeFieldType>& jac,
                         const SourceType& source,
                         const XT::Common::Parameter& /*param*/ = {}) const
  {
    typedef typename SourceType::SpaceType SpaceImp;
    typedef typename SpaceImp::BaseFunctionSetType BasisType;
    typedef LocalVolumeIntegralOperator<LocalFvRhsJacobianIntegrand<RhsEvaluationType, SourceType>, BasisType>
        LocalOperatorType;
    LocalOperatorType local_operator(rhs_evaluation_, source);
    SystemAssembler<SpaceImp> assembler(source.space());
    assembler.append(local_operator, jac);
    assembler.assemble(true);
  }

private:
  const RhsEvaluationType& rhs_evaluation_;
}; // class AdvectionRhsOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_RHS_HH
