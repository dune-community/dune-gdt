// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014, 2016 - 2017)
//   Tobias Leibner  (2014, 2016)

#ifndef DUNE_GDT_ASSEMBLER_WRAPPER_HH
#define DUNE_GDT_ASSEMBLER_WRAPPER_HH

#include <type_traits>

#include <dune/common/deprecated.hh>

#include <dune/xt/la/container/interfaces.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/walker/apply-on.hh>
#include <dune/xt/grid/walker/functors.hh>
#include <dune/xt/grid/walker/wrapper.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/local/assembler.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/constraints.hh>

namespace Dune {
namespace GDT {
namespace internal {

// //////////////////////
// // wrap constraints //
// //////////////////////

template <class TestSpaceType, class AnsatzSpaceType, class GridLayerType, class ConstraintsType>
class ConstraintsWrapper : public XT::Grid::internal::Codim0Object<GridLayerType>
{
  static_assert(AlwaysFalse<ConstraintsType>::value, "Please add a specialization for these Constraints!");
};


// given DirichletConstraints

template <class TestSpaceType, class AnsatzSpaceType, class GridLayerType>
class ConstraintsWrapper<TestSpaceType,
                         AnsatzSpaceType,
                         GridLayerType,
                         DirichletConstraints<XT::Grid::extract_intersection_t<GridLayerType>>>
    : public XT::Grid::internal::Codim0Object<GridLayerType>
{
  static_assert(is_space<TestSpaceType>::value, "TestSpaceType has to be derived from SpaceInterface!");
  static_assert(is_space<AnsatzSpaceType>::value, "AnsatzSpaceType has to be derived from SpaceInterface!");
  typedef XT::Grid::internal::Codim0Object<GridLayerType> BaseType;
  typedef DirichletConstraints<XT::Grid::extract_intersection_t<GridLayerType>> ConstraintsType;

public:
  using typename BaseType::EntityType;

  ConstraintsWrapper(const XT::Common::PerThreadValue<const TestSpaceType>& test_space,
                     const XT::Common::PerThreadValue<const AnsatzSpaceType>& ansatz_space,
                     const XT::Grid::ApplyOn::WhichEntity<GridLayerType>* where,
                     ConstraintsType& constraints)
    : test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , where_(where)
    , constraints_(constraints)
    , thread_local_constraints_(constraints_.boundary_info(), constraints_.size())
  {
  }

  bool apply_on(const GridLayerType& gv, const EntityType& entity) const override final
  {
    return where_->apply_on(gv, entity);
  }

  void apply_local(const EntityType& entity) override final
  {
    test_space_->local_constraints(*ansatz_space_, entity, *thread_local_constraints_);
  }

  void finalize() override final
  {
    DUNE_UNUSED std::lock_guard<std::mutex> mutex_guard(constraints_.mutex_);
    constraints_.dirichlet_DoFs_.insert(thread_local_constraints_->dirichlet_DoFs_.begin(),
                                        thread_local_constraints_->dirichlet_DoFs_.end());
  }

private:
  const XT::Common::PerThreadValue<const TestSpaceType>& test_space_;
  const XT::Common::PerThreadValue<const AnsatzSpaceType>& ansatz_space_;
  const std::unique_ptr<const XT::Grid::ApplyOn::WhichEntity<GridLayerType>> where_;
  ConstraintsType& constraints_;
  XT::Common::PerThreadValue<ConstraintsType> thread_local_constraints_;
}; // class ConstraintsWrapper


// //////////////////////////////////
// // wrap a local volume two-form //
// //////////////////////////////////

// given a local assembler

template <class AssemblerType, class MatrixType>
class DUNE_DEPRECATED_MSG("Use LocalVolumeTwoFormAssemblerFunctor instead or directly append the LocalVolumeTwoForm to "
                          "the SystemAssembler (13.05.2017)!") LocalVolumeTwoFormMatrixAssemblerWrapper
    : public XT::Grid::internal::Codim0Object<typename AssemblerType::GridLayerType>
{
  typedef XT::Grid::internal::Codim0Object<typename AssemblerType::GridLayerType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::AnsatzSpaceType AnsatzSpaceType;
  typedef typename AssemblerType::GridLayerType GridLayerType;
  using typename BaseType::EntityType;
  typedef LocalVolumeTwoFormAssembler<TestSpaceType, MatrixType, AnsatzSpaceType> LocalVolumeTwoFormAssemblerType;

  LocalVolumeTwoFormMatrixAssemblerWrapper(const XT::Common::PerThreadValue<const TestSpaceType>& test_space,
                                           const XT::Common::PerThreadValue<const AnsatzSpaceType>& ansatz_space,
                                           const XT::Grid::ApplyOn::WhichEntity<GridLayerType>* where,
                                           const LocalVolumeTwoFormAssemblerType& local_assembler,
                                           MatrixType& matrix)
    : test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , where_(where)
    , local_assembler_(local_assembler)
    , matrix_(matrix)
  {
  }

  bool apply_on(const GridLayerType& gv, const EntityType& entity) const override final
  {
    return where_->apply_on(gv, entity);
  }

  void apply_local(const EntityType& entity) override final
  {
    apply_local_imp(entity);
  }

private:
  const XT::Common::PerThreadValue<const TestSpaceType>& test_space_;
  const XT::Common::PerThreadValue<const AnsatzSpaceType>& ansatz_space_;
  const std::unique_ptr<const XT::Grid::ApplyOn::WhichEntity<GridLayerType>> where_;
  const LocalVolumeTwoFormAssemblerType& local_assembler_;
  MatrixType& matrix_;
}; // class LocalVolumeTwoFormMatrixAssemblerWrapper


// without a given local assembler

template <class AssemblerType, class MatrixType>
class DUNE_DEPRECATED_MSG("Use LocalVolumeTwoFormAssemblerFunctor instead or directly append the LocalVolumeTwoForm to "
                          "the SystemAssembler (13.05.2017)!") LocalVolumeTwoFormWrapper
    : private XT::Common::ConstStorageProvider<LocalVolumeTwoFormAssembler<typename AssemblerType::TestSpaceType,
                                                                           MatrixType,
                                                                           typename AssemblerType::AnsatzSpaceType>>,
      public LocalVolumeTwoFormMatrixAssemblerWrapper<AssemblerType, MatrixType>
{
  typedef XT::Common::ConstStorageProvider<LocalVolumeTwoFormAssembler<typename AssemblerType::TestSpaceType,
                                                                       MatrixType,
                                                                       typename AssemblerType::AnsatzSpaceType>>
      LocalAssemblerProvider;
  typedef LocalVolumeTwoFormMatrixAssemblerWrapper<AssemblerType, MatrixType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::AnsatzSpaceType AnsatzSpaceType;
  typedef typename AssemblerType::GridLayerType GridLayerType;
  typedef LocalVolumeTwoFormInterface<typename TestSpaceType::BaseFunctionSetType,
                                      typename AnsatzSpaceType::BaseFunctionSetType,
                                      typename MatrixType::ScalarType>
      LocalVolumeTwoFormType;

  LocalVolumeTwoFormWrapper(const XT::Common::PerThreadValue<const TestSpaceType>& test_space,
                            const XT::Common::PerThreadValue<const AnsatzSpaceType>& ansatz_space,
                            const XT::Grid::ApplyOn::WhichEntity<GridLayerType>* where,
                            const LocalVolumeTwoFormType& local_twoform,
                            MatrixType& matrix)
    : LocalAssemblerProvider(local_twoform)
    , BaseType(test_space, ansatz_space, where, LocalAssemblerProvider::access(), matrix)
  {
  }
}; // class LocalVolumeTwoFormWrapper


// ///////////////////////////////// //
// // wrap a local coupling two-form //
// ///////////////////////////////// //

// given a local assembler

template <class AssemblerType, class MatrixType>
class DUNE_DEPRECATED_MSG("Use LocalCouplingTwoFormAssemblerFunctor instead or directly append the "
                          "LocalCouplingTwoForm to the SystemAssembler (13.05.2017)!")
    LocalCouplingTwoFormMatrixAssemblerWrapper
    : public XT::Grid::internal::Codim1Object<typename AssemblerType::GridLayerType>
{
  static_assert(XT::LA::is_matrix<MatrixType>::value, "");

  typedef XT::Grid::internal::Codim1Object<typename AssemblerType::GridLayerType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::AnsatzSpaceType AnsatzSpaceType;
  typedef typename AssemblerType::OuterTestSpaceType OuterTestSpaceType;
  typedef typename AssemblerType::OuterAnsatzSpaceType OuterAnsatzSpaceType;
  typedef typename AssemblerType::GridLayerType GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;
  typedef LocalCouplingTwoFormAssembler<TestSpaceType,
                                        IntersectionType,
                                        MatrixType,
                                        AnsatzSpaceType,
                                        OuterTestSpaceType,
                                        OuterAnsatzSpaceType>
      LocalCouplingTwoFormAssemblerType;

  template <typename TestSpace,
            typename AnsatzSpace,
            typename = typename std::enable_if<(std::is_same<TestSpace, OuterTestSpaceType>::value)
                                               && (std::is_same<AnsatzSpace, OuterAnsatzSpaceType>::value)
                                               && sizeof(TestSpace)
                                               && sizeof(AnsatzSpace)>::type>
  LocalCouplingTwoFormMatrixAssemblerWrapper(const XT::Common::PerThreadValue<const TestSpace>& test_space,
                                             const XT::Common::PerThreadValue<const AnsatzSpace>& ansatz_space,
                                             const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where,
                                             const LocalCouplingTwoFormAssemblerType& local_assembler,
                                             MatrixType& matrix)
    : inner_test_space_(test_space)
    , inner_ansatz_space_(ansatz_space)
    , outer_test_space_(test_space)
    , outer_ansatz_space_(ansatz_space)
    , where_(where)
    , local_assembler_(local_assembler)
    , in_in_matrix_(matrix)
    , out_out_matrix_(matrix)
    , in_out_matrix_(matrix)
    , out_in_matrix_(matrix)
  {
  }

  LocalCouplingTwoFormMatrixAssemblerWrapper(
      const XT::Common::PerThreadValue<const TestSpaceType>& inner_test_space,
      const XT::Common::PerThreadValue<const AnsatzSpaceType>& inner_ansatz_space,
      const XT::Common::PerThreadValue<const OuterTestSpaceType>& outer_test_space,
      const XT::Common::PerThreadValue<const OuterAnsatzSpaceType>& outer_ansatz_space,
      const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where,
      const LocalCouplingTwoFormAssemblerType& local_assembler,
      MatrixType& in_in_matrix,
      MatrixType& out_out_matrix,
      MatrixType& in_out_matrix,
      MatrixType& out_in_matrix)
    : inner_test_space_(inner_test_space)
    , inner_ansatz_space_(inner_ansatz_space)
    , outer_test_space_(outer_test_space)
    , outer_ansatz_space_(outer_ansatz_space)
    , where_(where)
    , local_assembler_(local_assembler)
    , in_in_matrix_(in_in_matrix)
    , out_out_matrix_(out_out_matrix)
    , in_out_matrix_(in_out_matrix)
    , out_in_matrix_(out_in_matrix)
  {
  }

  bool apply_on(const GridLayerType& gv, const IntersectionType& intersection) const override final
  {
    return where_->apply_on(gv, intersection);
  }

  void apply_local(const IntersectionType& intersection,
                   const EntityType& /*inside_entity*/,
                   const EntityType& /*outside_entity*/) override final
  {
    local_assembler_.assemble(*inner_test_space_,
                              *inner_ansatz_space_,
                              *outer_test_space_,
                              *outer_ansatz_space_,
                              intersection,
                              in_in_matrix_,
                              out_out_matrix_,
                              in_out_matrix_,
                              out_in_matrix_);
  } // ... apply_local(...)

private:
  const XT::Common::PerThreadValue<const TestSpaceType>& inner_test_space_;
  const XT::Common::PerThreadValue<const AnsatzSpaceType>& inner_ansatz_space_;
  const XT::Common::PerThreadValue<const OuterTestSpaceType>& outer_test_space_;
  const XT::Common::PerThreadValue<const OuterAnsatzSpaceType>& outer_ansatz_space_;
  const std::unique_ptr<const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>> where_;
  const LocalCouplingTwoFormAssemblerType& local_assembler_;
  MatrixType& in_in_matrix_;
  MatrixType& out_out_matrix_;
  MatrixType& in_out_matrix_;
  MatrixType& out_in_matrix_;
}; // class LocalCouplingTwoFormMatrixAssemblerWrapper

// without a given local assembler

template <class AssemblerType, class MatrixType>
class DUNE_DEPRECATED_MSG("Use LocalCouplingTwoFormAssemblerFunctor instead or directly append the "
                          "LocalCouplingTwoForm to the SystemAssembler (13.05.2017)!") LocalCouplingTwoFormWrapper
    : private XT::Common::
          ConstStorageProvider<LocalCouplingTwoFormAssembler<typename AssemblerType::TestSpaceType,
                                                             typename AssemblerType::IntersectionType,
                                                             MatrixType,
                                                             typename AssemblerType::AnsatzSpaceType,
                                                             typename AssemblerType::OuterTestSpaceType,
                                                             typename AssemblerType::OuterAnsatzSpaceType>>,
      public LocalCouplingTwoFormMatrixAssemblerWrapper<AssemblerType, MatrixType>
{
  typedef XT::Common::ConstStorageProvider<LocalCouplingTwoFormAssembler<typename AssemblerType::TestSpaceType,
                                                                         typename AssemblerType::IntersectionType,
                                                                         MatrixType,
                                                                         typename AssemblerType::AnsatzSpaceType,
                                                                         typename AssemblerType::OuterTestSpaceType,
                                                                         typename AssemblerType::OuterAnsatzSpaceType>>
      LocalAssemblerProvider;
  typedef LocalCouplingTwoFormMatrixAssemblerWrapper<AssemblerType, MatrixType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::AnsatzSpaceType AnsatzSpaceType;
  typedef typename AssemblerType::OuterTestSpaceType OuterTestSpaceType;
  typedef typename AssemblerType::OuterAnsatzSpaceType OuterAnsatzSpaceType;
  typedef typename AssemblerType::GridLayerType GridLayerType;
  typedef LocalCouplingTwoFormInterface<typename TestSpaceType::BaseFunctionSetType,
                                        XT::Grid::extract_intersection_t<GridLayerType>,
                                        typename AnsatzSpaceType::BaseFunctionSetType,
                                        typename OuterTestSpaceType::BaseFunctionSetType,
                                        typename OuterAnsatzSpaceType::BaseFunctionSetType,
                                        typename MatrixType::ScalarType>
      LocalCouplingTwoFormType;

  template <typename TestSpace,
            typename AnsatzSpace,
            typename = typename std::enable_if<(std::is_same<TestSpace, OuterTestSpaceType>::value)
                                               && (std::is_same<AnsatzSpace, OuterAnsatzSpaceType>::value)
                                               && sizeof(TestSpace)
                                               && sizeof(AnsatzSpace)>::type>
  LocalCouplingTwoFormWrapper(const XT::Common::PerThreadValue<const TestSpace>& test_space,
                              const XT::Common::PerThreadValue<const AnsatzSpace>& ansatz_space,
                              const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where,
                              const LocalCouplingTwoFormType& local_twoform,
                              MatrixType& matrix)
    : LocalAssemblerProvider(local_twoform)
    , BaseType(test_space, ansatz_space, where, LocalAssemblerProvider::access(), matrix)
  {
  }

  LocalCouplingTwoFormWrapper(const XT::Common::PerThreadValue<const TestSpaceType>& inner_test_space,
                              const XT::Common::PerThreadValue<const AnsatzSpaceType>& inner_ansatz_space,
                              const XT::Common::PerThreadValue<const OuterTestSpaceType>& outer_test_space,
                              const XT::Common::PerThreadValue<const OuterAnsatzSpaceType>& outer_ansatz_space,
                              const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where,
                              const LocalCouplingTwoFormType& local_twoform,
                              MatrixType& matrix_in_in,
                              MatrixType& matrix_out_out,
                              MatrixType& matrix_in_out,
                              MatrixType& matrix_out_in)
    : LocalAssemblerProvider(local_twoform)
    , BaseType(inner_test_space,
               inner_ansatz_space,
               outer_test_space,
               outer_ansatz_space,
               where,
               LocalAssemblerProvider::access(),
               matrix_in_in,
               matrix_out_out,
               matrix_in_out,
               matrix_out_in)
  {
  }
}; // class LocalCouplingTwoFormWrapper


// ///////////////////////////////// //
// // wrap a local boundary two-form //
// ///////////////////////////////// //

// given a local assembler

template <class AssemblerType, class MatrixType>
class DUNE_DEPRECATED_MSG("Use LocalBoundaryTwoFormAssemblerFunctor instead or directly append the "
                          "LocalBoundaryTwoForm to the SystemAssembler (13.05.2017)!")
    LocalBoundaryTwoFormMatrixAssemblerWrapper
    : public XT::Grid::internal::Codim1Object<typename AssemblerType::GridLayerType>
{
  typedef XT::Grid::internal::Codim1Object<typename AssemblerType::GridLayerType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::AnsatzSpaceType AnsatzSpaceType;
  typedef typename AssemblerType::GridLayerType GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;
  typedef LocalBoundaryTwoFormAssembler<TestSpaceType, IntersectionType, MatrixType, AnsatzSpaceType>
      LocalBoundaryTwoFormAssemblerType;

  LocalBoundaryTwoFormMatrixAssemblerWrapper(const XT::Common::PerThreadValue<const TestSpaceType>& test_space,
                                             const XT::Common::PerThreadValue<const AnsatzSpaceType>& ansatz_space,
                                             const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where,
                                             const LocalBoundaryTwoFormAssemblerType& local_assembler,
                                             MatrixType& matrix)
    : test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , where_(where)
    , local_assembler_(local_assembler)
    , matrix_(matrix)
  {
  }

  bool apply_on(const GridLayerType& gv, const IntersectionType& intersection) const override final
  {
    return where_->apply_on(gv, intersection);
  }

  void apply_local(const IntersectionType& intersection,
                   const EntityType& /*inside_entity*/,
                   const EntityType& /*outside_entity*/) override final
  {
    local_assembler_.assemble(*test_space_, *ansatz_space_, intersection, matrix_);
  } // ... apply_local(...)

private:
  const XT::Common::PerThreadValue<const TestSpaceType>& test_space_;
  const XT::Common::PerThreadValue<const AnsatzSpaceType>& ansatz_space_;
  const std::unique_ptr<const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>> where_;
  const LocalBoundaryTwoFormAssemblerType& local_assembler_;
  MatrixType& matrix_;
}; // class LocalBoundaryTwoFormMatrixAssemblerWrapper

// without a given local assembler

template <class AssemblerType, class MatrixType>
class DUNE_DEPRECATED_MSG("Use LocalBoundaryTwoFormAssemblerFunctor instead or directly append the "
                          "LocalBoundaryTwoForm to the SystemAssembler (13.05.2017)!") LocalBoundaryTwoFormWrapper
    : private XT::Common::ConstStorageProvider<LocalBoundaryTwoFormAssembler<typename AssemblerType::TestSpaceType,
                                                                             typename AssemblerType::IntersectionType,
                                                                             MatrixType,
                                                                             typename AssemblerType::AnsatzSpaceType>>,
      public LocalBoundaryTwoFormMatrixAssemblerWrapper<AssemblerType, MatrixType>
{
  typedef XT::Common::ConstStorageProvider<LocalBoundaryTwoFormAssembler<typename AssemblerType::TestSpaceType,
                                                                         typename AssemblerType::IntersectionType,
                                                                         MatrixType,
                                                                         typename AssemblerType::AnsatzSpaceType>>
      LocalAssemblerProvider;
  typedef LocalBoundaryTwoFormMatrixAssemblerWrapper<AssemblerType, MatrixType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::AnsatzSpaceType AnsatzSpaceType;
  typedef typename AssemblerType::GridLayerType GridLayerType;
  typedef LocalBoundaryTwoFormInterface<typename TestSpaceType::BaseFunctionSetType,
                                        XT::Grid::extract_intersection_t<GridLayerType>,
                                        typename AnsatzSpaceType::BaseFunctionSetType,
                                        typename MatrixType::ScalarType>
      LocalBoundaryTwoFormType;

  LocalBoundaryTwoFormWrapper(const XT::Common::PerThreadValue<const TestSpaceType>& test_space,
                              const XT::Common::PerThreadValue<const AnsatzSpaceType>& ansatz_space,
                              const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where,
                              const LocalBoundaryTwoFormType& local_twoform,
                              MatrixType& matrix)
    : LocalAssemblerProvider(local_twoform)
    , BaseType(test_space, ansatz_space, where, LocalAssemblerProvider::access(), matrix)
  {
  }
}; // class LocalBoundaryTwoFormWrapper


// ////////////////////////////////////
// // wrap a local volume functional //
// ////////////////////////////////////

// given a local assembler

template <class AssemblerType, class VectorType>
class LocalVolumeFunctionalVectorAssemblerWrapper
    : public XT::Grid::internal::Codim0Object<typename AssemblerType::GridLayerType>
{
  typedef XT::Grid::internal::Codim0Object<typename AssemblerType::GridLayerType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::GridLayerType GridLayerType;
  using typename BaseType::EntityType;
  typedef LocalVolumeFunctionalAssembler<TestSpaceType, VectorType> LocalVolumeFunctionalAssemblerType;

  LocalVolumeFunctionalVectorAssemblerWrapper(const XT::Common::PerThreadValue<const TestSpaceType>& space,
                                              const XT::Grid::ApplyOn::WhichEntity<GridLayerType>* where,
                                              const LocalVolumeFunctionalAssemblerType& local_assembler,
                                              VectorType& vector)
    : space_(space)
    , where_(where)
    , local_assembler_(local_assembler)
    , vector_(vector)
  {
  }

  bool apply_on(const GridLayerType& gv, const EntityType& entity) const override final
  {
    return where_->apply_on(gv, entity);
  }

  void apply_local(const EntityType& entity) override final
  {
    local_assembler_.assemble(*space_, entity, vector_);
  }

private:
  const XT::Common::PerThreadValue<const TestSpaceType>& space_;
  const std::unique_ptr<const XT::Grid::ApplyOn::WhichEntity<GridLayerType>> where_;
  const LocalVolumeFunctionalAssemblerType& local_assembler_;
  VectorType& vector_;
}; // class LocalVolumeVectorAssemblerWrapper


// without a given local assembler

template <class AssemblerType, class VectorType>
class LocalVolumeFunctionalWrapper
    : private XT::Common::ConstStorageProvider<LocalVolumeFunctionalAssembler<typename AssemblerType::TestSpaceType,
                                                                              VectorType>>,
      public LocalVolumeFunctionalVectorAssemblerWrapper<AssemblerType, VectorType>
{
  typedef XT::Common::ConstStorageProvider<LocalVolumeFunctionalAssembler<typename AssemblerType::TestSpaceType,
                                                                          VectorType>>
      LocalAssemblerProvider;
  typedef LocalVolumeFunctionalVectorAssemblerWrapper<AssemblerType, VectorType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::GridLayerType GridLayerType;
  typedef LocalVolumeFunctionalInterface<typename TestSpaceType::BaseFunctionSetType, typename VectorType::ScalarType>
      LocalVolumeFunctionalType;

  LocalVolumeFunctionalWrapper(const XT::Common::PerThreadValue<const TestSpaceType>& test_space,
                               const XT::Grid::ApplyOn::WhichEntity<GridLayerType>* where,
                               const LocalVolumeFunctionalType& local_functional,
                               VectorType& vector)
    : LocalAssemblerProvider(local_functional)
    , BaseType(test_space, where, LocalAssemblerProvider::access(), vector)
  {
  }
}; // class LocalVolumeFunctionalWrapper


// /////////////////////////////// //
// // wrap a local face functional //
// /////////////////////////////// //

// given a local assembler

template <class AssemblerType, class VectorType>
class LocalFaceFunctionalVectorAssemblerWrapper
    : public XT::Grid::internal::Codim1Object<typename AssemblerType::GridLayerType>
{
  typedef XT::Grid::internal::Codim1Object<typename AssemblerType::GridLayerType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::GridLayerType GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;
  typedef LocalFaceFunctionalAssembler<TestSpaceType, IntersectionType, VectorType> LocalFaceFunctionalAssemblerType;

  LocalFaceFunctionalVectorAssemblerWrapper(const XT::Common::PerThreadValue<const TestSpaceType>& space,
                                            const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where,
                                            const LocalFaceFunctionalAssemblerType& local_assembler,
                                            VectorType& vector)
    : space_(space)
    , where_(where)
    , local_assembler_(local_assembler)
    , vector_(vector)
  {
  }

  bool apply_on(const GridLayerType& gv, const IntersectionType& intersection) const override final
  {
    return where_->apply_on(gv, intersection);
  }

  void apply_local(const IntersectionType& intersection,
                   const EntityType& /*inside_entity*/,
                   const EntityType& /*outside_entity*/) override final
  {
    local_assembler_.assemble(*space_, intersection, vector_);
  }

private:
  const XT::Common::PerThreadValue<const TestSpaceType>& space_;
  const std::unique_ptr<const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>> where_;
  const LocalFaceFunctionalAssemblerType& local_assembler_;
  VectorType& vector_;
}; // class LocalFaceFunctionalVectorAssemblerWrapper


// wihtout a given local assembler

template <class AssemblerType, class VectorType>
class LocalFaceFunctionalWrapper
    : private XT::Common::ConstStorageProvider<LocalFaceFunctionalAssembler<typename AssemblerType::TestSpaceType,
                                                                            typename AssemblerType::IntersectionType,
                                                                            VectorType>>,
      public LocalFaceFunctionalVectorAssemblerWrapper<AssemblerType, VectorType>
{
  typedef XT::Common::ConstStorageProvider<LocalFaceFunctionalAssembler<typename AssemblerType::TestSpaceType,
                                                                        typename AssemblerType::IntersectionType,
                                                                        VectorType>>
      LocalAssemblerProvider;
  typedef LocalFaceFunctionalVectorAssemblerWrapper<AssemblerType, VectorType> BaseType;

public:
  using typename BaseType::TestSpaceType;
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;
  typedef typename VectorType::ScalarType FieldType;
  typedef LocalFaceFunctionalInterface<typename TestSpaceType::BaseFunctionSetType, IntersectionType, FieldType>
      LocalFaceFunctionalType;

  LocalFaceFunctionalWrapper(const XT::Common::PerThreadValue<const TestSpaceType>& test_space,
                             const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where,
                             const LocalFaceFunctionalType& local_face_functional,
                             VectorType& vector)
    : LocalAssemblerProvider(local_face_functional)
    , BaseType(test_space, where, LocalAssemblerProvider::access(), vector)
  {
  }
}; // class LocalFaceFunctionalWrapper


} // namespace internal
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_WRAPPER_HH
