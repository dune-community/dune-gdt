#ifndef DUNE_GDT_ASSEMBLER_WRAPPER_HH
#define DUNE_GDT_ASSEMBLER_WRAPPER_HH

#include <type_traits>

#include <dune/common/deprecated.hh>

#include <dune/stuff/common/tmp-storage.hh>
#include <dune/stuff/la/container/interfaces.hh>
#include <dune/stuff/grid/walker.hh>
#include <dune/stuff/grid/walker/apply-on.hh>
#include <dune/stuff/grid/walker/functors.hh>
#include <dune/stuff/grid/walker/wrapper.hh>

#include <dune/gdt/local/assembler.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/constraints.hh>

namespace Dune {
namespace GDT {
namespace internal {

// //////////////////////
// // wrap constraints //
// //////////////////////

template <class TestSpaceType, class AnsatzSpaceType, class GridViewType, class ConstraintsType>
class ConstraintsWrapper : public Stuff::Grid::internal::Codim0Object<GridViewType>
{
  static_assert(AlwaysFalse<ConstraintsType>::value, "Please add a specialization for these Constraints!");
};


// given DirichletConstraints

template <class TestSpaceType, class AnsatzSpaceType, class GridViewType>
class ConstraintsWrapper<TestSpaceType, AnsatzSpaceType, GridViewType,
                         Spaces::DirichletConstraints<typename GridViewType::Intersection>>
    : public Stuff::Grid::internal::Codim0Object<GridViewType>
{
  static_assert(is_space<TestSpaceType>::value, "TestSpaceType has to be derived from SpaceInterface!");
  static_assert(is_space<AnsatzSpaceType>::value, "AnsatzSpaceType has to be derived from SpaceInterface!");
  typedef Stuff::Grid::internal::Codim0Object<GridViewType> BaseType;
  typedef Spaces::DirichletConstraints<typename GridViewType::Intersection> ConstraintsType;

public:
  using typename BaseType::EntityType;

  ConstraintsWrapper(const DS::PerThreadValue<const TestSpaceType>& test_space,
                     const DS::PerThreadValue<const AnsatzSpaceType>& ansatz_space,
                     const Stuff::Grid::ApplyOn::WhichEntity<GridViewType>* where, ConstraintsType& constraints)
    : test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , where_(where)
    , constraints_(constraints)
    , thread_local_constraints_(constraints_.boundary_info(), constraints_.size())
  {
  }

  virtual ~ConstraintsWrapper() = default;

  virtual bool apply_on(const GridViewType& gv, const EntityType& entity) const override final
  {
    return where_->apply_on(gv, entity);
  }

  virtual void apply_local(const EntityType& entity) override final
  {
    test_space_->local_constraints(*ansatz_space_, entity, *thread_local_constraints_);
  }

  virtual void finalize() override final
  {
    std::lock_guard<std::mutex> DUNE_UNUSED(mutex_guard)(constraints_.mutex_);
    constraints_.dirichlet_DoFs_.insert(thread_local_constraints_->dirichlet_DoFs_.begin(),
                                        thread_local_constraints_->dirichlet_DoFs_.end());
  }

private:
  const DS::PerThreadValue<const TestSpaceType>& test_space_;
  const DS::PerThreadValue<const AnsatzSpaceType>& ansatz_space_;
  const std::unique_ptr<const Stuff::Grid::ApplyOn::WhichEntity<GridViewType>> where_;
  ConstraintsType& constraints_;
  DS::PerThreadValue<ConstraintsType> thread_local_constraints_;
}; // class ConstraintsWrapper


// //////////////////////////////////
// // wrap a local volume two-form //
// //////////////////////////////////

// given a local assembler

template <class AssemblerType, class LocalVolumeTwoFormAssemblerType, class MatrixType>
class LocalVolumeTwoFormMatrixAssemblerWrapper
    : public Stuff::Grid::internal::Codim0Object<typename AssemblerType::GridViewType>
{
  typedef Stuff::Grid::internal::Codim0Object<typename AssemblerType::GridViewType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::AnsatzSpaceType AnsatzSpaceType;
  typedef typename AssemblerType::GridViewType GridViewType;
  using typename BaseType::EntityType;

  LocalVolumeTwoFormMatrixAssemblerWrapper(const DS::PerThreadValue<const TestSpaceType>& test_space,
                                           const DS::PerThreadValue<const AnsatzSpaceType>& ansatz_space,
                                           const Stuff::Grid::ApplyOn::WhichEntity<GridViewType>* where,
                                           const LocalVolumeTwoFormAssemblerType& local_assembler, MatrixType& matrix)
    : test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , where_(where)
    , local_assembler_(local_assembler)
    , matrix_(matrix)
  {
  }

  virtual ~LocalVolumeTwoFormMatrixAssemblerWrapper() = default;

  virtual bool apply_on(const GridViewType& gv, const EntityType& entity) const override final
  {
    return where_->apply_on(gv, entity);
  }

  virtual void apply_local(const EntityType& entity) override final
  {
    local_assembler_.assemble(*test_space_, *ansatz_space_, entity, matrix_);
  }

private:
  const DS::PerThreadValue<const TestSpaceType>& test_space_;
  const DS::PerThreadValue<const AnsatzSpaceType>& ansatz_space_;
  const std::unique_ptr<const Stuff::Grid::ApplyOn::WhichEntity<GridViewType>> where_;
  const LocalVolumeTwoFormAssemblerType& local_assembler_;
  MatrixType& matrix_;
}; // class LocalVolumeTwoFormMatrixAssemblerWrapper


// without a given local assembler

template <class AssemblerType, class LocalVolumeTwoFormType, class MatrixType>
class LocalVolumeTwoFormWrapper
    : private DSC::ConstStorageProvider<LocalVolumeTwoFormAssembler<LocalVolumeTwoFormType>>,
      public LocalVolumeTwoFormMatrixAssemblerWrapper<AssemblerType,
                                                      LocalVolumeTwoFormAssembler<LocalVolumeTwoFormType>, MatrixType>
{
  typedef DSC::ConstStorageProvider<LocalVolumeTwoFormAssembler<LocalVolumeTwoFormType>> LocalAssemblerProvider;
  typedef LocalVolumeTwoFormMatrixAssemblerWrapper<AssemblerType, LocalVolumeTwoFormAssembler<LocalVolumeTwoFormType>,
                                                   MatrixType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::AnsatzSpaceType AnsatzSpaceType;
  typedef typename AssemblerType::GridViewType GridViewType;

  LocalVolumeTwoFormWrapper(const DS::PerThreadValue<const TestSpaceType>& test_space,
                            const DS::PerThreadValue<const AnsatzSpaceType>& ansatz_space,
                            const Stuff::Grid::ApplyOn::WhichEntity<GridViewType>* where,
                            const LocalVolumeTwoFormType& local_twoform, MatrixType& matrix)
    : LocalAssemblerProvider(local_twoform)
    , BaseType(test_space, ansatz_space, where, LocalAssemblerProvider::access(), matrix)
  {
  }
}; // class LocalVolumeTwoFormWrapper


// ///////////////////////////////// //
// // wrap a local coupling two-form //
// ///////////////////////////////// //

// given a local assembler

template <class AssemblerType, class LocalCouplingTwoFormAssemblerType, class MatrixType>
class LocalCouplingTwoFormMatrixAssemblerWrapper
    : public Stuff::Grid::internal::Codim1Object<typename AssemblerType::GridViewType>
{
  typedef Stuff::Grid::internal::Codim1Object<typename AssemblerType::GridViewType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::AnsatzSpaceType AnsatzSpaceType;
  typedef typename AssemblerType::GridViewType GridViewType;
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;

  LocalCouplingTwoFormMatrixAssemblerWrapper(const DS::PerThreadValue<const TestSpaceType>& test_space,
                                             const DS::PerThreadValue<const AnsatzSpaceType>& ansatz_space,
                                             const Stuff::Grid::ApplyOn::WhichIntersection<GridViewType>* where,
                                             const LocalCouplingTwoFormAssemblerType& local_assembler,
                                             MatrixType& matrix)
    : test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , where_(where)
    , local_assembler_(local_assembler)
    , matrix_(matrix)
  {
  }

  virtual ~LocalCouplingTwoFormMatrixAssemblerWrapper() = default;

  virtual bool apply_on(const GridViewType& gv, const IntersectionType& intersection) const override final
  {
    return where_->apply_on(gv, intersection);
  }

  virtual void apply_local(const IntersectionType& intersection, const EntityType& /*inside_entity*/,
                           const EntityType& /*outside_entity*/) override final
  {
    local_assembler_.assemble(*test_space_, *ansatz_space_, *test_space_, *ansatz_space_, intersection, matrix_);
  } // ... apply_local(...)

private:
  const DS::PerThreadValue<const TestSpaceType>& test_space_;
  const DS::PerThreadValue<const AnsatzSpaceType>& ansatz_space_;
  const std::unique_ptr<const Stuff::Grid::ApplyOn::WhichIntersection<GridViewType>> where_;
  const LocalCouplingTwoFormAssemblerType& local_assembler_;
  MatrixType& matrix_;
}; // class LocalCouplingTwoFormMatrixAssemblerWrapper

// without a given local assembler

template <class AssemblerType, class LocalCouplingTwoFormType, class MatrixType>
class LocalCouplingTwoFormWrapper
    : private DSC::ConstStorageProvider<LocalCouplingTwoFormAssembler<LocalCouplingTwoFormType>>,
      public LocalCouplingTwoFormMatrixAssemblerWrapper<AssemblerType,
                                                        LocalCouplingTwoFormAssembler<LocalCouplingTwoFormType>,
                                                        MatrixType>
{
  typedef DSC::ConstStorageProvider<LocalCouplingTwoFormAssembler<LocalCouplingTwoFormType>> LocalAssemblerProvider;
  typedef LocalCouplingTwoFormMatrixAssemblerWrapper<AssemblerType,
                                                     LocalCouplingTwoFormAssembler<LocalCouplingTwoFormType>,
                                                     MatrixType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::AnsatzSpaceType AnsatzSpaceType;
  typedef typename AssemblerType::GridViewType GridViewType;

  LocalCouplingTwoFormWrapper(const DS::PerThreadValue<const TestSpaceType>& test_space,
                              const DS::PerThreadValue<const AnsatzSpaceType>& ansatz_space,
                              const Stuff::Grid::ApplyOn::WhichIntersection<GridViewType>* where,
                              const LocalCouplingTwoFormType& local_twoform, MatrixType& matrix)
    : LocalAssemblerProvider(local_twoform)
    , BaseType(test_space, ansatz_space, where, LocalAssemblerProvider::access(), matrix)
  {
  }
}; // class LocalCouplingTwoFormWrapper


// ///////////////////////////////// //
// // wrap a local boundary two-form //
// ///////////////////////////////// //

// given a local assembler

template <class AssemblerType, class LocalBoundaryTwoFormAssemblerType, class MatrixType>
class LocalBoundaryTwoFormMatrixAssemblerWrapper
    : public Stuff::Grid::internal::Codim1Object<typename AssemblerType::GridViewType>
{
  typedef Stuff::Grid::internal::Codim1Object<typename AssemblerType::GridViewType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::AnsatzSpaceType AnsatzSpaceType;
  typedef typename AssemblerType::GridViewType GridViewType;
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;

  LocalBoundaryTwoFormMatrixAssemblerWrapper(const DS::PerThreadValue<const TestSpaceType>& test_space,
                                             const DS::PerThreadValue<const AnsatzSpaceType>& ansatz_space,
                                             const Stuff::Grid::ApplyOn::WhichIntersection<GridViewType>* where,
                                             const LocalBoundaryTwoFormAssemblerType& local_assembler,
                                             MatrixType& matrix)
    : test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , where_(where)
    , local_assembler_(local_assembler)
    , matrix_(matrix)
  {
  }

  virtual ~LocalBoundaryTwoFormMatrixAssemblerWrapper() = default;

  virtual bool apply_on(const GridViewType& gv, const IntersectionType& intersection) const override final
  {
    return where_->apply_on(gv, intersection);
  }

  virtual void apply_local(const IntersectionType& intersection, const EntityType& /*inside_entity*/,
                           const EntityType& /*outside_entity*/) override final
  {
    local_assembler_.assemble(*test_space_, *ansatz_space_, intersection, matrix_);
  } // ... apply_local(...)

private:
  const DS::PerThreadValue<const TestSpaceType>& test_space_;
  const DS::PerThreadValue<const AnsatzSpaceType>& ansatz_space_;
  const std::unique_ptr<const Stuff::Grid::ApplyOn::WhichIntersection<GridViewType>> where_;
  const LocalBoundaryTwoFormAssemblerType& local_assembler_;
  MatrixType& matrix_;
}; // class LocalBoundaryTwoFormMatrixAssemblerWrapper

// without a given local assembler

template <class AssemblerType, class LocalBoundaryTwoFormType, class MatrixType>
class LocalBoundaryTwoFormWrapper
    : private DSC::ConstStorageProvider<LocalBoundaryTwoFormAssembler<LocalBoundaryTwoFormType>>,
      public LocalBoundaryTwoFormMatrixAssemblerWrapper<AssemblerType,
                                                        LocalBoundaryTwoFormAssembler<LocalBoundaryTwoFormType>,
                                                        MatrixType>
{
  typedef DSC::ConstStorageProvider<LocalBoundaryTwoFormAssembler<LocalBoundaryTwoFormType>> LocalAssemblerProvider;
  typedef LocalBoundaryTwoFormMatrixAssemblerWrapper<AssemblerType,
                                                     LocalBoundaryTwoFormAssembler<LocalBoundaryTwoFormType>,
                                                     MatrixType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::AnsatzSpaceType AnsatzSpaceType;
  typedef typename AssemblerType::GridViewType GridViewType;

  LocalBoundaryTwoFormWrapper(const DS::PerThreadValue<const TestSpaceType>& test_space,
                              const DS::PerThreadValue<const AnsatzSpaceType>& ansatz_space,
                              const Stuff::Grid::ApplyOn::WhichIntersection<GridViewType>* where,
                              const LocalBoundaryTwoFormType& local_twoform, MatrixType& matrix)
    : LocalAssemblerProvider(local_twoform)
    , BaseType(test_space, ansatz_space, where, LocalAssemblerProvider::access(), matrix)
  {
  }
}; // class LocalBoundaryTwoFormWrapper


// ////////////////////////////////////
// // wrap a local volume functional //
// ////////////////////////////////////

// given a local assembler

template <class AssemblerType, class LocalVolumeFunctionalAssemblerType, class VectorType>
class LocalVolumeFunctionalVectorAssemblerWrapper
    : public Stuff::Grid::internal::Codim0Object<typename AssemblerType::GridViewType>
{
  typedef Stuff::Grid::internal::Codim0Object<typename AssemblerType::GridViewType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::GridViewType GridViewType;
  using typename BaseType::EntityType;

  LocalVolumeFunctionalVectorAssemblerWrapper(const DS::PerThreadValue<const TestSpaceType>& space,
                                              const Stuff::Grid::ApplyOn::WhichEntity<GridViewType>* where,
                                              const LocalVolumeFunctionalAssemblerType& local_assembler,
                                              VectorType& vector)
    : space_(space)
    , where_(where)
    , local_assembler_(local_assembler)
    , vector_(vector)
  {
  }

  virtual ~LocalVolumeFunctionalVectorAssemblerWrapper() = default;

  virtual bool apply_on(const GridViewType& gv, const EntityType& entity) const override final
  {
    return where_->apply_on(gv, entity);
  }

  virtual void apply_local(const EntityType& entity) override final
  {
    local_assembler_.assemble(*space_, entity, vector_);
  }

private:
  const DS::PerThreadValue<const TestSpaceType>& space_;
  const std::unique_ptr<const Stuff::Grid::ApplyOn::WhichEntity<GridViewType>> where_;
  const LocalVolumeFunctionalAssemblerType& local_assembler_;
  VectorType& vector_;
}; // class LocalVolumeVectorAssemblerWrapper


// without a given local assembler

template <class AssemblerType, class LocalFunctionalType, class VectorType>
class LocalVolumeFunctionalWrapper
    : private DSC::ConstStorageProvider<LocalVolumeFunctionalAssembler<LocalFunctionalType>>,
      public LocalVolumeFunctionalVectorAssemblerWrapper<AssemblerType,
                                                         LocalVolumeFunctionalAssembler<LocalFunctionalType>,
                                                         VectorType>
{
  typedef DSC::ConstStorageProvider<LocalVolumeFunctionalAssembler<LocalFunctionalType>> LocalAssemblerProvider;
  typedef LocalVolumeFunctionalVectorAssemblerWrapper<AssemblerType,
                                                      LocalVolumeFunctionalAssembler<LocalFunctionalType>,
                                                      VectorType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::GridViewType GridViewType;

  LocalVolumeFunctionalWrapper(const DS::PerThreadValue<const TestSpaceType>& test_space,
                               const Stuff::Grid::ApplyOn::WhichEntity<GridViewType>* where,
                               const LocalFunctionalType& local_functional, VectorType& vector)
    : LocalAssemblerProvider(local_functional)
    , BaseType(test_space, where, LocalAssemblerProvider::access(), vector)
  {
  }
}; // class LocalVolumeFunctionalWrapper


// /////////////////////////////// //
// // wrap a local face functional //
// /////////////////////////////// //

// given a local assembler

template <class AssemblerType, class LocalFaceFunctionalAssemblerType, class VectorType>
class LocalFaceFunctionalVectorAssemblerWrapper
    : public Stuff::Grid::internal::Codim1Object<typename AssemblerType::GridViewType>
{
  typedef Stuff::Grid::internal::Codim1Object<typename AssemblerType::GridViewType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::GridViewType GridViewType;
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;

  LocalFaceFunctionalVectorAssemblerWrapper(const DS::PerThreadValue<const TestSpaceType>& space,
                                            const Stuff::Grid::ApplyOn::WhichIntersection<GridViewType>* where,
                                            const LocalFaceFunctionalAssemblerType& local_assembler, VectorType& vector)
    : space_(space)
    , where_(where)
    , local_assembler_(local_assembler)
    , vector_(vector)
  {
  }

  virtual ~LocalFaceFunctionalVectorAssemblerWrapper() = default;

  virtual bool apply_on(const GridViewType& gv, const IntersectionType& intersection) const override final
  {
    return where_->apply_on(gv, intersection);
  }

  virtual void apply_local(const IntersectionType& intersection, const EntityType& /*inside_entity*/,
                           const EntityType& /*outside_entity*/) override final
  {
    local_assembler_.assemble(*space_, intersection, vector_);
  }

private:
  const DS::PerThreadValue<const TestSpaceType>& space_;
  const std::unique_ptr<const Stuff::Grid::ApplyOn::WhichIntersection<GridViewType>> where_;
  const LocalFaceFunctionalAssemblerType& local_assembler_;
  VectorType& vector_;
}; // class LocalFaceFunctionalVectorAssemblerWrapper


// wihtout a given local assembler

template <class AssemblerType, class LocalFaceFunctionalType, class VectorType>
class LocalFaceFunctionalWrapper
    : private DSC::ConstStorageProvider<LocalFaceFunctionalAssembler<LocalFaceFunctionalType>>,
      public LocalFaceFunctionalVectorAssemblerWrapper<AssemblerType,
                                                       LocalFaceFunctionalAssembler<LocalFaceFunctionalType>,
                                                       VectorType>
{
  typedef DSC::ConstStorageProvider<LocalFaceFunctionalAssembler<LocalFaceFunctionalType>> LocalAssemblerProvider;
  typedef LocalFaceFunctionalVectorAssemblerWrapper<AssemblerType,
                                                    LocalFaceFunctionalAssembler<LocalFaceFunctionalType>,
                                                    VectorType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::GridViewType GridViewType;

  LocalFaceFunctionalWrapper(const DS::PerThreadValue<const TestSpaceType>& test_space,
                             const Stuff::Grid::ApplyOn::WhichIntersection<GridViewType>* where,
                             const LocalFaceFunctionalType& local_functional, VectorType& vector)
    : LocalAssemblerProvider(local_functional)
    , BaseType(test_space, where, LocalAssemblerProvider::access(), vector)
  {
  }
}; // class LocalFaceFunctionalWrapper


} // namespace internal
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_WRAPPER_HH
