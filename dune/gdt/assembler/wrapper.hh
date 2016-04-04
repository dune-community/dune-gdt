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

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/constraints.hh>

#include "local.hh"
#include "local/codim0.hh"
#include "local/codim1.hh"
#include "tmp-storage.hh"

namespace Dune {
namespace GDT {
namespace internal {


template <class TestSpaceType, class AnsatzSpaceType, class GridViewType, class ConstraintsType>
class ConstraintsWrapper : public Stuff::Grid::internal::Codim0Object<GridViewType>
{
  static_assert(AlwaysFalse<ConstraintsType>::value, "Please add a specialization for these Constraints!");
};


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


template <class AssemblerType, class LocalTwoFormType, class MatrixType>
class LocalVolumeTwoFormWrapper : public Stuff::Grid::internal::Codim0Object<typename AssemblerType::GridViewType>
{
  typedef Stuff::Grid::internal::Codim0Object<typename AssemblerType::GridViewType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::AnsatzSpaceType AnsatzSpaceType;
  typedef typename AssemblerType::GridViewType GridViewType;
  using typename BaseType::EntityType;

  LocalVolumeTwoFormWrapper(const DS::PerThreadValue<const TestSpaceType>& test_space,
                            const DS::PerThreadValue<const AnsatzSpaceType>& ansatz_space,
                            const Stuff::Grid::ApplyOn::WhichEntity<GridViewType>* where,
                            const LocalTwoFormType& local_twoform, MatrixType& matrix)
    : test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , where_(where)
    , local_twoform_(local_twoform)
    , matrix_(matrix)
    , local_assembler_(local_twoform_)
  {
  }

  virtual ~LocalVolumeTwoFormWrapper() = default;

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
  const LocalTwoFormType& local_twoform_;
  MatrixType& matrix_;
  const LocalVolumeTwoFormAssembler<LocalTwoFormType> local_assembler_;
}; // class LocalVolumeTwoFormWrapper


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


template <class AssemblerType, class LocalVolumeFunctionalAssemblerType, class VectorType>
class LocalVolumeFunctionalVectorAssemblerWrapper
    : public Stuff::Grid::internal::Codim0Object<typename AssemblerType::GridViewType>
{
  typedef Stuff::Grid::internal::Codim0Object<typename AssemblerType::GridViewType> BaseType;
  typedef DSC::TmpVectorsStorage<typename AssemblerType::TestSpaceType::RangeFieldType> TmpVectorsProvider;

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


template <class AssemblerType, class LocalFunctionalType, class VectorType>
class LocalVolumeFunctionalWrapper : public Stuff::Grid::internal::Codim0Object<typename AssemblerType::GridViewType>
{
  typedef Stuff::Grid::internal::Codim0Object<typename AssemblerType::GridViewType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::GridViewType GridViewType;
  using typename BaseType::EntityType;

  LocalVolumeFunctionalWrapper(const DS::PerThreadValue<const TestSpaceType>& test_space,
                               const Stuff::Grid::ApplyOn::WhichEntity<GridViewType>* where,
                               const LocalFunctionalType& local_functional, VectorType& vector)
    : test_space_(test_space)
    , where_(where)
    , local_functional_(local_functional)
    , vector_(vector)
    , local_assembler_(local_functional_)
  {
  }

  virtual ~LocalVolumeFunctionalWrapper() = default;

  virtual bool apply_on(const GridViewType& gv, const EntityType& entity) const override final
  {
    return where_->apply_on(gv, entity);
  }

  virtual void apply_local(const EntityType& entity) override final
  {
    local_assembler_.assemble(*test_space_, entity, vector_);
  }

private:
  const DS::PerThreadValue<const TestSpaceType>& test_space_;
  const std::unique_ptr<const Stuff::Grid::ApplyOn::WhichEntity<GridViewType>> where_;
  const LocalFunctionalType& local_functional_;
  VectorType& vector_;
  const LocalVolumeFunctionalAssembler<LocalFunctionalType> local_assembler_;
}; // class LocalVolumeFunctionalWrapper


template <class AssemblerType, class LocalFunctionalType, class VectorType>
class LocalFaceFunctionalWrapper : public Stuff::Grid::internal::Codim1Object<typename AssemblerType::GridViewType>
{
  typedef Stuff::Grid::internal::Codim1Object<typename AssemblerType::GridViewType> BaseType;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::GridViewType GridViewType;
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;

  LocalFaceFunctionalWrapper(const DS::PerThreadValue<const TestSpaceType>& test_space,
                             const Stuff::Grid::ApplyOn::WhichIntersection<GridViewType>* where,
                             const LocalFunctionalType& local_functional, VectorType& vector)
    : test_space_(test_space)
    , where_(where)
    , local_functional_(local_functional)
    , vector_(vector)
    , local_assembler_(local_functional_)
  {
  }

  virtual ~LocalFaceFunctionalWrapper() = default;

  virtual bool apply_on(const GridViewType& gv, const IntersectionType& intersection) const override final
  {
    return where_->apply_on(gv, intersection);
  }

  virtual void apply_local(const IntersectionType& intersection, const EntityType& /*inside_entity*/,
                           const EntityType& /*outside_entity*/) override final
  {
    local_assembler_.assemble(*test_space_, intersection, vector_);
  }

private:
  const DS::PerThreadValue<const TestSpaceType>& test_space_;
  const std::unique_ptr<const Stuff::Grid::ApplyOn::WhichIntersection<GridViewType>> where_;
  const LocalFunctionalType& local_functional_;
  VectorType& vector_;
  const LocalFaceFunctionalAssembler<LocalFunctionalType> local_assembler_;
}; // class LocalFaceFunctionalWrapper


template <class AssemblerType, class LocalVolumeMatrixAssembler, class MatrixType>
class LocalVolumeMatrixAssemblerWrapper
    : public Stuff::Grid::internal::Codim0Object<typename AssemblerType::GridViewType>,
      DSC::TmpMatricesStorage<typename AssemblerType::TestSpaceType::RangeFieldType>
{
  typedef Stuff::Grid::internal::Codim0Object<typename AssemblerType::GridViewType> BaseType;
  typedef DSC::TmpMatricesStorage<typename AssemblerType::TestSpaceType::RangeFieldType> TmpMatricesProvider;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::AnsatzSpaceType AnsatzSpaceType;
  typedef typename AssemblerType::GridViewType GridViewType;
  using typename BaseType::EntityType;

  LocalVolumeMatrixAssemblerWrapper(const DS::PerThreadValue<const TestSpaceType>& test_space,
                                    const DS::PerThreadValue<const AnsatzSpaceType>& ansatz_space,
                                    const Stuff::Grid::ApplyOn::WhichEntity<GridViewType>* where,
                                    const LocalVolumeMatrixAssembler& localAssembler, MatrixType& matrix)
    : TmpMatricesProvider(localAssembler.numTmpObjectsRequired(), test_space->mapper().maxNumDofs(),
                          ansatz_space->mapper().maxNumDofs())
    , test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , where_(where)
    , localMatrixAssembler_(localAssembler)
    , matrix_(matrix)
  {
  }

  virtual ~LocalVolumeMatrixAssemblerWrapper() = default;

  virtual bool apply_on(const GridViewType& gv, const EntityType& entity) const override final
  {
    return where_->apply_on(gv, entity);
  }

  virtual void apply_local(const EntityType& entity) override final
  {
    localMatrixAssembler_.assembleLocal(
        *test_space_, *ansatz_space_, entity, matrix_, this->matrices(), this->indices());
  }

private:
  const DS::PerThreadValue<const TestSpaceType>& test_space_;
  const DS::PerThreadValue<const AnsatzSpaceType>& ansatz_space_;
  const std::unique_ptr<const Stuff::Grid::ApplyOn::WhichEntity<GridViewType>> where_;
  const LocalVolumeMatrixAssembler& localMatrixAssembler_;
  MatrixType& matrix_;
}; // class LocalVolumeMatrixAssemblerWrapper


template <class AssemblerType, class LocalFaceMatrixAssembler, class MatrixType>
class LocalFaceMatrixAssemblerWrapper
    : public Stuff::Grid::internal::Codim1Object<typename AssemblerType::GridViewType>,
      DSC::TmpMatricesStorage<typename AssemblerType::TestSpaceType::RangeFieldType>
{
  typedef Stuff::Grid::internal::Codim1Object<typename AssemblerType::GridViewType> BaseType;
  typedef DSC::TmpMatricesStorage<typename AssemblerType::TestSpaceType::RangeFieldType> TmpMatricesProvider;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::AnsatzSpaceType AnsatzSpaceType;
  typedef typename AssemblerType::GridViewType GridViewType;
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;

  LocalFaceMatrixAssemblerWrapper(const DS::PerThreadValue<const TestSpaceType>& test_space,
                                  const DS::PerThreadValue<const AnsatzSpaceType>& ansatz_space,
                                  const Stuff::Grid::ApplyOn::WhichIntersection<GridViewType>* where,
                                  const LocalFaceMatrixAssembler& localAssembler, MatrixType& matrix)
    : TmpMatricesProvider(localAssembler.numTmpObjectsRequired(), test_space->mapper().maxNumDofs(),
                          ansatz_space->mapper().maxNumDofs())
    , test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , where_(where)
    , localMatrixAssembler_(localAssembler)
    , matrix_(matrix)
  {
  }

  virtual ~LocalFaceMatrixAssemblerWrapper() = default;

  virtual bool apply_on(const GridViewType& gv, const IntersectionType& intersection) const override final
  {
    return where_->apply_on(gv, intersection);
  }

  virtual void apply_local(const IntersectionType& intersection, const EntityType& /*inside_entity*/,
                           const EntityType& /*outside_entity*/) override final
  {
    localMatrixAssembler_.assembleLocal(
        *test_space_, *ansatz_space_, intersection, matrix_, this->matrices(), this->indices());
  } // ... apply_local(...)

private:
  const DS::PerThreadValue<const TestSpaceType>& test_space_;
  const DS::PerThreadValue<const AnsatzSpaceType>& ansatz_space_;
  const std::unique_ptr<const Stuff::Grid::ApplyOn::WhichIntersection<GridViewType>> where_;
  const LocalFaceMatrixAssembler& localMatrixAssembler_;
  MatrixType& matrix_;
}; // class LocalFaceMatrixAssemblerWrapper


template <class AssemblerType, class LocalVolumeVectorAssembler, class VectorType>
class LocalVolumeVectorAssemblerWrapper
    : public Stuff::Grid::internal::Codim0Object<typename AssemblerType::GridViewType>,
      DSC::TmpVectorsStorage<typename AssemblerType::TestSpaceType::RangeFieldType>
{
  typedef Stuff::Grid::internal::Codim0Object<typename AssemblerType::GridViewType> BaseType;
  typedef DSC::TmpVectorsStorage<typename AssemblerType::TestSpaceType::RangeFieldType> TmpVectorsProvider;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::GridViewType GridViewType;
  using typename BaseType::EntityType;

  LocalVolumeVectorAssemblerWrapper(const DS::PerThreadValue<const TestSpaceType>& space,
                                    const Stuff::Grid::ApplyOn::WhichEntity<GridViewType>* where,
                                    const LocalVolumeVectorAssembler& localAssembler, VectorType& vector)
    : TmpVectorsProvider(localAssembler.numTmpObjectsRequired(), space->mapper().maxNumDofs())
    , space_(space)
    , where_(where)
    , localVectorAssembler_(localAssembler)
    , vector_(vector)
  {
  }

  virtual ~LocalVolumeVectorAssemblerWrapper() = default;

  virtual bool apply_on(const GridViewType& gv, const EntityType& entity) const override final
  {
    return where_->apply_on(gv, entity);
  }

  virtual void apply_local(const EntityType& entity) override final
  {
    localVectorAssembler_.assembleLocal(*space_, entity, vector_, this->vectors(), this->indices());
  }

private:
  const DS::PerThreadValue<const TestSpaceType>& space_;
  const std::unique_ptr<const Stuff::Grid::ApplyOn::WhichEntity<GridViewType>> where_;
  const LocalVolumeVectorAssembler& localVectorAssembler_;
  VectorType& vector_;
}; // class LocalVolumeVectorAssemblerWrapper


template <class AssemblerType, class LocalFaceVectorAssembler, class VectorType>
class LocalFaceVectorAssemblerWrapper
    : public Stuff::Grid::internal::Codim1Object<typename AssemblerType::GridViewType>,
      DSC::TmpVectorsStorage<typename AssemblerType::TestSpaceType::RangeFieldType>
{
  typedef Stuff::Grid::internal::Codim1Object<typename AssemblerType::GridViewType> BaseType;
  typedef DSC::TmpVectorsStorage<typename AssemblerType::TestSpaceType::RangeFieldType> TmpVectorsProvider;

public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::GridViewType GridViewType;
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;

  LocalFaceVectorAssemblerWrapper(const DS::PerThreadValue<const TestSpaceType>& space,
                                  const Stuff::Grid::ApplyOn::WhichIntersection<GridViewType>* where,
                                  const LocalFaceVectorAssembler& localAssembler, VectorType& vector)
    : TmpVectorsProvider(localAssembler.numTmpObjectsRequired(), space->mapper().maxNumDofs())
    , space_(space)
    , where_(where)
    , localVectorAssembler_(localAssembler)
    , vector_(vector)
  {
  }

  virtual ~LocalFaceVectorAssemblerWrapper() = default;

  virtual bool apply_on(const GridViewType& gv, const IntersectionType& intersection) const override final
  {
    return where_->apply_on(gv, intersection);
  }

  virtual void apply_local(const IntersectionType& intersection, const EntityType& /*inside_entity*/,
                           const EntityType& /*outside_entity*/) override final
  {
    localVectorAssembler_.assembleLocal(*space_, intersection, vector_, this->vectors(), this->indices());
  }

private:
  const DS::PerThreadValue<const TestSpaceType>& space_;
  const std::unique_ptr<const Stuff::Grid::ApplyOn::WhichIntersection<GridViewType>> where_;
  const LocalFaceVectorAssembler& localVectorAssembler_;
  VectorType& vector_;
}; // class LocalFaceVectorAssemblerWrapper


} // namespace internal
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_WRAPPER_HH
