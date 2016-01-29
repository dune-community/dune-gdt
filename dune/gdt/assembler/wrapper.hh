#ifndef DUNE_GDT_ASSEMBLER_WRAPPER_HH
#define DUNE_GDT_ASSEMBLER_WRAPPER_HH

#include <type_traits>

#include <dune/stuff/common/tmp-storage.hh>
#include <dune/stuff/la/container/interfaces.hh>
#include <dune/stuff/grid/walker.hh>
#include <dune/stuff/grid/walker/apply-on.hh>
#include <dune/stuff/grid/walker/functors.hh>
#include <dune/stuff/grid/walker/wrapper.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/constraints.hh>

#include "local/codim0.hh"
#include "local/codim1.hh"
#include "tmp-storage.hh"

namespace Dune {
namespace GDT {
namespace internal {


template< class TestSpaceType, class AnsatzSpaceType, class GridViewType, class ConstraintsType >
class ConstraintsWrapper
  : public Stuff::Grid::internal::Codim0Object< GridViewType >
{
  static_assert(is_space< TestSpaceType >::value,   "TestSpaceType has to be derived from SpaceInterface!");
  static_assert(is_space< AnsatzSpaceType >::value, "AnsatzSpaceType has to be derived from SpaceInterface!");
  static_assert(std::is_base_of< Spaces::ConstraintsInterface< typename ConstraintsType::Traits >,
                                 ConstraintsType >::value, "");
  typedef Stuff::Grid::internal::Codim0Object< GridViewType > BaseType;
public:
  typedef typename BaseType::EntityType EntityType;

  ConstraintsWrapper(const DS::PerThreadValue< const TestSpaceType >& test_space,
                     const DS::PerThreadValue< const AnsatzSpaceType >& ansatz_space,
                     const Stuff::Grid::ApplyOn::WhichEntity< GridViewType >* where,
                     ConstraintsType& constraints)
    : test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , where_(where)
    , constraints_(constraints)
  {}

  virtual ~ConstraintsWrapper() {}

  virtual bool apply_on(const GridViewType& gv, const EntityType& entity) const override final
  {
    return where_->apply_on(gv, entity);
  }

  virtual void apply_local(const EntityType& entity) override final
  {
    test_space_->local_constraints(*ansatz_space_, entity, constraints_);
  }

private:
  const DS::PerThreadValue< const TestSpaceType >& test_space_;
  const DS::PerThreadValue< const AnsatzSpaceType >& ansatz_space_;
  const std::unique_ptr< const Stuff::Grid::ApplyOn::WhichEntity< GridViewType > > where_;
  ConstraintsType& constraints_;
}; // class ConstraintsWrapper


template< class AssemblerType, class LocalVolumeMatrixAssembler, class MatrixType >
class LocalVolumeMatrixAssemblerWrapper
  : public Stuff::Grid::internal::Codim0Object<typename AssemblerType::GridViewType>
  , DSC::TmpMatricesStorage< typename AssemblerType::TestSpaceType::RangeFieldType >
{
  typedef DSC::TmpMatricesStorage< typename AssemblerType::TestSpaceType::RangeFieldType > TmpMatricesProvider;
public:
  typedef typename AssemblerType::TestSpaceType   TestSpaceType;
  typedef typename AssemblerType::AnsatzSpaceType AnsatzSpaceType;
  typedef typename AssemblerType::GridViewType    GridViewType;
  typedef typename AssemblerType::EntityType      EntityType;

  LocalVolumeMatrixAssemblerWrapper(const DS::PerThreadValue< const TestSpaceType >& test_space,
                                    const DS::PerThreadValue< const AnsatzSpaceType >& ansatz_space,
                                    const Stuff::Grid::ApplyOn::WhichEntity< GridViewType >* where,
                                    const LocalVolumeMatrixAssembler& localAssembler,
                                    MatrixType& matrix)
    : TmpMatricesProvider(localAssembler.numTmpObjectsRequired(),
                          test_space->mapper().maxNumDofs(),
                          ansatz_space->mapper().maxNumDofs())
    , test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , where_(where)
    , localMatrixAssembler_(localAssembler)
    , matrix_(matrix)
  {}

  virtual ~LocalVolumeMatrixAssemblerWrapper() {}

  virtual bool apply_on(const GridViewType& gv, const EntityType& entity) const override final
  {
    return where_->apply_on(gv, entity);
  }

  virtual void apply_local(const EntityType& entity) override final
  {
    localMatrixAssembler_.assembleLocal(*test_space_, *ansatz_space_, entity, matrix_, this->matrices(), this->indices());
  }

private:
  const DS::PerThreadValue< const TestSpaceType >& test_space_;
  const DS::PerThreadValue< const AnsatzSpaceType >& ansatz_space_;
  const std::unique_ptr< const Stuff::Grid::ApplyOn::WhichEntity< GridViewType > > where_;
  const LocalVolumeMatrixAssembler& localMatrixAssembler_;
  MatrixType& matrix_;
}; // class LocalVolumeMatrixAssemblerWrapper


template< class AssemblerType, class LocalFaceMatrixAssembler, class MatrixType >
class LocalFaceMatrixAssemblerWrapper
  : public Stuff::Grid::internal::Codim1Object< typename AssemblerType::GridViewType >
  , DSC::TmpMatricesStorage< typename AssemblerType::TestSpaceType::RangeFieldType >
{
  typedef DSC::TmpMatricesStorage< typename AssemblerType::TestSpaceType::RangeFieldType > TmpMatricesProvider;
public:
  typedef typename AssemblerType::TestSpaceType                                            TestSpaceType;
  typedef typename AssemblerType::AnsatzSpaceType                                          AnsatzSpaceType;
  typedef typename AssemblerType::GridViewType                                             GridViewType;
  typedef typename AssemblerType::EntityType                                               EntityType;
  typedef typename Stuff::Grid::internal::Codim1Object< GridViewType >::IntersectionType   IntersectionType;

  LocalFaceMatrixAssemblerWrapper(const DS::PerThreadValue< const TestSpaceType >& test_space,
                                  const DS::PerThreadValue< const AnsatzSpaceType >& ansatz_space,
                                  const Stuff::Grid::ApplyOn::WhichIntersection< GridViewType >* where,
                                  const LocalFaceMatrixAssembler& localAssembler,
                                  MatrixType& matrix)
    : TmpMatricesProvider(localAssembler.numTmpObjectsRequired(),
                          test_space->mapper().maxNumDofs(),
                          ansatz_space->mapper().maxNumDofs())
    , test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , where_(where)
    , localMatrixAssembler_(localAssembler)
    , matrix_(matrix)
  {}

  virtual ~LocalFaceMatrixAssemblerWrapper() {}

  virtual bool apply_on(const GridViewType& gv, const IntersectionType& intersection) const override final
  {
    return where_->apply_on(gv, intersection);
  }

  virtual void apply_local(const IntersectionType& intersection,
                           const EntityType& /*inside_entity*/,
                           const EntityType& /*outside_entity*/) override final
  {
    localMatrixAssembler_.assembleLocal(*test_space_, *ansatz_space_,
                                        intersection,
                                        matrix_,
                                        this->matrices(), this->indices());
  } // ... apply_local(...)

private:
  const DS::PerThreadValue< const TestSpaceType >& test_space_;
  const DS::PerThreadValue< const AnsatzSpaceType >& ansatz_space_;
  const std::unique_ptr< const Stuff::Grid::ApplyOn::WhichIntersection< GridViewType > > where_;
  const LocalFaceMatrixAssembler& localMatrixAssembler_;
  MatrixType& matrix_;
}; // class LocalFaceMatrixAssemblerWrapper


template< class AssemblerType, class LocalVolumeVectorAssembler, class VectorType >
class LocalVolumeVectorAssemblerWrapper
  : public Stuff::Grid::internal::Codim0Object< typename AssemblerType::GridViewType >
  , DSC::TmpVectorsStorage< typename AssemblerType::TestSpaceType::RangeFieldType >
{
  typedef DSC::TmpVectorsStorage< typename AssemblerType::TestSpaceType::RangeFieldType > TmpVectorsProvider;
public:
  typedef typename AssemblerType::TestSpaceType TestSpaceType;
  typedef typename AssemblerType::GridViewType  GridViewType;
  typedef typename AssemblerType::EntityType    EntityType;

  LocalVolumeVectorAssemblerWrapper(const DS::PerThreadValue< const TestSpaceType >& space,
                                    const Stuff::Grid::ApplyOn::WhichEntity< GridViewType >* where,
                                    const LocalVolumeVectorAssembler& localAssembler,
                                    VectorType& vector)
    : TmpVectorsProvider(localAssembler.numTmpObjectsRequired(), space->mapper().maxNumDofs())
    , space_(space)
    , where_(where)
    , localVectorAssembler_(localAssembler)
    , vector_(vector)
  {}

  virtual ~LocalVolumeVectorAssemblerWrapper() {}

  virtual bool apply_on(const GridViewType& gv, const EntityType& entity) const override final
  {
    return where_->apply_on(gv, entity);
  }

  virtual void apply_local(const EntityType& entity) override final
  {
    localVectorAssembler_.assembleLocal(*space_, entity, vector_, this->vectors(), this->indices());
  }

private:
  const DS::PerThreadValue< const TestSpaceType >& space_;
  const std::unique_ptr< const Stuff::Grid::ApplyOn::WhichEntity< GridViewType > > where_;
  const LocalVolumeVectorAssembler& localVectorAssembler_;
  VectorType& vector_;
}; // class LocalVolumeVectorAssemblerWrapper


template< class AssemblerType, class LocalFaceVectorAssembler, class VectorType >
class LocalFaceVectorAssemblerWrapper
  : public Stuff::Grid::internal::Codim1Object< typename AssemblerType::GridViewType >
  , DSC::TmpVectorsStorage< typename AssemblerType::TestSpaceType::RangeFieldType >
{
  typedef DSC::TmpVectorsStorage< typename AssemblerType::TestSpaceType::RangeFieldType > TmpVectorsProvider;
public:
  typedef typename AssemblerType::TestSpaceType                                          TestSpaceType;
  typedef typename AssemblerType::GridViewType                                           GridViewType;
  typedef typename AssemblerType::EntityType                                             EntityType;
  typedef typename Stuff::Grid::internal::Codim1Object< GridViewType >::IntersectionType IntersectionType;

  LocalFaceVectorAssemblerWrapper(const DS::PerThreadValue< const TestSpaceType >& space,
                                  const Stuff::Grid::ApplyOn::WhichIntersection< GridViewType >* where,
                                  const LocalFaceVectorAssembler& localAssembler,
                                  VectorType& vector)
    : TmpVectorsProvider(localAssembler.numTmpObjectsRequired(), space->mapper().maxNumDofs())
    , space_(space)
    , where_(where)
    , localVectorAssembler_(localAssembler)
    , vector_(vector)
  {}

  virtual ~LocalFaceVectorAssemblerWrapper() {}

  virtual bool apply_on(const GridViewType& gv, const IntersectionType& intersection) const override final
  {
    return where_->apply_on(gv, intersection);
  }

  virtual void apply_local(const IntersectionType& intersection,
                           const EntityType& /*inside_entity*/,
                           const EntityType& /*outside_entity*/) override final
  {
    localVectorAssembler_.assembleLocal(*space_, intersection, vector_, this->vectors(), this->indices());
  }

private:
  const DS::PerThreadValue< const TestSpaceType >& space_;
  const std::unique_ptr< const Stuff::Grid::ApplyOn::WhichIntersection< GridViewType > > where_;
  const LocalFaceVectorAssembler& localVectorAssembler_;
  VectorType& vector_;
}; // class LocalFaceVectorAssemblerWrapper


} // namespace internal
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_WRAPPER_HH
