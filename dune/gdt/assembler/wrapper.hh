#ifndef DUNE_GDT_ASSEMBLER_WRAPPER_HH
#define DUNE_GDT_ASSEMBLER_WRAPPER_HH

#include <type_traits>

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


template< class TestSpaceType, class AnsatzSpaceType, class GridViewType, class ConstraintsType, class MatrixType >
class LocalMatrixConstraintsWrapper
  : public Stuff::Grid::internal::Codim0Object< GridViewType >
{
  static_assert(std::is_base_of< SpaceInterface< typename TestSpaceType::Traits >, TestSpaceType >::value, "");
  static_assert(std::is_base_of< SpaceInterface< typename AnsatzSpaceType::Traits >, AnsatzSpaceType >::value, "");
  static_assert(std::is_base_of< Stuff::LA::MatrixInterface< typename MatrixType::Traits >, MatrixType >::value, "");
  static_assert(std::is_base_of< Spaces::ConstraintsInterface< typename ConstraintsType::Traits,
                                                               typename ConstraintsType::ValueType >,
                                 ConstraintsType >::value, "");
  typedef Stuff::Grid::internal::Codim0Object< GridViewType > BaseType;
public:
  typedef typename BaseType::EntityType EntityType;

  LocalMatrixConstraintsWrapper(const TestSpaceType& test_space,
                                const AnsatzSpaceType& ansatz_space,
                                const Stuff::Grid::ApplyOn::WhichEntity< GridViewType >* where,
                                ConstraintsType& constraints,
                                MatrixType& matrix)
    : test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , where_(where)
    , constraints_(constraints)
    , matrix_(matrix)
  {}

  virtual ~LocalMatrixConstraintsWrapper() {}

  virtual bool apply_on(const GridViewType& gv, const EntityType& entity) const DS_OVERRIDE DS_FINAL
  {
    return where_->apply_on(gv, entity);
  }

  virtual void apply_local(const EntityType& entity) DS_OVERRIDE DS_FINAL
  {
    test_space_.local_constraints(ansatz_space_, entity, constraints_);
    for (size_t ii = 0; ii < constraints_.rows(); ++ii) {
      const size_t row = constraints_.global_row(ii);
      for (size_t jj = 0; jj < constraints_.cols(); ++jj) {
        matrix_.set_entry(row, constraints_.global_col(jj), constraints_.value(ii, jj));
      }
    }
  } // ... apply_local(...)

private:
  const TestSpaceType& test_space_;
  const AnsatzSpaceType& ansatz_space_;
  std::unique_ptr< const Stuff::Grid::ApplyOn::WhichEntity< GridViewType > > where_;
  ConstraintsType& constraints_;
  MatrixType& matrix_;
}; // class LocalMatrixConstraintsWrapper


template< class AssemblerType, class ConstraintsType, class VectorType >
class LocalVectorConstraintsWrapper
    : public Stuff::Grid::internal::Codim0Object<typename AssemblerType::GridViewType>
{
public:
  LocalVectorConstraintsWrapper(const typename AssemblerType::TestSpaceType& test_space,
                                const Stuff::Grid::ApplyOn::WhichEntity< typename AssemblerType::GridViewType >* where,
                                ConstraintsType& constraints,
                                VectorType& vector)
    : test_space_(test_space)
    , where_(where)
    , constraints_(constraints)
    , vector_(vector)
  {}

  virtual ~LocalVectorConstraintsWrapper() {}

  virtual bool apply_on(const typename AssemblerType::GridViewType& gv,
                        const typename AssemblerType::EntityType& entity) const DS_OVERRIDE DS_FINAL
  {
    return where_->apply_on(gv, entity);
  }

  virtual void apply_local(const typename AssemblerType::EntityType& entity) DS_OVERRIDE DS_FINAL
  {
    test_space_.local_constraints(entity, constraints_);
    for (size_t ii = 0; ii < constraints_.rows(); ++ii) {
      vector_.set_entry(constraints_.global_row(ii), typename AssemblerType::TestSpaceType::RangeFieldType(0));
    }
  }

private:
  const typename AssemblerType::TestSpaceType& test_space_;
  std::unique_ptr< const Stuff::Grid::ApplyOn::WhichEntity< typename AssemblerType::GridViewType > > where_;
  ConstraintsType& constraints_;
  VectorType& vector_;
}; // class LocalVectorConstraintsWrapper


template< class AssemblerType, class LocalVolumeMatrixAssembler, class MatrixType >
class LocalVolumeMatrixAssemblerWrapper
  : public Stuff::Grid::internal::Codim0Object<typename AssemblerType::GridViewType>
  , TmpStorageProvider::Matrices< typename AssemblerType::TestSpaceType::RangeFieldType >
{
  typedef TmpStorageProvider::Matrices< typename AssemblerType::TestSpaceType::RangeFieldType > TmpMatricesProvider;
public:
  LocalVolumeMatrixAssemblerWrapper(const typename AssemblerType::TestSpaceType& test_space,
                                    const typename AssemblerType::AnsatzSpaceType& ansatz_space,
                                    const Stuff::Grid::ApplyOn::WhichEntity< typename AssemblerType::GridViewType >* where,
                                    const LocalVolumeMatrixAssembler& localAssembler,
                                    MatrixType& matrix)
    : TmpMatricesProvider(localAssembler.numTmpObjectsRequired(),
                          test_space.mapper().maxNumDofs(),
                          ansatz_space.mapper().maxNumDofs())
    , test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , where_(where)
    , localMatrixAssembler_(localAssembler)
    , matrix_(matrix)
  {}

  virtual ~LocalVolumeMatrixAssemblerWrapper() {}

  virtual bool apply_on(const typename AssemblerType::GridViewType& gv,
                        const typename AssemblerType::EntityType& entity) const DS_OVERRIDE DS_FINAL
  {
    return where_->apply_on(gv, entity);
  }

  virtual void apply_local(const typename AssemblerType::EntityType& entity) DS_OVERRIDE DS_FINAL
  {
    localMatrixAssembler_.assembleLocal(test_space_, ansatz_space_, entity, matrix_, this->matrices(), this->indices());
  }

private:
  const typename AssemblerType::TestSpaceType& test_space_;
  const typename AssemblerType::AnsatzSpaceType& ansatz_space_;
  std::unique_ptr< const Stuff::Grid::ApplyOn::WhichEntity< typename AssemblerType::GridViewType > > where_;
  const LocalVolumeMatrixAssembler& localMatrixAssembler_;
  MatrixType& matrix_;
}; // class LocalVolumeMatrixAssemblerWrapper


template< class AssemblerType, class LocalFaceMatrixAssembler, class MatrixType >
class LocalFaceMatrixAssemblerWrapper
  : public Stuff::Grid::internal::Codim1Object<typename AssemblerType::GridViewType>
  , TmpStorageProvider::Matrices< typename AssemblerType::TestSpaceType::RangeFieldType >
{
  typedef TmpStorageProvider::Matrices< typename AssemblerType::TestSpaceType::RangeFieldType > TmpMatricesProvider;
public:
  LocalFaceMatrixAssemblerWrapper(const typename AssemblerType::TestSpaceType& test_space,
                                      const typename AssemblerType::AnsatzSpaceType& ansatz_space,
                                      const Stuff::Grid::ApplyOn::WhichIntersection< typename AssemblerType::GridViewType >* where,
                                      const LocalFaceMatrixAssembler& localAssembler,
                                      MatrixType& matrix)
    : TmpMatricesProvider(localAssembler.numTmpObjectsRequired(),
                          test_space.mapper().maxNumDofs(),
                          ansatz_space.mapper().maxNumDofs())
    , test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , where_(where)
    , localMatrixAssembler_(localAssembler)
    , matrix_(matrix)
  {}

  virtual ~LocalFaceMatrixAssemblerWrapper() {}

  virtual bool apply_on(const typename AssemblerType::GridViewType& gv,
                        const typename AssemblerType::IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
  {
    return where_->apply_on(gv, intersection);
  }

  virtual void apply_local(const typename AssemblerType::IntersectionType& intersection,
                           const typename AssemblerType::EntityType& /*inside_entity*/,
                           const typename AssemblerType::EntityType& /*outside_entity*/) DS_OVERRIDE DS_FINAL
  {
    localMatrixAssembler_.assembleLocal(test_space_, ansatz_space_, intersection, matrix_, this->matrices(), this->indices());
  }

private:
  const typename AssemblerType::TestSpaceType& test_space_;
  const typename AssemblerType::AnsatzSpaceType& ansatz_space_;
  std::unique_ptr< const Stuff::Grid::ApplyOn::WhichIntersection< typename AssemblerType::GridViewType > > where_;
  const LocalFaceMatrixAssembler& localMatrixAssembler_;
  MatrixType& matrix_;
}; // class LocalFaceMatrixAssemblerWrapper


template< class AssemblerType, class LocalVolumeVectorAssembler, class VectorType >
class LocalVolumeVectorAssemblerWrapper
  : public Stuff::Grid::internal::Codim0Object<typename AssemblerType::GridViewType>
  , TmpStorageProvider::Vectors< typename AssemblerType::TestSpaceType::RangeFieldType >
{
  typedef TmpStorageProvider::Vectors< typename AssemblerType::TestSpaceType::RangeFieldType > TmpVectorsProvider;
public:
  LocalVolumeVectorAssemblerWrapper(const typename AssemblerType::TestSpaceType& space,
                                    const Stuff::Grid::ApplyOn::WhichEntity< typename AssemblerType::GridViewType >* where,
                                    const LocalVolumeVectorAssembler& localAssembler,
                                    VectorType& vector)
    : TmpVectorsProvider(localAssembler.numTmpObjectsRequired(), space.mapper().maxNumDofs())
    , space_(space)
    , where_(where)
    , localVectorAssembler_(localAssembler)
    , vector_(vector)
  {}

  virtual ~LocalVolumeVectorAssemblerWrapper() {}

  virtual bool apply_on(const typename AssemblerType::GridViewType& gv,
                        const typename AssemblerType::EntityType& entity) const DS_OVERRIDE DS_FINAL
  {
    return where_->apply_on(gv, entity);
  }

  virtual void apply_local(const typename AssemblerType::EntityType& entity) DS_OVERRIDE DS_FINAL
  {
    localVectorAssembler_.assembleLocal(space_, entity, vector_, this->vectors(), this->indices());
  }

private:
  const typename AssemblerType::TestSpaceType& space_;
  std::unique_ptr< const Stuff::Grid::ApplyOn::WhichEntity< typename AssemblerType::GridViewType > > where_;
  const LocalVolumeVectorAssembler& localVectorAssembler_;
  VectorType& vector_;
}; // class LocalVolumeVectorAssemblerWrapper


template< class AssemblerType, class LocalFaceVectorAssembler, class VectorType >
class LocalFaceVectorAssemblerWrapper
  : public Stuff::Grid::internal::Codim1Object<typename AssemblerType::GridViewType>
  , TmpStorageProvider::Vectors< typename AssemblerType::TestSpaceType::RangeFieldType >
{
  typedef TmpStorageProvider::Vectors< typename AssemblerType::TestSpaceType::RangeFieldType > TmpVectorsProvider;
public:
  LocalFaceVectorAssemblerWrapper(const typename AssemblerType::TestSpaceType& space,
                                  const Stuff::Grid::ApplyOn::WhichIntersection< typename AssemblerType::GridViewType >* where,
                                  const LocalFaceVectorAssembler& localAssembler,
                                  VectorType& vector)
    : TmpVectorsProvider(localAssembler.numTmpObjectsRequired(), space.mapper().maxNumDofs())
    , space_(space)
    , where_(where)
    , localVectorAssembler_(localAssembler)
    , vector_(vector)
  {}

  virtual ~LocalFaceVectorAssemblerWrapper() {}

  virtual bool apply_on(const typename AssemblerType::GridViewType& gv,
                        const typename AssemblerType::IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
  {
    return where_->apply_on(gv, intersection);
  }

  virtual void apply_local(const typename AssemblerType::IntersectionType& intersection,
                           const typename AssemblerType::EntityType& /*inside_entity*/,
                           const typename AssemblerType::EntityType& /*outside_entity*/) DS_OVERRIDE DS_FINAL
  {
    localVectorAssembler_.assembleLocal(space_, intersection, vector_, this->vectors(), this->indices());
  }

private:
  const typename AssemblerType::TestSpaceType& space_;
  std::unique_ptr< const Stuff::Grid::ApplyOn::WhichIntersection< typename AssemblerType::GridViewType > > where_;
  const LocalFaceVectorAssembler& localVectorAssembler_;
  VectorType& vector_;
}; // class LocalFaceVectorAssemblerWrapper


} // namespace internal
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_WRAPPER_HH
