#ifndef DUNE_GDT_ASSEMBLER_WALKER_WRAPPER_HH
#define DUNE_GDT_ASSEMBLER_WALKER_WRAPPER_HH

#include <dune/gdt/assembler/apply_on.hh>
#include <dune/gdt/assembler/functors.hh>

namespace Dune {
namespace GDT {
namespace internal {


/** \defgroup movetostuff ``These classes should move to dune-stuff alongside the gridwalker asap!´´ */
/* \{ */


template <class GridViewType>
class Codim0Object : public Functor::Codim0<GridViewType>
{
  typedef typename GridViewType::template Codim<0>::Entity EntityType;

public:
  ~Codim0Object()
  {
  }
  virtual bool apply_on(const GridViewType& grid_view, const EntityType& entity) const = 0;
};


template <class GridViewType, class Codim0FunctorType>
class Codim0FunctorWrapper : public Codim0Object<GridViewType>
{
  typedef typename GridViewType::template Codim<0>::Entity EntityType;

public:
  Codim0FunctorWrapper(Codim0FunctorType& wrapped_functor, const ApplyOn::WhichEntity<GridViewType>* where)
    : wrapped_functor_(wrapped_functor)
    , where_(where)
  {
  }

  virtual ~Codim0FunctorWrapper()
  {
  }

  virtual void prepare() DS_OVERRIDE DS_FINAL
  {
    wrapped_functor_.prepare();
  }

  virtual bool apply_on(const GridViewType& grid_view, const EntityType& entity) const DS_OVERRIDE DS_FINAL
  {
    return where_->apply_on(grid_view, entity);
  }

  virtual void apply_local(const EntityType& entity) DS_OVERRIDE DS_FINAL
  {
    wrapped_functor_.apply_local(entity);
  }

  virtual void finalize() DS_OVERRIDE DS_FINAL
  {
    wrapped_functor_.finalize();
  }

private:
  Codim0FunctorType& wrapped_functor_;
  std::unique_ptr<const ApplyOn::WhichEntity<GridViewType>> where_;
}; // class Codim0FunctorWrapper


template <class GridViewType>
class Codim1Object : public Functor::Codim1<GridViewType>
{
  typedef typename GridViewType::Intersection IntersectionType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;

public:
  ~Codim1Object()
  {
  }
  virtual bool apply_on(const GridViewType& grid_view, const IntersectionType& intersection) const = 0;
};


template <class GridViewType, class Codim1FunctorType>
class Codim1FunctorWrapper : public Codim1Object<GridViewType>
{
  typedef typename GridViewType::Intersection IntersectionType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;

public:
  Codim1FunctorWrapper(Codim1FunctorType& wrapped_functor, const ApplyOn::WhichIntersection<GridViewType>* where)
    : wrapped_functor_(wrapped_functor)
    , where_(where)
  {
  }

  virtual void prepare() DS_OVERRIDE DS_FINAL
  {
    wrapped_functor_.prepare();
  }

  virtual bool apply_on(const GridViewType& grid_view, const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
  {
    return where_->apply_on(grid_view, intersection);
  }

  virtual void apply_local(const IntersectionType& intersection, const EntityType& inside_entity,
                           const EntityType& outside_entity) DS_OVERRIDE DS_FINAL
  {
    wrapped_functor_.apply_local(intersection, inside_entity, outside_entity);
  }

  virtual void finalize() DS_OVERRIDE DS_FINAL
  {
    wrapped_functor_.finalize();
  }

private:
  Codim1FunctorType& wrapped_functor_;
  std::unique_ptr<const ApplyOn::WhichIntersection<GridViewType>> where_;
}; // class Codim1FunctorWrapper


template <class GridViewType, class WalkerType>
class GridWalkerWrapper : public Codim0Object<GridViewType>, public Codim1Object<GridViewType>
{
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::Intersection IntersectionType;

public:
  GridWalkerWrapper(WalkerType& grid_walker, const ApplyOn::WhichEntity<GridViewType>* which_entities)
    : grid_walker_(grid_walker)
    , which_entities_(which_entities)
    , which_intersections_(new ApplyOn::AllIntersections<GridViewType>())
  {
  }

  GridWalkerWrapper(WalkerType& grid_walker, const ApplyOn::WhichIntersection<GridViewType>* which_intersections)
    : grid_walker_(grid_walker)
    , which_entities_(new ApplyOn::AllEntities<GridViewType>())
    , which_intersections_(which_intersections)
  {
  }

  virtual ~GridWalkerWrapper()
  {
  }

  virtual void prepare() DS_OVERRIDE DS_FINAL
  {
    grid_walker_.prepare();
  }

  virtual bool apply_on(const GridViewType& grid_view, const EntityType& entity) const DS_OVERRIDE DS_FINAL
  {
    return which_entities_->apply_on(grid_view, entity) && grid_walker_.apply_on(entity);
  }

  virtual bool apply_on(const GridViewType& grid_view, const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
  {
    return which_intersections_->apply_on(grid_view, intersection) && grid_walker_.apply_on(intersection);
  }

  virtual void apply_local(const EntityType& entity) DS_OVERRIDE DS_FINAL
  {
    grid_walker_.apply_local(entity);
  }

  virtual void apply_local(const IntersectionType& intersection, const EntityType& inside_entity,
                           const EntityType& outside_entity) DS_OVERRIDE DS_FINAL
  {
    grid_walker_.apply_local(intersection, inside_entity, outside_entity);
  }

  virtual void finalize() DS_OVERRIDE DS_FINAL
  {
    grid_walker_.finalize();
  }

private:
  WalkerType& grid_walker_;
  std::unique_ptr<const ApplyOn::WhichEntity<GridViewType>> which_entities_;
  std::unique_ptr<const ApplyOn::WhichIntersection<GridViewType>> which_intersections_;
}; // class GridWalkerWrapper


template <class GridViewType>
class Codim0LambdaWrapper : public Codim0Object<GridViewType>
{
  typedef Codim0Object<GridViewType> BaseType;

public:
  typedef typename BaseType::EntityType EntityType;
  typedef std::function<void(const EntityType&)> LambdaType;

  Codim0LambdaWrapper(LambdaType lambda, const ApplyOn::WhichEntity<GridViewType>* where)
    : lambda_(lambda)
    , where_(where)
  {
  }

  virtual ~Codim0LambdaWrapper()
  {
  }

  virtual bool apply_on(const GridViewType& grid_view, const EntityType& entity) const DS_OVERRIDE DS_FINAL
  {
    return where_->apply_on(grid_view, entity);
  }

  virtual void apply_local(const EntityType& entity) DS_OVERRIDE DS_FINAL
  {
    lambda_(entity);
  }

private:
  LambdaType lambda_;
  std::unique_ptr<const ApplyOn::WhichEntity<GridViewType>> where_;
}; // class Codim0LambdaWrapper


/* \} */ // movetostuff


template <class AssemblerType, class ConstraintsType, class MatrixType>
class LocalMatrixConstraintsWrapper : public Codim0Object<typename AssemblerType::GridViewType>
{
public:
  LocalMatrixConstraintsWrapper(const typename AssemblerType::TestSpaceType& t_space,
                                const typename AssemblerType::AnsatzSpaceType& a_space,
                                const ApplyOn::WhichEntity<typename AssemblerType::GridViewType>* where,
                                ConstraintsType& constraints, MatrixType& matrix)
    : t_space_(t_space)
    , a_space_(a_space)
    , where_(where)
    , constraints_(constraints)
    , matrix_(matrix)
  {
  }

  virtual ~LocalMatrixConstraintsWrapper()
  {
  }

  virtual bool apply_on(const typename AssemblerType::GridViewType& gv,
                        const typename AssemblerType::EntityType& entity) const DS_OVERRIDE DS_FINAL
  {
    return where_->apply_on(gv, entity);
  }

  virtual void apply_local(const typename AssemblerType::EntityType& entity) DS_OVERRIDE DS_FINAL
  {
    t_space_.local_constraints(a_space_, entity, constraints_);
    for (size_t ii = 0; ii < constraints_.rows(); ++ii) {
      const size_t row = constraints_.globalRow(ii);
      for (size_t jj = 0; jj < constraints_.cols(); ++jj) {
        matrix_.set_entry(row, constraints_.globalCol(jj), constraints_.value(ii, jj));
      }
    }
  } // ... apply_local(...)

private:
  const typename AssemblerType::TestSpaceType& t_space_;
  const typename AssemblerType::AnsatzSpaceType& a_space_;
  std::unique_ptr<const ApplyOn::WhichEntity<typename AssemblerType::GridViewType>> where_;
  ConstraintsType& constraints_;
  MatrixType& matrix_;
}; // class LocalMatrixConstraintsWrapper


template <class AssemblerType, class ConstraintsType, class VectorType>
class LocalVectorConstraintsWrapper : public Codim0Object<typename AssemblerType::GridViewType>
{
public:
  LocalVectorConstraintsWrapper(const typename AssemblerType::TestSpaceType& t_space,
                                const ApplyOn::WhichEntity<typename AssemblerType::GridViewType>* where,
                                ConstraintsType& constraints, VectorType& vector)
    : t_space_(t_space)
    , where_(where)
    , constraints_(constraints)
    , vector_(vector)
  {
  }

  virtual ~LocalVectorConstraintsWrapper()
  {
  }

  virtual bool apply_on(const typename AssemblerType::GridViewType& gv,
                        const typename AssemblerType::EntityType& entity) const DS_OVERRIDE DS_FINAL
  {
    return where_->apply_on(gv, entity);
  }

  virtual void apply_local(const typename AssemblerType::EntityType& entity) DS_OVERRIDE DS_FINAL
  {
    t_space_.local_constraints(entity, constraints_);
    for (size_t ii = 0; ii < constraints_.rows(); ++ii) {
      vector_.set_entry(constraints_.globalRow(ii), typename AssemblerType::TestSpaceType::RangeFieldType(0));
    }
  }

private:
  const typename AssemblerType::TestSpaceType& t_space_;
  std::unique_ptr<const ApplyOn::WhichEntity<typename AssemblerType::GridViewType>> where_;
  ConstraintsType& constraints_;
  VectorType& vector_;
}; // class LocalVectorConstraintsWrapper


template <class AssemblerType, class LocalVolumeMatrixAssembler, class MatrixType>
class LocalVolumeMatrixAssemblerWrapper
    : public Codim0Object<typename AssemblerType::GridViewType>,
      TmpStorageProvider::Matrices<typename AssemblerType::TestSpaceType::RangeFieldType>
{
  typedef TmpStorageProvider::Matrices<typename AssemblerType::TestSpaceType::RangeFieldType> TmpMatricesProvider;

public:
  LocalVolumeMatrixAssemblerWrapper(const typename AssemblerType::TestSpaceType& t_space,
                                    const typename AssemblerType::AnsatzSpaceType& a_space,
                                    const ApplyOn::WhichEntity<typename AssemblerType::GridViewType>* where,
                                    const LocalVolumeMatrixAssembler& localAssembler, MatrixType& matrix)
    : TmpMatricesProvider(localAssembler.numTmpObjectsRequired(), t_space.mapper().maxNumDofs(),
                          a_space.mapper().maxNumDofs())
    , t_space_(t_space)
    , a_space_(a_space)
    , where_(where)
    , localMatrixAssembler_(localAssembler)
    , matrix_(matrix)
  {
  }

  virtual ~LocalVolumeMatrixAssemblerWrapper()
  {
  }

  virtual bool apply_on(const typename AssemblerType::GridViewType& gv,
                        const typename AssemblerType::EntityType& entity) const DS_OVERRIDE DS_FINAL
  {
    return where_->apply_on(gv, entity);
  }

  virtual void apply_local(const typename AssemblerType::EntityType& entity) DS_OVERRIDE DS_FINAL
  {
    localMatrixAssembler_.assembleLocal(t_space_, a_space_, entity, matrix_, this->matrices(), this->indices());
  }

private:
  const typename AssemblerType::TestSpaceType& t_space_;
  const typename AssemblerType::AnsatzSpaceType& a_space_;
  std::unique_ptr<const ApplyOn::WhichEntity<typename AssemblerType::GridViewType>> where_;
  const LocalVolumeMatrixAssembler& localMatrixAssembler_;
  MatrixType& matrix_;
}; // class LocalVolumeMatrixAssemblerWrapper


template <class AssemblerType, class LocalFaceMatrixAssembler, class MatrixType>
class LocalFaceMatrixAssemblerWrapper
    : public Codim1Object<typename AssemblerType::GridViewType>,
      TmpStorageProvider::Matrices<typename AssemblerType::TestSpaceType::RangeFieldType>
{
  typedef TmpStorageProvider::Matrices<typename AssemblerType::TestSpaceType::RangeFieldType> TmpMatricesProvider;

public:
  LocalFaceMatrixAssemblerWrapper(const typename AssemblerType::TestSpaceType& t_space,
                                  const typename AssemblerType::AnsatzSpaceType& a_space,
                                  const ApplyOn::WhichIntersection<typename AssemblerType::GridViewType>* where,
                                  const LocalFaceMatrixAssembler& localAssembler, MatrixType& matrix)
    : TmpMatricesProvider(localAssembler.numTmpObjectsRequired(), t_space.mapper().maxNumDofs(),
                          a_space.mapper().maxNumDofs())
    , t_space_(t_space)
    , a_space_(a_space)
    , where_(where)
    , localMatrixAssembler_(localAssembler)
    , matrix_(matrix)
  {
  }

  virtual ~LocalFaceMatrixAssemblerWrapper()
  {
  }

  virtual bool apply_on(const typename AssemblerType::GridViewType& gv,
                        const typename AssemblerType::IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
  {
    return where_->apply_on(gv, intersection);
  }

  virtual void apply_local(const typename AssemblerType::IntersectionType& intersection,
                           const typename AssemblerType::EntityType& /*inside_entity*/,
                           const typename AssemblerType::EntityType& /*outside_entity*/) DS_OVERRIDE DS_FINAL
  {
    localMatrixAssembler_.assembleLocal(t_space_, a_space_, intersection, matrix_, this->matrices(), this->indices());
  }

private:
  const typename AssemblerType::TestSpaceType& t_space_;
  const typename AssemblerType::AnsatzSpaceType& a_space_;
  std::unique_ptr<const ApplyOn::WhichIntersection<typename AssemblerType::GridViewType>> where_;
  const LocalFaceMatrixAssembler& localMatrixAssembler_;
  MatrixType& matrix_;
}; // class LocalFaceMatrixAssemblerWrapper


template <class AssemblerType, class LocalVolumeVectorAssembler, class VectorType>
class LocalVolumeVectorAssemblerWrapper
    : public Codim0Object<typename AssemblerType::GridViewType>,
      TmpStorageProvider::Vectors<typename AssemblerType::TestSpaceType::RangeFieldType>
{
  typedef TmpStorageProvider::Vectors<typename AssemblerType::TestSpaceType::RangeFieldType> TmpVectorsProvider;

public:
  LocalVolumeVectorAssemblerWrapper(const typename AssemblerType::TestSpaceType& space,
                                    const ApplyOn::WhichEntity<typename AssemblerType::GridViewType>* where,
                                    const LocalVolumeVectorAssembler& localAssembler, VectorType& vector)
    : TmpVectorsProvider(localAssembler.numTmpObjectsRequired(), space.mapper().maxNumDofs())
    , space_(space)
    , where_(where)
    , localVectorAssembler_(localAssembler)
    , vector_(vector)
  {
  }

  virtual ~LocalVolumeVectorAssemblerWrapper()
  {
  }

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
  std::unique_ptr<const ApplyOn::WhichEntity<typename AssemblerType::GridViewType>> where_;
  const LocalVolumeVectorAssembler& localVectorAssembler_;
  VectorType& vector_;
}; // class LocalVolumeVectorAssemblerWrapper


template <class AssemblerType, class LocalFaceVectorAssembler, class VectorType>
class LocalFaceVectorAssemblerWrapper
    : public Codim1Object<typename AssemblerType::GridViewType>,
      TmpStorageProvider::Vectors<typename AssemblerType::TestSpaceType::RangeFieldType>
{
  typedef TmpStorageProvider::Vectors<typename AssemblerType::TestSpaceType::RangeFieldType> TmpVectorsProvider;

public:
  LocalFaceVectorAssemblerWrapper(const typename AssemblerType::TestSpaceType& space,
                                  const ApplyOn::WhichIntersection<typename AssemblerType::GridViewType>* where,
                                  const LocalFaceVectorAssembler& localAssembler, VectorType& vector)
    : TmpVectorsProvider(localAssembler.numTmpObjectsRequired(), space.mapper().maxNumDofs())
    , space_(space)
    , where_(where)
    , localVectorAssembler_(localAssembler)
    , vector_(vector)
  {
  }

  virtual ~LocalFaceVectorAssemblerWrapper()
  {
  }

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
  std::unique_ptr<const ApplyOn::WhichIntersection<typename AssemblerType::GridViewType>> where_;
  const LocalFaceVectorAssembler& localVectorAssembler_;
  VectorType& vector_;
}; // class LocalFaceVectorAssemblerWrapper


} // namespace internal
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_WALKER_WRAPPER_HH
