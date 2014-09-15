// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_ASSEMBLER_SYSTEM_HH
#define DUNE_GDT_ASSEMBLER_SYSTEM_HH

#include <type_traits>
#include <memory>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/constraints.hh>

#if DUNE_VERSION_NEWER(DUNE_COMMON,3,9) //&& HAVE_TBB //EXADUNE
# include <dune/grid/utility/partitioning/seedlist.hh>
#endif

#include "local/codim0.hh"
#include "local/codim1.hh"
#include "gridwalker.hh"
#include "tmp-storage.hh"

namespace Dune {
namespace GDT {


template< class TestSpaceImp,
          class GridViewImp = typename TestSpaceImp::GridViewType,
          class AnsatzSpaceImp = TestSpaceImp >
class SystemAssembler
  : public GridWalker< GridViewImp >
{
  static_assert(std::is_base_of< SpaceInterface< typename TestSpaceImp::Traits >, TestSpaceImp >::value,
                "TestSpaceImp has to be derived from SpaceInterface!");
  static_assert(std::is_base_of< SpaceInterface< typename AnsatzSpaceImp::Traits >, AnsatzSpaceImp >::value,
                "AnsatzSpaceImp has to be derived from SpaceInterface!");
  typedef GridWalker< GridViewImp > BaseType;
public:
  typedef TestSpaceImp    TestSpaceType;
  typedef AnsatzSpaceImp  AnsatzSpaceType;

  typedef typename BaseType::GridViewType     GridViewType;
  typedef typename BaseType::EntityType       EntityType;
  typedef typename BaseType::IntersectionType IntersectionType;
  typedef typename BaseType::BoundaryInfoType BoundaryInfoType;

private:
  typedef typename TestSpaceType::RangeFieldType RangeFieldType;

public:
  SystemAssembler(const TestSpaceType& test, const AnsatzSpaceType& ansatz, const GridViewType& grid_view)
    : BaseType(grid_view)
    , test_space_(test)
    , ansatz_space_(ansatz)
  {}

  SystemAssembler(const TestSpaceType& test, const AnsatzSpaceType& ansatz)
    : BaseType(*(test.grid_view()))
    , test_space_(test)
    , ansatz_space_(ansatz)
  {}

  SystemAssembler(const TestSpaceType& test)
    : BaseType(*(test.grid_view()))
    , test_space_(test)
    , ansatz_space_(test_space_)
  {}

  SystemAssembler(const TestSpaceType& test, const GridViewType& grid_view)
    : BaseType(grid_view)
    , test_space_(test)
    , ansatz_space_(test_space_)
  {}

  const TestSpaceType& test_space() const
  {
    return test_space_;
  }

  const AnsatzSpaceType& ansatz_space() const
  {
    return ansatz_space_;
  }

  using BaseType::add;

private:
  template< class ConstraintsType, class MatrixType >
  class LocalMatrixConstraintsWrapper
    : public BaseType::Codim0Object
  {
  public:
    LocalMatrixConstraintsWrapper(const TestSpaceType& t_space,
                                  const AnsatzSpaceType& a_space,
                                  const ApplyOn::WhichEntity< GridViewType >* where,
                                  ConstraintsType& constraints,
                                  MatrixType& matrix)
      : t_space_(t_space)
      , a_space_(a_space)
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
      t_space_.local_constraints(a_space_, entity, constraints_);
      for (size_t ii = 0; ii < constraints_.rows(); ++ii) {
        const size_t row = constraints_.globalRow(ii);
        for (size_t jj = 0; jj < constraints_.cols(); ++jj) {
          matrix_.set_entry(row, constraints_.globalCol(jj), constraints_.value(ii, jj));
        }
      }
    } // ... apply_local(...)

  private:
    const TestSpaceType& t_space_;
    const AnsatzSpaceType& a_space_;
    std::unique_ptr< const ApplyOn::WhichEntity< GridViewType > > where_;
    ConstraintsType& constraints_;
    MatrixType& matrix_;
  }; // class LocalMatrixConstraintsWrapper


  template< class ConstraintsType, class VectorType >
  class LocalVectorConstraintsWrapper
      : public BaseType::Codim0Object
  {
  public:
    LocalVectorConstraintsWrapper(const TestSpaceType& t_space,
                                  const ApplyOn::WhichEntity< GridViewType >* where,
                                  ConstraintsType& constraints,
                                  VectorType& vector)
      : t_space_(t_space)
      , where_(where)
      , constraints_(constraints)
      , vector_(vector)
    {}

    virtual ~LocalVectorConstraintsWrapper() {}

    virtual bool apply_on(const GridViewType& gv, const EntityType& entity) const DS_OVERRIDE DS_FINAL
    {
      return where_->apply_on(gv, entity);
    }

    virtual void apply_local(const EntityType& entity) DS_OVERRIDE DS_FINAL
    {
      t_space_.local_constraints(entity, constraints_);
      for (size_t ii = 0; ii < constraints_.rows(); ++ii) {
        vector_.set_entry(constraints_.globalRow(ii), RangeFieldType(0));
      }
    }

  private:
    const TestSpaceType& t_space_;
    std::unique_ptr< const ApplyOn::WhichEntity< GridViewType > > where_;
    ConstraintsType& constraints_;
    VectorType& vector_;
  }; // class LocalVectorConstraintsWrapper


  template< class LocalVolumeMatrixAssembler, class MatrixType >
  class LocalVolumeMatrixAssemblerWrapper
    : public BaseType::Codim0Object
    , TmpStorageProvider::Matrices< RangeFieldType >
  {
    typedef TmpStorageProvider::Matrices< RangeFieldType > TmpMatricesProvider;
  public:
    LocalVolumeMatrixAssemblerWrapper(const TestSpaceType& t_space,
                                      const AnsatzSpaceType& a_space,
                                      const ApplyOn::WhichEntity< GridViewType >* where,
                                      const LocalVolumeMatrixAssembler& localAssembler,
                                      MatrixType& matrix)
      : TmpMatricesProvider(localAssembler.numTmpObjectsRequired(),
                            t_space.mapper().maxNumDofs(),
                            a_space.mapper().maxNumDofs())
      , t_space_(t_space)
      , a_space_(a_space)
      , where_(where)
      , localMatrixAssembler_(localAssembler)
      , matrix_(matrix)
    {}

    virtual ~LocalVolumeMatrixAssemblerWrapper() {}

    virtual bool apply_on(const GridViewType& gv, const EntityType& entity) const DS_OVERRIDE DS_FINAL
    {
      return where_->apply_on(gv, entity);
    }

    virtual void apply_local(const EntityType& entity) DS_OVERRIDE DS_FINAL
    {
      localMatrixAssembler_.assembleLocal(t_space_, a_space_, entity, matrix_, this->matrices(), this->indices());
    }

  private:
    const TestSpaceType& t_space_;
    const AnsatzSpaceType& a_space_;
    std::unique_ptr< const ApplyOn::WhichEntity< GridViewType > > where_;
    const LocalVolumeMatrixAssembler& localMatrixAssembler_;
    MatrixType& matrix_;
  }; // class LocalVolumeMatrixAssemblerWrapper


  template< class LocalFaceMatrixAssembler, class MatrixType >
  class LocalFaceMatrixAssemblerWrapper
    : public BaseType::Codim1Object
    , TmpStorageProvider::Matrices< RangeFieldType >
  {
    typedef TmpStorageProvider::Matrices< RangeFieldType > TmpMatricesProvider;
  public:
    LocalFaceMatrixAssemblerWrapper(const TestSpaceType& t_space,
                                        const AnsatzSpaceType& a_space,
                                        const ApplyOn::WhichIntersection< GridViewType >* where,
                                        const LocalFaceMatrixAssembler& localAssembler,
                                        MatrixType& matrix)
      : TmpMatricesProvider(localAssembler.numTmpObjectsRequired(),
                            t_space.mapper().maxNumDofs(),
                            a_space.mapper().maxNumDofs())
      , t_space_(t_space)
      , a_space_(a_space)
      , where_(where)
      , localMatrixAssembler_(localAssembler)
      , matrix_(matrix)
    {}

    virtual ~LocalFaceMatrixAssemblerWrapper() {}

    virtual bool apply_on(const GridViewType& gv, const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
    {
      return where_->apply_on(gv, intersection);
    }

    virtual void apply_local(const IntersectionType& intersection,
                             const EntityType& /*inside_entity*/,
                             const EntityType& /*outside_entity*/) DS_OVERRIDE DS_FINAL
    {
      localMatrixAssembler_.assembleLocal(t_space_, a_space_, intersection, matrix_, this->matrices(), this->indices());
    }

  private:
    const TestSpaceType& t_space_;
    const AnsatzSpaceType& a_space_;
    std::unique_ptr< const ApplyOn::WhichIntersection< GridViewType > > where_;
    const LocalFaceMatrixAssembler& localMatrixAssembler_;
    MatrixType& matrix_;
  }; // class LocalFaceMatrixAssemblerWrapper


  template< class LocalVolumeVectorAssembler, class VectorType >
  class LocalVolumeVectorAssemblerWrapper
    : public BaseType::Codim0Object
    , TmpStorageProvider::Vectors< RangeFieldType >
  {
    typedef TmpStorageProvider::Vectors< RangeFieldType > TmpVectorsProvider;
  public:
    LocalVolumeVectorAssemblerWrapper(const TestSpaceType& space,
                                      const ApplyOn::WhichEntity< GridViewType >* where,
                                      const LocalVolumeVectorAssembler& localAssembler,
                                      VectorType& vector)
      : TmpVectorsProvider(localAssembler.numTmpObjectsRequired(), space.mapper().maxNumDofs())
      , space_(space)
      , where_(where)
      , localVectorAssembler_(localAssembler)
      , vector_(vector)
    {}

    virtual ~LocalVolumeVectorAssemblerWrapper() {}

    virtual bool apply_on(const GridViewType& gv, const EntityType& entity) const DS_OVERRIDE DS_FINAL
    {
      return where_->apply_on(gv, entity);
    }

    virtual void apply_local(const EntityType& entity) DS_OVERRIDE DS_FINAL
    {
      localVectorAssembler_.assembleLocal(space_, entity, vector_, this->vectors(), this->indices());
    }

  private:
    const TestSpaceType& space_;
    std::unique_ptr< const ApplyOn::WhichEntity< GridViewType > > where_;
    const LocalVolumeVectorAssembler& localVectorAssembler_;
    VectorType& vector_;
  }; // class LocalVolumeVectorAssemblerWrapper


  template< class LocalFaceVectorAssembler, class VectorType >
  class LocalFaceVectorAssemblerWrapper
    : public BaseType::Codim1Object
    , TmpStorageProvider::Vectors< RangeFieldType >
  {
    typedef TmpStorageProvider::Vectors< RangeFieldType > TmpVectorsProvider;
  public:
    LocalFaceVectorAssemblerWrapper(const TestSpaceType& space,
                                    const ApplyOn::WhichIntersection< GridViewType >* where,
                                    const LocalFaceVectorAssembler& localAssembler,
                                    VectorType& vector)
      : TmpVectorsProvider(localAssembler.numTmpObjectsRequired(), space.mapper().maxNumDofs())
      , space_(space)
      , where_(where)
      , localVectorAssembler_(localAssembler)
      , vector_(vector)
    {}

    virtual ~LocalFaceVectorAssemblerWrapper() {}

    virtual bool apply_on(const GridViewType& gv, const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
    {
      return where_->apply_on(gv, intersection);
    }

    virtual void apply_local(const IntersectionType& intersection,
                             const EntityType& /*inside_entity*/,
                             const EntityType& /*outside_entity*/) DS_OVERRIDE DS_FINAL
    {
      localVectorAssembler_.assembleLocal(space_, intersection, vector_, this->vectors(), this->indices());
    }

  private:
    const TestSpaceType& space_;
    std::unique_ptr< const ApplyOn::WhichIntersection< GridViewType > > where_;
    const LocalFaceVectorAssembler& localVectorAssembler_;
    VectorType& vector_;
  }; // class LocalFaceVectorAssemblerWrapper

public:
  template< class ConstraintsType, class M >
  void add(ConstraintsType& constraints,
           Dune::Stuff::LA::MatrixInterface< M >& matrix,
           const ApplyOn::WhichEntity< GridViewType >* where = new ApplyOn::AllEntities< GridViewType >())
  {
    typedef typename M::derived_type MatrixType;
    MatrixType& matrix_imp = static_cast< MatrixType& >(matrix);
    assert(matrix_imp.rows() == test_space_.mapper().size());
    assert(matrix_imp.cols() == ansatz_space_.mapper().size());
    typedef LocalMatrixConstraintsWrapper< ConstraintsType, MatrixType > WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, ansatz_space_, where, constraints, matrix_imp));
  }  // ... add(...)

  template< class ConstraintsType, class V >
  void add(ConstraintsType& constraints,
           Dune::Stuff::LA::VectorInterface< V >& vector,
           const ApplyOn::WhichEntity< GridViewType >* where = new ApplyOn::AllEntities< GridViewType >())
  {
    typedef typename V::derived_type VectorType;
    VectorType& vector_imp = static_cast< VectorType& >(vector);
    assert(vector_imp.size() == test_space_.mapper().size());
    typedef LocalVectorConstraintsWrapper< ConstraintsType, VectorType > WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, where, constraints, vector_imp));
  }  // ... add(...)

  template< class L, class M >
  void add(const LocalAssembler::Codim0Matrix< L >& local_assembler,
           Dune::Stuff::LA::MatrixInterface< M >& matrix,
           const ApplyOn::WhichEntity< GridViewType >* where = new ApplyOn::AllEntities< GridViewType >())
  {
    typedef typename M::derived_type MatrixType;
    MatrixType& matrix_imp = static_cast< MatrixType& >(matrix);
    assert(matrix_imp.rows() == test_space_.mapper().size());
    assert(matrix_imp.cols() == ansatz_space_.mapper().size());
    typedef LocalVolumeMatrixAssemblerWrapper< LocalAssembler::Codim0Matrix< L >, MatrixType > WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix_imp));
  }  // ... add(...)

  template< class Codim0Assembler, class M >
  void add_codim0_assembler(const Codim0Assembler& local_assembler,
           Dune::Stuff::LA::MatrixInterface< M >& matrix,
           const ApplyOn::WhichEntity< GridViewType >* where = new ApplyOn::AllEntities< GridViewType >())
  {
    typedef typename M::derived_type MatrixType;
    MatrixType& matrix_imp = static_cast< MatrixType& >(matrix);
    assert(matrix_imp.rows() == test_space_.mapper().size());
    assert(matrix_imp.cols() == ansatz_space_.mapper().size());
    typedef LocalVolumeMatrixAssemblerWrapper< Codim0Assembler, MatrixType > WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix_imp));
  }  // ... add(...)

  template< class Codim0Assembler, class V >
  void add_codim0_assembler(const Codim0Assembler& local_assembler,
           Dune::Stuff::LA::VectorInterface< V >& vector,
           const ApplyOn::WhichEntity< GridViewType >* where = new ApplyOn::AllEntities< GridViewType >())
  {
    typedef typename V::derived_type VectorType;
    VectorType& vector_imp = static_cast< VectorType& >(vector);
    assert(vector_imp.size() == test_space_.mapper().size());
    typedef LocalVolumeVectorAssemblerWrapper< Codim0Assembler, VectorType > WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, where, local_assembler, vector_imp));
  }  // ... add(...)

  template< class L, class M >
  void add(const LocalAssembler::Codim1CouplingMatrix< L >& local_assembler,
           Dune::Stuff::LA::MatrixInterface< M >& matrix,
           const ApplyOn::WhichIntersection< GridViewType >* where = new ApplyOn::AllIntersections< GridViewType >())
  {
    typedef typename M::derived_type MatrixType;
    MatrixType& matrix_imp = static_cast< MatrixType& >(matrix);
    assert(matrix_imp.rows() == test_space_.mapper().size());
    assert(matrix_imp.cols() == ansatz_space_.mapper().size());
    typedef LocalFaceMatrixAssemblerWrapper< LocalAssembler::Codim1CouplingMatrix< L >, MatrixType > WrapperType;
    this->codim1_functors_.emplace_back(new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix_imp));
  }  // ... add(...)

  template< class L, class M >
  void add(const LocalAssembler::Codim1BoundaryMatrix< L >& local_assembler,
           Dune::Stuff::LA::MatrixInterface< M >& matrix,
           const ApplyOn::WhichIntersection< GridViewType >* where = new ApplyOn::AllIntersections< GridViewType >())
  {
    typedef typename M::derived_type MatrixType;
    MatrixType& matrix_imp = static_cast< MatrixType& >(matrix);
    assert(matrix_imp.rows() == test_space_.mapper().size());
    assert(matrix_imp.cols() == ansatz_space_.mapper().size());
    typedef LocalFaceMatrixAssemblerWrapper< LocalAssembler::Codim1BoundaryMatrix< L >, MatrixType > WrapperType;
    this->codim1_functors_.emplace_back(new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix_imp));
  }  // ... add(...)

  template< class L, class V >
  void add(const LocalAssembler::Codim0Vector< L >& local_assembler,
           Dune::Stuff::LA::VectorInterface< V >& vector,
           const ApplyOn::WhichEntity< GridViewType >* where = new ApplyOn::AllEntities< GridViewType >())
  {
    typedef typename V::derived_type VectorType;
    VectorType& vector_imp = static_cast< VectorType& >(vector);
    assert(vector_imp.size() == test_space_.mapper().size());
    typedef LocalVolumeVectorAssemblerWrapper< LocalAssembler::Codim0Vector< L >, VectorType > WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, where, local_assembler, vector_imp));
  }  // ... add(...)

  template< class L, class V >
  void add(const LocalAssembler::Codim1Vector< L >& local_assembler,
           Dune::Stuff::LA::VectorInterface< V >& vector,
           const ApplyOn::WhichIntersection< GridViewType >* where = new ApplyOn::AllIntersections< GridViewType >())
  {
    typedef typename V::derived_type VectorType;
    VectorType& vector_imp = static_cast< VectorType& >(vector);
    assert(vector_imp.size() == test_space_.mapper().size());
    typedef LocalFaceVectorAssemblerWrapper< LocalAssembler::Codim1Vector< L >, VectorType > WrapperType;
    this->codim1_functors_.emplace_back(new WrapperType(test_space_, where, local_assembler, vector_imp));
  }  // ... add(...)

  void assemble(const bool clear_stack = true)
  {
    this->walk(clear_stack);
  }

#if 1 //HAVE_TBB
  void tbb_assemble(const bool clear_stack = true)
  {
    struct IndexSetPartitioner {
      typedef typename GridViewType::IndexSet IndexSetType;
      IndexSetPartitioner(const IndexSetType& index_set)
        : index_set_(index_set)
      {}

      std::size_t partition(const EntityType &e) const
      {
        return index_set_.index(e);
      }

      std::size_t partitions() const
      {
        return index_set_.size(0);
      }

    private:
      const IndexSetType& index_set_;
    };
    IndexSetPartitioner partioner(this->grid_view_.indexSet());
    SeedListPartitioning<typename GridViewType::Grid, 0> partitioning(this->grid_view_, partioner);
    this->tbb_walk(partitioning, clear_stack);
  }
#endif

private:
  const TestSpaceType& test_space_;
  const AnsatzSpaceType& ansatz_space_;
}; // class SystemAssembler


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_SYSTEM_HH
