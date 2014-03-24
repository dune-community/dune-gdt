// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_ASSEMBLER_SYSTEM_HH
#define DUNE_GDT_ASSEMBLER_SYSTEM_HH

#include <type_traits>
//#include <vector>
#include <memory>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
//#include <dune/common/static_assert.hh>
//#include <dune/common/typetraits.hh>

//#include <dune/stuff/la/container/interfaces.hh>
//#include <dune/stuff/la/container/pattern.hh>
//#ifdef DUNE_STUFF_PROFILER_ENABLED
//# include <dune/stuff/common/profiler.hh>
//#endif
#include <dune/gdt/space/interface.hh>
#include <dune/gdt/space/constraints.hh>

#include "local/codim0.hh"
#include "local/codim1.hh"
#include "gridwalker.hh"

namespace Dune {
namespace GDT {


template <class TestSpaceImp, class GridViewImp = typename TestSpaceImp::GridViewType,
          class AnsatzSpaceImp                  = TestSpaceImp>
class SystemAssembler : public GridWalker<GridViewImp>
{
  static_assert(std::is_base_of<SpaceInterface<typename TestSpaceImp::Traits>, TestSpaceImp>::value,
                "TestSpaceImp has to be derived from SpaceInterface!");
  static_assert(std::is_base_of<SpaceInterface<typename AnsatzSpaceImp::Traits>, AnsatzSpaceImp>::value,
                "AnsatzSpaceImp has to be derived from SpaceInterface!");
  typedef GridWalker<GridViewImp> BaseType;

public:
  typedef TestSpaceImp TestSpaceType;
  typedef AnsatzSpaceImp AnsatzSpaceType;

  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::IntersectionType IntersectionType;
  typedef typename BaseType::BoundaryInfoType BoundaryInfoType;

private:
  typedef typename TestSpaceType::RangeFieldType RangeFieldType;
  typedef Dune::DynamicMatrix<RangeFieldType> LocalMatrixType;
  typedef Dune::DynamicVector<RangeFieldType> LocalVectorType;
  typedef std::vector<std::vector<LocalMatrixType>> LocalMatricesContainerType;
  typedef std::vector<std::vector<LocalVectorType>> LocalVectorsContainerType;
  typedef std::vector<Dune::DynamicVector<size_t>> IndicesContainer;

private:
#if 0
  class LocalCodim0MatrixAssemblerApplication
  {
  public:
    virtual ~LocalCodim0MatrixAssemblerApplication(){}

    virtual void apply(const TestSpaceType& /*_testSpace*/,
                       const AnsatzSpaceType& /*_ansatzSpace*/,
                       const EntityType& /*_entity*/,
                       LocalMatricesContainerType& /*_localMatricesContainer*/,
                       IndicesContainer& /*indicesContainer*/) const = 0;

    virtual std::vector< size_t > numTmpObjectsRequired() const = 0;
  };

  template< class L, class M >
  class LocalCodim0MatrixAssemblerWrapper
    : public LocalCodim0MatrixAssemblerApplication
  {
  public:
    LocalCodim0MatrixAssemblerWrapper(const LocalAssembler::Codim0Matrix< L >& localAssembler,
                                      Dune::Stuff::LA::MatrixInterface< M >& matrix)
      : localMatrixAssembler_(localAssembler)
      , matrix_(matrix)
    {}

    virtual void apply(const TestSpaceType& testSpace,
                       const AnsatzSpaceType& ansatzSpace,
                       const EntityType& entity,
                       LocalMatricesContainerType& localMatricesContainer,
                       IndicesContainer& indicesContainer) const DS_OVERRIDE
    {
      localMatrixAssembler_.assembleLocal(testSpace, ansatzSpace, entity, matrix_, localMatricesContainer, indicesContainer);
    }

    virtual std::vector< size_t > numTmpObjectsRequired() const DS_OVERRIDE
    {
      return localMatrixAssembler_.numTmpObjectsRequired();
    }

  private:
    const LocalAssembler::Codim0Matrix< L >& localMatrixAssembler_;
    Dune::Stuff::LA::MatrixInterface< M >& matrix_;
  }; // class LocalCodim0MatrixAssemblerWrapper

  class LocalCodim1MatrixAssemblerApplication
  {
  public:
    virtual ~LocalCodim1MatrixAssemblerApplication(){}

    virtual void apply(const GridViewType& /*gridView*/,
                       const TestSpaceType& /*_testSpace*/,
                       const AnsatzSpaceType& /*_ansatzSpace*/,
                       const IntersectionType& /*_intersection*/,
                       LocalMatricesContainerType& /*_localMatricesContainer*/,
                       IndicesContainer& /*indicesContainer*/) const = 0;

    virtual std::vector< size_t > numTmpObjectsRequired() const = 0;
  };

  template< class LocalAssemblerType, class AssembleOnFunctorType, class M >
  class LocalCodim1MatrixAssemblerWrapper
    : public LocalCodim1MatrixAssemblerApplication
  {
  public:
    LocalCodim1MatrixAssemblerWrapper(const LocalAssemblerType& localAssembler,
                                      const AssembleOnFunctorType functor,
                                      Dune::Stuff::LA::MatrixInterface< M >& matrix)
      : localMatrixAssembler_(localAssembler)
      , functor_(functor)
      , matrix_(matrix)
    {}

    virtual void apply(const GridViewType& gridView,
                       const TestSpaceType& testSpace,
                       const AnsatzSpaceType& ansatzSpace,
                       const IntersectionType& intersection,
                       LocalMatricesContainerType& localMatricesContainer,
                       IndicesContainer& indicesContainer) const DS_OVERRIDE
    {
      if (functor_.assembleOn(gridView, intersection))
        localMatrixAssembler_.assembleLocal(testSpace, ansatzSpace,
                                            intersection,
                                            matrix_,
                                            localMatricesContainer, indicesContainer);
    }

    virtual std::vector< size_t > numTmpObjectsRequired() const DS_OVERRIDE
    {
      return localMatrixAssembler_.numTmpObjectsRequired();
    }

  private:
    const LocalAssemblerType& localMatrixAssembler_;
    const AssembleOnFunctorType functor_;
    Dune::Stuff::LA::MatrixInterface< M >& matrix_;
  }; // class LocalCodim1MatrixAssemblerWrapper

  class LocalCodim0VectorAssemblerApplication
  {
  public:
    virtual ~LocalCodim0VectorAssemblerApplication(){}

    virtual void apply(const TestSpaceType& /*_testSpace*/,
                       const EntityType& /*_entity*/,
                       LocalVectorsContainerType& /*_localVectorsContainer*/,
                       Dune::DynamicVector< size_t >& /*indices*/) const = 0;

    virtual std::vector< size_t > numTmpObjectsRequired() const = 0;
  };

  template< class L, class V >
  class LocalCodim0VectorAssemblerWrapper
    : public LocalCodim0VectorAssemblerApplication
  {
  public:
    LocalCodim0VectorAssemblerWrapper(const LocalAssembler::Codim0Vector< L >& localAssembler,
                                      Dune::Stuff::LA::VectorInterface< V >& vector)
      : localVectorAssembler_(localAssembler)
      , vector_(vector)
    {}

    virtual void apply(const TestSpaceType& testSpace,
                       const EntityType& entity,
                       LocalVectorsContainerType& localVectorsContainer,
                       Dune::DynamicVector< size_t >& indices) const DS_OVERRIDE
    {
      localVectorAssembler_.assembleLocal(testSpace, entity, vector_, localVectorsContainer, indices);
    }

    virtual std::vector< size_t > numTmpObjectsRequired() const DS_OVERRIDE
    {
      return localVectorAssembler_.numTmpObjectsRequired();
    }

  private:
    const LocalAssembler::Codim0Vector< L >& localVectorAssembler_;
    Dune::Stuff::LA::VectorInterface< V >& vector_;
  }; // class LocalCodim0VectorAssemblerWrapper

  class LocalCodim1VectorAssemblerApplication
  {
  public:
    virtual ~LocalCodim1VectorAssemblerApplication(){}

    virtual void apply(const GridViewType& /*gridView*/,
                       const TestSpaceType& /*_testSpace*/,
                       const IntersectionType& /*_intersection*/,
                       LocalVectorsContainerType& /*_localVectorsContainer*/,
                       Dune::DynamicVector< size_t >& /*_indices*/) const = 0;

    virtual std::vector< size_t > numTmpObjectsRequired() const = 0;
  };

  template< class LocalAssemblerType, class AssembleOnFunctorType, class V >
  class LocalCodim1VectorAssemblerWrapper
    : public LocalCodim1VectorAssemblerApplication
  {
  public:
    LocalCodim1VectorAssemblerWrapper(const LocalAssemblerType& localAssembler,
                                      const AssembleOnFunctorType functor,
                                      Dune::Stuff::LA::VectorInterface< V >& vector)
      : localVectorAssembler_(localAssembler)
      , functor_(functor)
      , vector_(vector)
    {}

    virtual void apply(const GridViewType& gridView,
                       const TestSpaceType& testSpace,
                       const IntersectionType& intersection,
                       LocalVectorsContainerType& localVectorsContainer,
                       Dune::DynamicVector< size_t >& indices) const DS_OVERRIDE
    {
      if (functor_.assembleOn(gridView, intersection))
        localVectorAssembler_.assembleLocal(testSpace, intersection, vector_, localVectorsContainer, indices);
    }

    virtual std::vector< size_t > numTmpObjectsRequired() const DS_OVERRIDE
    {
      return localVectorAssembler_.numTmpObjectsRequired();
    }

  private:
    const LocalAssemblerType& localVectorAssembler_;
    const AssembleOnFunctorType functor_;
    Dune::Stuff::LA::VectorInterface< V >& vector_;
  }; // class LocalCodim1VectorAssemblerWrapper
#endif
public:
  SystemAssembler(const TestSpaceType& test, const AnsatzSpaceType& ansatz, const GridViewType& grid_view)
    : BaseType(grid_view)
    , test_space_(test)
    , ansatz_space_(ansatz)
  {
  }

  SystemAssembler(const TestSpaceType& test, const AnsatzSpaceType& ansatz)
    : BaseType(*(test.grid_view()))
    , test_space_(test)
    , ansatz_space_(ansatz)
  {
  }

  SystemAssembler(const TestSpaceType& test)
    : BaseType(*(test.grid_view()))
    , test_space_(test)
    , ansatz_space_(test_space_)
  {
  }

  SystemAssembler(const TestSpaceType& test, const GridViewType& grid_view)
    : BaseType(grid_view)
    , test_space_(test)
    , ansatz_space_(test_space_)
  {
  }

  const TestSpaceType& test_space() const
  {
    return test_space_;
  }

  const AnsatzSpaceType& ansatz_space() const
  {
    return ansatz_space_;
  }

  using BaseType::add;

#if 0
  template< class L, class M >
  void addLocalAssembler(const LocalAssembler::Codim0Matrix< L >& localAssembler,
                         Dune::Stuff::LA::MatrixInterface< M >& matrix)
  {
    assert(matrix.rows() == test_space_.mapper().size());
    assert(matrix.cols() == ansatz_space_.mapper().size());
    localCodim0MatrixAssemblers_.push_back(new LocalCodim0MatrixAssemblerWrapper< L, M >(localAssembler, matrix));
  }

  template< class L, class AssembleOnFunctorType, class M >
  void addLocalAssembler(const LocalAssembler::Codim1CouplingMatrix< L >& localAssembler,
                         const AssembleOnFunctorType functor,
                         Dune::Stuff::LA::MatrixInterface< M >& matrix)
  {
    dune_static_assert((Dune::IsBaseOf< AssembleOnFunctorInterface, AssembleOnFunctorType >::value),
                       "ERROR: AssembleOnFunctorType is not a AssembleOnFunctorInterface");
    assert(matrix.rows() == test_space_.mapper().size());
    assert(matrix.cols() == ansatz_space_.mapper().size());
    localCodim1MatrixAssemblers_.push_back(
          new LocalCodim1MatrixAssemblerWrapper<  LocalAssembler::Codim1CouplingMatrix< L >,
                                                  AssembleOnFunctorType,
                                                  M >(localAssembler, functor, matrix));
  }

  template< class L, class AssembleOnFunctorType, class M >
  void addLocalAssembler(const LocalAssembler::Codim1BoundaryMatrix< L >& localAssembler,
                         const AssembleOnFunctorType functor,
                         Dune::Stuff::LA::MatrixInterface< M >& matrix)
  {
    dune_static_assert((Dune::IsBaseOf< AssembleOnFunctorInterface, AssembleOnFunctorType >::value),
                       "ERROR: AssembleOnFunctorType is not a AssembleOnFunctorInterface");
    assert(matrix.rows() == test_space_.mapper().size());
    assert(matrix.cols() == ansatz_space_.mapper().size());
    localCodim1MatrixAssemblers_.push_back(
          new LocalCodim1MatrixAssemblerWrapper<  LocalAssembler::Codim1BoundaryMatrix< L >,
                                                  AssembleOnFunctorType,
                                                  M >(localAssembler, functor, matrix));
  }

  template< class L, class V >
  void addLocalAssembler(const LocalAssembler::Codim0Vector< L >& localAssembler,
                         Dune::Stuff::LA::VectorInterface< V >& vector)
  {
    assert(vector.size() == test_space_.mapper().size());
    localCodim0VectorAssemblers_.push_back(new LocalCodim0VectorAssemblerWrapper< L, V >(localAssembler, vector));
  }

  template< class L, class AssembleOnFunctorType, class V >
  void addLocalAssembler(const LocalAssembler::Codim1Vector< L >& localAssembler,
                         const AssembleOnFunctorType functor,
                         Dune::Stuff::LA::VectorInterface< V >& vector)
  {
    dune_static_assert((Dune::IsBaseOf< AssembleOnFunctorInterface, AssembleOnFunctorType >::value),
                       "ERROR: AssembleOnFunctorType is not a AssembleOnFunctorInterface");
    assert(vector.size() == test_space_.mapper().size());
    localCodim1VectorAssemblers_.push_back(
          new LocalCodim1VectorAssemblerWrapper<  LocalAssembler::Codim1Vector< L >,
                                                  AssembleOnFunctorType,
                                                  V >(localAssembler, functor, vector));
  }
#endif

  template <class ConstraintsType, class MatrixType>
  class LocalMatrixConstraintsWrapper : public BaseType::Codim0Object
  {
  public:
    LocalMatrixConstraintsWrapper(const TestSpaceType& t_space, const AnsatzSpaceType& a_space,
                                  const ApplyOn::WhichEntity<GridViewType>* where, ConstraintsType& constraints,
                                  MatrixType& matrix)
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

    virtual bool apply_on(const GridViewType& gv, const EntityType& entity) const DS_FINAL
    {
      return where_->apply_on(gv, entity);
    }

    virtual void apply_local(const EntityType& entity) DS_FINAL
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
    std::unique_ptr<const ApplyOn::WhichEntity<GridViewType>> where_;
    ConstraintsType& constraints_;
    MatrixType& matrix_;
  }; // class LocalMatrixConstraintsWrapper

  template <class ConstraintsType, class VectorType>
  class LocalVectorConstraintsWrapper : public BaseType::Codim0Object
  {
  public:
    LocalVectorConstraintsWrapper(const TestSpaceType& t_space, const ApplyOn::WhichEntity<GridViewType>* where,
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

    virtual bool apply_on(const GridViewType& gv, const EntityType& entity) const DS_FINAL
    {
      return where_->apply_on(gv, entity);
    }

    virtual void apply_local(const EntityType& entity) DS_FINAL
    {
      t_space_.local_constraints(entity, constraints_);
      for (size_t ii = 0; ii < constraints_.rows(); ++ii) {
        vector_.set_entry(constraints_.globalRow(ii), RangeFieldType(0));
      }
    }

  private:
    const TestSpaceType& t_space_;
    std::unique_ptr<const ApplyOn::WhichEntity<GridViewType>> where_;
    ConstraintsType& constraints_;
    VectorType& vector_;
  }; // class LocalVectorConstraintsWrapper

public:
  template <class ConstraintsType, class M>
  void add(ConstraintsType& constraints, Dune::Stuff::LA::MatrixInterface<M>& matrix,
           const ApplyOn::WhichEntity<GridViewType>* where = new ApplyOn::AllEntities<GridViewType>())
  {
    typedef typename M::derived_type MatrixType;
    MatrixType& matrix_imp = static_cast<MatrixType&>(matrix);
    assert(matrix_imp.rows() == test_space_.mapper().size());
    assert(matrix_imp.cols() == ansatz_space_.mapper().size());
    typedef LocalMatrixConstraintsWrapper<ConstraintsType, MatrixType> WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, ansatz_space_, where, constraints, matrix_imp));
  } // ... add(...)

  template <class ConstraintsType, class V>
  void add(ConstraintsType& constraints, Dune::Stuff::LA::VectorInterface<V>& vector,
           const ApplyOn::WhichEntity<GridViewType>* where = new ApplyOn::AllEntities<GridViewType>())
  {
    typedef typename V::derived_type VectorType;
    VectorType& vector_imp = static_cast<VectorType&>(vector);
    assert(vector_imp.size() == test_space_.mapper().size());
    typedef LocalVectorConstraintsWrapper<ConstraintsType, VectorType> WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(test_space_, where, constraints, vector_imp));
  } // ... add(...)

private:
  const TestSpaceType& test_space_;
  const AnsatzSpaceType& ansatz_space_;
}; // class SystemAssembler


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_SYSTEM_HH
