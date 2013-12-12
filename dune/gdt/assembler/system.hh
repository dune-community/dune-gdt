// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_ASSEMBLER_SYSTEM_HH
#define DUNE_GDT_ASSEMBLER_SYSTEM_HH

#include <vector>
#include <memory>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/typetraits.hh>

#include <dune/stuff/la/container/interfaces.hh>
#include <dune/stuff/la/container/pattern.hh>
#ifdef DUNE_STUFF_PROFILER_ENABLED
# include <dune/stuff/common/profiler.hh>
#endif
#include <dune/gdt/space/interface.hh>
#include <dune/gdt/space/constraints.hh>

#include "local/codim0.hh"
#include "local/codim1.hh"

namespace Dune {
namespace GDT {


template< class TestSpaceImp,
          class GridPartImp = typename TestSpaceImp::GridPartType,
          class AnsatzSpaceImp = TestSpaceImp >
class SystemAssembler
{
public:
  typedef GridPartImp GridPartType;
  typedef SpaceInterface< typename TestSpaceImp::Traits >     TestSpaceType;
  typedef SpaceInterface< typename AnsatzSpaceImp::Traits >   AnsatzSpaceType;
  typedef typename GridPartType::IntersectionType             IntersectionType;
  typedef Dune::Stuff::GridboundaryInterface< IntersectionType > GridBoundaryType;

private:
  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  typedef typename TestSpaceType::RangeFieldType RangeFieldType;
  typedef Dune::DynamicMatrix< RangeFieldType > LocalMatrixType;
  typedef Dune::DynamicVector< RangeFieldType > LocalVectorType;
  typedef std::vector< std::vector< LocalMatrixType > > LocalMatricesContainerType;
  typedef std::vector< std::vector< LocalVectorType > > LocalVectorsContainerType;
  typedef std::vector< Dune::DynamicVector< size_t > > IndicesContainer;

public:
  /**
   *  \brief Interface for functors to tell on which intersection to assemble on.
   *
   *  The derived class has to provide a method with the following signature:
   *  \code
  bool assembleOn(const GridPartType& gridPart, const IntersectionType& intersection) const
  {
    ...
  }
  \endcode
   */
  class AssembleOnFunctorInterface
  {
  public:
    virtual ~ AssembleOnFunctorInterface() {}
    virtual bool assembleOn(const GridPartType& /*gridPart*/, const IntersectionType& /*intersection*/) const = 0;
  };

  /**
   *  \brief Selects each inner intersection only once.
   *
   *  To decide if this in an inner intersection,
  \code
  intersection.neighbor() && !intersection.boundary()
  \endcode
   *  is used.
   */
  class AssembleOnInner
    : public AssembleOnFunctorInterface
  {
  public:
    virtual bool assembleOn(const GridPartType& /*gridPart*/, const IntersectionType& intersection) const DS_OVERRIDE
    {
      return intersection.neighbor() && !intersection.boundary();
    }
  }; // class AssembleOnInner

  /**
   *  \brief Selects each inner intersection only once.
   *
   *  To decide if this in an inner intersection,
  \code
  intersection.neighbor() && !intersection.boundary()
  \endcode
   *  is used, and true is returned, if the index of the inside() entity is smaller than the index of the outside()
   *  entity.
   */
  class AssembleOnInnerPrimally
    : public AssembleOnFunctorInterface
  {
  public:
    virtual bool assembleOn(const GridPartType& gridPart, const IntersectionType& intersection) const DS_OVERRIDE
    {
      if (intersection.neighbor() && !intersection.boundary()) {
        const auto insideEntityPtr = intersection.inside();
        const auto& insideEntity = *insideEntityPtr;
        const auto outsideNeighborPtr = intersection.outside();
        const auto& outsideNeighbor = *outsideNeighborPtr;
        return gridPart.indexSet().index(insideEntity) < gridPart.indexSet().index(outsideNeighbor);
      } else
        return false;
    }
  }; // class AssembleOnInnerPrimally

  class AssembleOnBoundary
    : public AssembleOnFunctorInterface
  {
  public:
    virtual bool assembleOn(const GridPartType& /*gridPart*/, const IntersectionType& intersection) const DS_OVERRIDE
    {
      return intersection.boundary();
    }
  }; // class AssembleOnBoundary

  class AssembleOnDirichlet
    : public AssembleOnFunctorInterface
  {
  public:
    AssembleOnDirichlet(const GridBoundaryType& boundaryInfo)
      : boundaryInfo_(boundaryInfo)
    {}

    virtual bool assembleOn(const GridPartType& /*gridPart*/, const IntersectionType& intersection) const DS_OVERRIDE
    {
      return boundaryInfo_.dirichlet(intersection);
    }

  private:
    const GridBoundaryType& boundaryInfo_;
  }; // class AssembleOnDirichlet

  class AssembleOnNeumann
    : public AssembleOnFunctorInterface
  {
  public:
    AssembleOnNeumann(const GridBoundaryType& boundaryInfo)
      : boundaryInfo_(boundaryInfo)
    {}

    virtual bool assembleOn(const GridPartType& /*gridPart*/, const IntersectionType& intersection) const DS_OVERRIDE
    {
      return boundaryInfo_.neumann(intersection);
    }

  private:
    const GridBoundaryType& boundaryInfo_;
  }; // class AssembleOnNeumann

private:
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

    virtual void apply(const GridPartType& /*gridPart*/,
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

    virtual void apply(const GridPartType& gridPart,
                       const TestSpaceType& testSpace,
                       const AnsatzSpaceType& ansatzSpace,
                       const IntersectionType& intersection,
                       LocalMatricesContainerType& localMatricesContainer,
                       IndicesContainer& indicesContainer) const DS_OVERRIDE
    {
      if (functor_.assembleOn(gridPart, intersection))
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

    virtual void apply(const GridPartType& /*gridPart*/,
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

    virtual void apply(const GridPartType& gridPart,
                       const TestSpaceType& testSpace,
                       const IntersectionType& intersection,
                       LocalVectorsContainerType& localVectorsContainer,
                       Dune::DynamicVector< size_t >& indices) const DS_OVERRIDE
    {
      if (functor_.assembleOn(gridPart, intersection))
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

public:
  SystemAssembler(const TestSpaceType& test, const AnsatzSpaceType& ansatz, const GridPartType& grid_part)
    : testSpace_(test)
    , ansatzSpace_(ansatz)
    , grid_part_(grid_part)
  {}

  SystemAssembler(const TestSpaceType& test, const AnsatzSpaceType& ansatz)
    : testSpace_(test)
    , ansatzSpace_(ansatz)
    , grid_part_(*(testSpace_.gridPart()))
  {}

  SystemAssembler(const TestSpaceType& test)
    : testSpace_(test)
    , ansatzSpace_(testSpace_)
    , grid_part_(*(testSpace_.gridPart()))
  {}

  SystemAssembler(const TestSpaceType& test, const GridPartType& grid_part)
    : testSpace_(test)
    , ansatzSpace_(testSpace_)
    , grid_part_(grid_part)
  {}

  ~SystemAssembler()
  {
    clearAll();
  }

  void clearLocalAssemblers()
  {
    for (auto& element : localCodim0MatrixAssemblers_)  delete element;
    for (auto& element : localCodim0VectorAssemblers_)  delete element;
    for (auto& element : localCodim1MatrixAssemblers_)  delete element;
    for (auto& element : localCodim1VectorAssemblers_)  delete element;
  }

  void clearLocalConstraints()
  {
    for (auto& element : localConstraints_) delete element;
  }

  void clearAll()
  {
    clearLocalAssemblers();
    clearLocalConstraints();
  }

  const GridPartType& gridPart() const
  {
    return grid_part_;
  }

  const TestSpaceType& testSpace() const
  {
    return testSpace_;
  }

  const AnsatzSpaceType& ansatzSpace() const
  {
    return ansatzSpace_;
  }

  template< class L, class M >
  void addLocalAssembler(const LocalAssembler::Codim0Matrix< L >& localAssembler,
                         Dune::Stuff::LA::MatrixInterface< M >& matrix)
  {
    assert(matrix.rows() == testSpace_.mapper().size());
    assert(matrix.cols() == ansatzSpace_.mapper().size());
    localCodim0MatrixAssemblers_.push_back(new LocalCodim0MatrixAssemblerWrapper< L, M >(localAssembler, matrix));
  }

  template< class L, class AssembleOnFunctorType, class M >
  void addLocalAssembler(const LocalAssembler::Codim1CouplingMatrix< L >& localAssembler,
                         const AssembleOnFunctorType functor,
                         Dune::Stuff::LA::MatrixInterface< M >& matrix)
  {
    dune_static_assert((Dune::IsBaseOf< AssembleOnFunctorInterface, AssembleOnFunctorType >::value),
                       "ERROR: AssembleOnFunctorType is not a AssembleOnFunctorInterface");
    assert(matrix.rows() == testSpace_.mapper().size());
    assert(matrix.cols() == ansatzSpace_.mapper().size());
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
    assert(matrix.rows() == testSpace_.mapper().size());
    assert(matrix.cols() == ansatzSpace_.mapper().size());
    localCodim1MatrixAssemblers_.push_back(
          new LocalCodim1MatrixAssemblerWrapper<  LocalAssembler::Codim1BoundaryMatrix< L >,
                                                  AssembleOnFunctorType,
                                                  M >(localAssembler, functor, matrix));
  }

  template< class L, class V >
  void addLocalAssembler(const LocalAssembler::Codim0Vector< L >& localAssembler,
                         Dune::Stuff::LA::VectorInterface< V >& vector)
  {
    assert(vector.size() == testSpace_.mapper().size());
    localCodim0VectorAssemblers_.push_back(new LocalCodim0VectorAssemblerWrapper< L, V >(localAssembler, vector));
  }

  template< class L, class AssembleOnFunctorType, class V >
  void addLocalAssembler(const LocalAssembler::Codim1Vector< L >& localAssembler,
                         const AssembleOnFunctorType functor,
                         Dune::Stuff::LA::VectorInterface< V >& vector)
  {
    dune_static_assert((Dune::IsBaseOf< AssembleOnFunctorInterface, AssembleOnFunctorType >::value),
                       "ERROR: AssembleOnFunctorType is not a AssembleOnFunctorInterface");
    assert(vector.size() == testSpace_.mapper().size());
    localCodim1VectorAssemblers_.push_back(
          new LocalCodim1VectorAssemblerWrapper<  LocalAssembler::Codim1Vector< L >,
                                                  AssembleOnFunctorType,
                                                  V >(localAssembler, functor, vector));
  }

  void assemble() const
  {
    // only do something, if there are local assemblers
    const bool codim1_assemblers_present =
        (localCodim1MatrixAssemblers_.size() + localCodim1VectorAssemblers_.size()) > 0;
    if ((localCodim0MatrixAssemblers_.size() + localCodim0VectorAssemblers_.size()) > 0
        || codim1_assemblers_present) {
      // common tmp storage for all entities
      // * for the matrix assemblers
      std::vector< size_t > numberOfTmpMatricesNeeded(2, 0);
      for (auto& localCodim0MatrixAssembler : localCodim0MatrixAssemblers_) {
        const auto tmp = localCodim0MatrixAssembler->numTmpObjectsRequired();
        assert(tmp.size() == 2);
        numberOfTmpMatricesNeeded[0] = std::max(numberOfTmpMatricesNeeded[0], tmp[0]);
        numberOfTmpMatricesNeeded[1] = std::max(numberOfTmpMatricesNeeded[1], tmp[1]);
      }
      for (auto& localCodim1MatrixAssembler : localCodim1MatrixAssemblers_) {
        const auto tmp = localCodim1MatrixAssembler->numTmpObjectsRequired();
        assert(tmp.size() == 2);
        numberOfTmpMatricesNeeded[0] = std::max(numberOfTmpMatricesNeeded[0], tmp[0]);
        numberOfTmpMatricesNeeded[1] = std::max(numberOfTmpMatricesNeeded[1], tmp[1]);
      }
      const size_t maxLocalSize = std::max(testSpace_.mapper().maxNumDofs(), ansatzSpace_.mapper().maxNumDofs());
      std::vector< LocalMatrixType > tmpLocalAssemblerMatrices( numberOfTmpMatricesNeeded[0],
                                                                LocalMatrixType(maxLocalSize,
                                                                                maxLocalSize,
                                                                                RangeFieldType(0)));
      std::vector< LocalMatrixType > tmpLocalOperatorMatrices(numberOfTmpMatricesNeeded[1],
                                                              LocalMatrixType(maxLocalSize,
                                                                              maxLocalSize,
                                                                              RangeFieldType(0)));
      std::vector< std::vector< LocalMatrixType > > tmpLocalMatricesContainer;
      tmpLocalMatricesContainer.push_back(tmpLocalAssemblerMatrices);
      tmpLocalMatricesContainer.push_back(tmpLocalOperatorMatrices);
      // * for the vector assemblers
      std::vector< size_t > numberOfTmpVectorsNeeded(2, 0);
      for (auto& localCodim0VectorAssembler : localCodim0VectorAssemblers_) {
        const auto tmp = localCodim0VectorAssembler->numTmpObjectsRequired();
        assert(tmp.size() == 2);
        numberOfTmpVectorsNeeded[0] = std::max(numberOfTmpVectorsNeeded[0], tmp[0]);
        numberOfTmpVectorsNeeded[1] = std::max(numberOfTmpVectorsNeeded[1], tmp[1]);
      }
      for (auto& localCodim1VectorAssembler : localCodim1VectorAssemblers_) {
        const auto tmp = localCodim1VectorAssembler->numTmpObjectsRequired();
        assert(tmp.size() == 2);
        numberOfTmpVectorsNeeded[0] = std::max(numberOfTmpVectorsNeeded[0], tmp[0]);
        numberOfTmpVectorsNeeded[1] = std::max(numberOfTmpVectorsNeeded[1], tmp[1]);
      }
      std::vector< LocalVectorType > tmpLocalAssemblerVectors(numberOfTmpVectorsNeeded[0],
                                                              LocalVectorType(maxLocalSize,
                                                                              RangeFieldType(0)));
      std::vector< LocalVectorType > tmpLocalFunctionalVectors( numberOfTmpVectorsNeeded[1],
                                                                LocalVectorType(maxLocalSize,
                                                                                RangeFieldType(0)));
      std::vector< std::vector< LocalVectorType > > tmpLocalVectorsContainer;
      tmpLocalVectorsContainer.push_back(tmpLocalAssemblerVectors);
      tmpLocalVectorsContainer.push_back(tmpLocalFunctionalVectors);
      // * for the global indices
      std::vector< Dune::DynamicVector< size_t > > tmpIndices = {
          Dune::DynamicVector< size_t >(maxLocalSize)
        , Dune::DynamicVector< size_t >(maxLocalSize)
        , Dune::DynamicVector< size_t >(maxLocalSize)
        , Dune::DynamicVector< size_t >(maxLocalSize)
      };

      // walk the grid
      const auto entityEndIt = grid_part_.template end< 0 >();
      for(auto entityIt = grid_part_.template begin< 0 >(); entityIt != entityEndIt; ++entityIt ) {
        const EntityType& entity = *entityIt;
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.startTiming("GDT.SystemAssembler.assemble.local_codim0_matrix_assemblers");
#endif
        // call local matrix assemblers
        for (auto& localCodim0MatrixAssembler : localCodim0MatrixAssemblers_) {
          localCodim0MatrixAssembler->apply(testSpace_, ansatzSpace_, entity, tmpLocalMatricesContainer, tmpIndices);
        }
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.stopTiming("GDT.SystemAssembler.assemble.local_codim0_matrix_assemblers");
#endif
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.startTiming("GDT.SystemAssembler.assemble.local_codim0_vector_assemblers");
#endif
        // call local vector assemblers
        for (auto& localCodim0VectorAssembler : localCodim0VectorAssemblers_) {
          localCodim0VectorAssembler->apply(testSpace_, entity, tmpLocalVectorsContainer, tmpIndices[0]);
        }
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.stopTiming("GDT.SystemAssembler.assemble.local_codim0_vector_assemblers");
#endif
        // only walk the intersections, if there are local assemblers present
        if (codim1_assemblers_present) {
          // walk the intersections
          const auto intersectionEndIt = grid_part_.iend(entity);
          for (auto intersectionIt = grid_part_.ibegin(entity);
               intersectionIt != intersectionEndIt;
               ++intersectionIt) {
            const auto& intersection = *intersectionIt;
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.startTiming("GDT.SystemAssembler.assemble.local_codim1_matrix_assemblers");
#endif
            // call local matrix assemblers
            for (auto& localCodim1MatrixAssembler : localCodim1MatrixAssemblers_) {
              localCodim1MatrixAssembler->apply(grid_part_, testSpace_, ansatzSpace_,
                                                intersection,
                                                tmpLocalMatricesContainer, tmpIndices);
            }
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.stopTiming("GDT.SystemAssembler.assemble.local_codim1_matrix_assemblers");
#endif
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.startTiming("GDT.SystemAssembler.assemble.local_codim1_vector_assemblers");
#endif
            // call local vector assemblers
            for (auto& localCodim1VectorAssembler : localCodim1VectorAssemblers_) {
              localCodim1VectorAssembler->apply(grid_part_, testSpace_, intersection, tmpLocalVectorsContainer, tmpIndices[0]);
            }
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.stopTiming("GDT.SystemAssembler.assemble.local_codim1_vector_assemblers");
#endif
          } // walk the intersections
        } // only walk the intersections, if there are local assemblers present
      } // walk the grid
    } // only do something, if there are local assemblers
  } // void assemble() const

private:
  class LocalConstraintsApplication
  {
  public:
    virtual ~LocalConstraintsApplication() {}
    virtual void apply(const TestSpaceType& testSpace, const EntityType& entity) = 0;
  };

  template< class ConstraintsType, class M >
  class LocalMatrixConstraintsApplicationWrapper
    : public LocalConstraintsApplication
  {
  public:
    LocalMatrixConstraintsApplicationWrapper(ConstraintsType& constraints,
                                             Dune::Stuff::LA::MatrixInterface< M >& matrix)
      : constraints_(constraints)
      , matrix_(matrix)
    {}

    virtual void apply(const TestSpaceType& testSpace, const EntityType& entity)
    {
      testSpace.localConstraints(entity, constraints_);
      for (size_t ii = 0; ii < constraints_.rows(); ++ii) {
        const size_t row = constraints_.globalRow(ii);
        for (size_t jj = 0; jj < constraints_.cols(); ++jj) {
          matrix_.set_entry(row, constraints_.globalCol(jj), constraints_.value(ii, jj));
        }
      }
    }

  private:
    ConstraintsType& constraints_;
    Dune::Stuff::LA::MatrixInterface< M >& matrix_;
  }; // class LocalMatrixConstraintsApplicationWrapper

  template< class ConstraintsType, class V >
  class LocalVectorConstraintsApplicationWrapper
    : public LocalConstraintsApplication
  {
  public:
    LocalVectorConstraintsApplicationWrapper(ConstraintsType& constraints,
                                             Dune::Stuff::LA::VectorInterface< V >& vector)
      : constraints_(constraints)
      , vector_(vector)
    {}

    virtual void apply(const TestSpaceType& testSpace, const EntityType& entity)
    {
      testSpace.localConstraints(entity, constraints_);
      for (size_t ii = 0; ii < constraints_.rows(); ++ii) {
        vector_.set_entry(constraints_.globalRow(ii), RangeFieldType(0));
      }
    }

  private:
    ConstraintsType& constraints_;
    Dune::Stuff::LA::VectorInterface< V >& vector_;
  }; // class LocalVectorConstraintsApplicationWrapper

public:
  template< class ConstraintsType, class M >
  void addLocalConstraints(ConstraintsType& constraints,
                           Dune::Stuff::LA::MatrixInterface< M >& matrix)
  {
    assert(matrix.rows() == testSpace_.mapper().size());
    assert(matrix.cols() == ansatzSpace_.mapper().size());
    localConstraints_.push_back(new LocalMatrixConstraintsApplicationWrapper< ConstraintsType, M >(constraints, matrix));
  }

  template< class ConstraintsType, class V >
  void addLocalConstraints(ConstraintsType& constraints,
                           Dune::Stuff::LA::VectorInterface< V >& vector)
  {
    assert(vector.size() == testSpace_.mapper().size());
    localConstraints_.push_back(new LocalVectorConstraintsApplicationWrapper< ConstraintsType, V >(constraints, vector));
  }

  void applyConstraints() const
  {
    // walk the grid
    const auto entityEndIt = grid_part_.template end< 0 >();
    for (auto entityIt = grid_part_.template begin< 0 >(); entityIt != entityEndIt; ++entityIt) {
      const EntityType& entity = *entityIt;
      for (auto& localConstraints : localConstraints_)
        localConstraints->apply(testSpace_, entity);
    } // walk the grid
  } // ... applyConstraints()

private:
  const TestSpaceType& testSpace_;
  const AnsatzSpaceType& ansatzSpace_;
  const GridPartType& grid_part_;
  std::vector< LocalCodim0MatrixAssemblerApplication* > localCodim0MatrixAssemblers_;
  std::vector< LocalCodim0VectorAssemblerApplication* > localCodim0VectorAssemblers_;
  std::vector< LocalCodim1MatrixAssemblerApplication* > localCodim1MatrixAssemblers_;
  std::vector< LocalCodim1VectorAssemblerApplication* > localCodim1VectorAssemblers_;
  std::vector< LocalConstraintsApplication* > localConstraints_;
}; // class SystemAssembler


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_SYSTEM_HH
