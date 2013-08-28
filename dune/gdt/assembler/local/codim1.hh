// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_ASSEMLBER_LOCAL_CODIM1_HH
#define DUNE_GDT_ASSEMLBER_LOCAL_CODIM1_HH

#include <vector>

#include <dune/stuff/common/matrix.hh>
#include <dune/stuff/la/container/interface.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#ifdef DUNE_STUFF_PROFILER_ENABLED
#include <dune/stuff/common/profiler.hh>
#endif

#include <dune/gdt/localoperator/interface.hh>
#include <dune/gdt/localfunctional/interface.hh>
#include <dune/gdt/space/interface.hh>

namespace Dune {
namespace GDT {
namespace LocalAssembler {


template< class LocalOperatorImp >
class Codim1CouplingMatrix;


template< class LocalOperatorImp >
class Codim1CouplingMatrixTraits
{
public:
  typedef Codim1CouplingMatrix< LocalOperatorImp > derived_type;
  typedef LocalOperator::Codim1CouplingInterface< typename LocalOperatorImp::Traits > LocalOperatorType;
};


template< class LocalOperatorImp >
class Codim1CouplingMatrix
{
public:
  typedef Codim1CouplingMatrixTraits< LocalOperatorImp > Traits;
  typedef typename Traits::LocalOperatorType LocalOperatorType;

  Codim1CouplingMatrix(const LocalOperatorType& op)
    : localOperator_(op)
  {}

  const LocalOperatorType& localOperator() const
  {
    return localOperator_;
  }

private:
  static const size_t numTmpObjectsRequired_ = 4;

public:
  std::vector< size_t > numTmpObjectsRequired() const
  {
    return {numTmpObjectsRequired_, localOperator_.numTmpObjectsRequired()};
  }

  template< class TE, class AE, class TN, class AN,
            class IntersectionType, class MEE, class MNN, class MEN, class MNE, class R >
  void assembleLocal(const SpaceInterface< TE >& testSpaceEntity,
                     const SpaceInterface< AE >& ansatzSpaceEntity,
                     const SpaceInterface< TN >& testSpaceNeighbor,
                     const SpaceInterface< AN >& ansatzSpaceNeighbor,
                     const IntersectionType& intersection,
                     Dune::Stuff::LA::MatrixInterface< MEE >& entityEntityMatrix,
                     Dune::Stuff::LA::MatrixInterface< MNN >& neighborNeighborMatrix,
                     Dune::Stuff::LA::MatrixInterface< MEN >& entityNeighborMatrix,
                     Dune::Stuff::LA::MatrixInterface< MNE >& neighborEntityMatrix,
                     std::vector< std::vector< Dune::DynamicMatrix< R > > >& tmpLocalMatricesContainer,
                     std::vector< Dune::DynamicVector< size_t > >& tmpIndicesContainer) const
  {
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.startTiming("GDT.LocalAssembler.Codim1CouplingMatrix.assembleLocal");
#endif
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.startTiming("GDT.LocalAssembler.Codim1CouplingMatrix.assembleLocal.1_check_and_clear");
#endif
    // check
    assert(tmpLocalMatricesContainer.size() >= 2);
    assert(tmpLocalMatricesContainer[0].size() >= numTmpObjectsRequired_);
    assert(tmpLocalMatricesContainer[1].size() >= localOperator_.numTmpObjectsRequired());
    assert(tmpIndicesContainer.size() >= 4);
    // get and clear matrix
    Dune::DynamicMatrix< R >& localEntityEntityMatrix = tmpLocalMatricesContainer[0][0];
    Dune::DynamicMatrix< R >& localNeighborNeighborMatrix = tmpLocalMatricesContainer[0][1];
    Dune::DynamicMatrix< R >& localEntityNeighborMatrix = tmpLocalMatricesContainer[0][2];
    Dune::DynamicMatrix< R >& localNeighborEntityMatrix = tmpLocalMatricesContainer[0][3];
    Dune::Stuff::Common::clear(localEntityEntityMatrix);
    Dune::Stuff::Common::clear(localNeighborNeighborMatrix);
    Dune::Stuff::Common::clear(localEntityNeighborMatrix);
    Dune::Stuff::Common::clear(localNeighborEntityMatrix);
    auto& tmpOperatorMatrices = tmpLocalMatricesContainer[1];
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.stopTiming("GDT.LocalAssembler.Codim1CouplingMatrix.assembleLocal.1_check_and_clear");
#endif
    // get entities
    const auto entityPtr = intersection.inside();
    const auto& entity = *entityPtr;
    const auto neighborPtr = intersection.outside();
    const auto& neighbor = *neighborPtr;
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.startTiming("GDT.LocalAssembler.Codim1CouplingMatrix.assembleLocal.2_apply_local_operator");
#endif
    // apply local operator (results are in local*Matrix)
    localOperator_.apply(testSpaceEntity.baseFunctionSet(entity), ansatzSpaceEntity.baseFunctionSet(entity),
                         testSpaceNeighbor.baseFunctionSet(neighbor), ansatzSpaceNeighbor.baseFunctionSet(neighbor),
                         intersection,
                         localEntityEntityMatrix,
                         localNeighborNeighborMatrix,
                         localEntityNeighborMatrix,
                         localNeighborEntityMatrix,
                         tmpOperatorMatrices);
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.stopTiming("GDT.LocalAssembler.Codim1CouplingMatrix.assembleLocal.2_apply_local_operator");
#endif
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.startTiming("GDT.LocalAssembler.Codim1CouplingMatrix.assembleLocal.3_map_indices");
#endif
    // write local matrices to global
    const size_t rowsEn = testSpaceEntity.mapper().numDofs(entity);
    const size_t colsEn = ansatzSpaceEntity.mapper().numDofs(entity);
    const size_t rowsNe = testSpaceNeighbor.mapper().numDofs(neighbor);
    const size_t colsNe = ansatzSpaceNeighbor.mapper().numDofs(neighbor);
    Dune::DynamicVector< size_t >& globalRowsEn = tmpIndicesContainer[0];
    Dune::DynamicVector< size_t >& globalColsEn = tmpIndicesContainer[1];
    Dune::DynamicVector< size_t >& globalRowsNe = tmpIndicesContainer[2];
    Dune::DynamicVector< size_t >& globalColsNe = tmpIndicesContainer[3];
    assert(globalRowsEn.size() >= rowsEn);
    assert(globalColsEn.size() >= colsEn);
    assert(globalRowsNe.size() >= rowsNe);
    assert(globalColsNe.size() >= colsNe);
    testSpaceEntity.mapper().globalIndices(entity, globalRowsEn);
    ansatzSpaceEntity.mapper().globalIndices(entity, globalColsEn);
    testSpaceNeighbor.mapper().globalIndices(neighbor, globalRowsNe);
    ansatzSpaceNeighbor.mapper().globalIndices(neighbor, globalColsNe);
    assert(localEntityEntityMatrix.rows() >= rowsEn);
    assert(localEntityEntityMatrix.cols() >= colsEn);
    assert(localNeighborNeighborMatrix.rows() >= rowsNe);
    assert(localNeighborNeighborMatrix.cols() >= colsNe);
    assert(localEntityNeighborMatrix.rows() >= rowsEn);
    assert(localEntityNeighborMatrix.cols() >= colsNe);
    assert(localNeighborEntityMatrix.rows() >= rowsNe);
    assert(localNeighborEntityMatrix.cols() >= colsEn);
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.stopTiming("GDT.LocalAssembler.Codim1CouplingMatrix.assembleLocal.3_map_indices");
#endif
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.startTiming("GDT.LocalAssembler.Codim1CouplingMatrix.assembleLocal.4_write_matrices");
#endif
    for (size_t ii = 0; ii < rowsEn; ++ii) {
      const auto& localEntityEntityMatrixRow = localEntityEntityMatrix[ii];
      const auto& localEntityNeighborMatrixRow = localEntityNeighborMatrix[ii];
      const size_t globalII = globalRowsEn[ii];
      for (size_t jj = 0; jj < colsEn; ++jj) {
        const size_t globalJJ = globalColsEn[jj];
        entityEntityMatrix.add(globalII, globalJJ, localEntityEntityMatrixRow[jj]);
      }
      for (size_t jj = 0; jj < colsNe; ++jj) {
        const size_t globalJJ = globalColsNe[jj];
        entityNeighborMatrix.add(globalII, globalJJ, localEntityNeighborMatrixRow[jj]);
      }
    }
    for (size_t ii = 0; ii < rowsNe; ++ii) {
      const auto& localNeighborEntityMatrixRow = localNeighborEntityMatrix[ii];
      const auto& localNeighborNeighborMatrixRow = localNeighborNeighborMatrix[ii];
      const size_t globalII = globalRowsNe[ii];
      for (size_t jj = 0; jj < colsEn; ++jj) {
        const size_t globalJJ = globalColsEn[jj];
        neighborEntityMatrix.add(globalII, globalJJ, localNeighborEntityMatrixRow[jj]);
      }
      for (size_t jj = 0; jj < colsNe; ++jj) {
        const size_t globalJJ = globalColsNe[jj];
        neighborNeighborMatrix.add(globalII, globalJJ, localNeighborNeighborMatrixRow[jj]);
      }
    }
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.stopTiming("GDT.LocalAssembler.Codim1CouplingMatrix.assembleLocal.4_write_matrices");
#endif
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.stopTiming("GDT.LocalAssembler.Codim1CouplingMatrix.assembleLocal");
#endif
  } // void assembleLocal(...) const

  template< class T, class A, class IntersectionType, class M, class R >
  void assembleLocal(const SpaceInterface< T >& testSpace,
                     const SpaceInterface< A >& ansatzSpace,
                     const IntersectionType& intersection,
                     Dune::Stuff::LA::MatrixInterface< M >& systemMatrix,
                     std::vector< std::vector< Dune::DynamicMatrix< R > > >& tmpLocalMatricesContainer,
                     std::vector< Dune::DynamicVector< size_t > >& tmpIndicesContainer) const
  {
    assembleLocal(testSpace, ansatzSpace, testSpace, ansatzSpace,
                  intersection,
                  systemMatrix, systemMatrix, systemMatrix, systemMatrix,
                  tmpLocalMatricesContainer,
                  tmpIndicesContainer);
  } // void assembleLocal(...) const

private:
  const LocalOperatorType& localOperator_;
}; // class Codim1CouplingMatrix


template< class LocalOperatorImp >
class Codim1BoundaryMatrix;


template< class LocalOperatorImp >
class Codim1BoundaryMatrixTraits
{
public:
  typedef Codim1BoundaryMatrix< LocalOperatorImp > derived_type;
  typedef LocalOperator::Codim1BoundaryInterface< typename LocalOperatorImp::Traits > LocalOperatorType;
};


template< class LocalOperatorImp >
class Codim1BoundaryMatrix
{
public:
  typedef Codim1BoundaryMatrixTraits< LocalOperatorImp > Traits;
  typedef typename Traits::LocalOperatorType LocalOperatorType;

  Codim1BoundaryMatrix(const LocalOperatorType& op)
    : localOperator_(op)
  {}

  const LocalOperatorType& localOperator() const
  {
    return localOperator_;
  }

private:
  static const size_t numTmpObjectsRequired_ = 1;

public:
  std::vector< size_t > numTmpObjectsRequired() const
  {
    return {numTmpObjectsRequired_, localOperator_.numTmpObjectsRequired()};
  }

  template< class T, class A, class IntersectionType, class M, class R >
  void assembleLocal(const SpaceInterface< T >& testSpace,
                     const SpaceInterface< A >& ansatzSpace,
                     const IntersectionType& intersection,
                     Dune::Stuff::LA::MatrixInterface< M >& systemMatrix,
                     std::vector< std::vector< Dune::DynamicMatrix< R > > >& tmpLocalMatricesContainer,
                     std::vector< Dune::DynamicVector< size_t > >& tmpIndicesContainer) const
  {
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.startTiming("GDT.LocalAssembler.Codim1BoundaryMatrix.assembleLocal");
#endif
    // check
    assert(tmpLocalMatricesContainer.size() >= 2);
    assert(tmpLocalMatricesContainer[0].size() >= numTmpObjectsRequired_);
    assert(tmpLocalMatricesContainer[1].size() >= localOperator_.numTmpObjectsRequired());
    assert(tmpIndicesContainer.size() >= 2);
    // get and clear matrix
    Dune::DynamicMatrix< R >& localMatrix = tmpLocalMatricesContainer[0][0];
    Dune::Stuff::Common::clear(localMatrix);
    auto& tmpOperatorMatrices = tmpLocalMatricesContainer[1];
    // get entity
    const auto entityPtr = intersection.inside();
    const auto& entity = *entityPtr;
    // apply local operator (results are in local*Matrix)
    localOperator_.apply(testSpace.baseFunctionSet(entity), ansatzSpace.baseFunctionSet(entity),
                         intersection,
                         localMatrix, tmpOperatorMatrices);
    // write local matrices to global
    const size_t rows = testSpace.mapper().numDofs(entity);
    const size_t cols = ansatzSpace.mapper().numDofs(entity);
    Dune::DynamicVector< size_t >& globalRows = tmpIndicesContainer[0];
    Dune::DynamicVector< size_t >& globalCols = tmpIndicesContainer[1];
    assert(globalRows.size() >= rows);
    assert(globalCols.size() >= cols);
    assert(localMatrix.size() >= rows);
    assert(localMatrix.size() >= cols);
    testSpace.mapper().globalIndices(entity, globalRows);
    ansatzSpace.mapper().globalIndices(entity, globalCols);
    for (size_t ii = 0; ii < rows; ++ii) {
      const auto& localMatrixRow = localMatrix[ii];
      const size_t globalII = globalRows[ii];
      for (size_t jj = 0; jj < cols; ++jj) {
        const size_t globalJJ = globalCols[jj];
        systemMatrix.add(globalII, globalJJ, localMatrixRow[jj]);
      }
    }
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.stopTiming("GDT.LocalAssembler.Codim1BoundaryMatrix.assembleLocal");
#endif
  } // void assembleLocal(...) const

private:
  const LocalOperatorType& localOperator_;
}; // class Codim1BoundaryMatrix


template< class LocalFunctionalImp >
class Codim1Vector;


template< class LocalFunctionalImp >
class Codim1VectorTraits
{
public:
  typedef Codim1Vector< LocalFunctionalImp > derived_type;
  typedef LocalFunctional::Codim1Interface< typename LocalFunctionalImp::Traits > LocalFunctionalType;
};


template< class LocalFunctionalImp >
class Codim1Vector
{
public:
  typedef Codim1VectorTraits< LocalFunctionalImp > Traits;
  typedef typename Traits::LocalFunctionalType LocalFunctionalType;

  Codim1Vector(const LocalFunctionalType& fu)
    : localFunctional_(fu)
  {}

  const LocalFunctionalType& localFunctional() const
  {
    return localFunctional_;
  }

private:
  static const size_t numTmpObjectsRequired_ = 1;

public:
  std::vector< size_t > numTmpObjectsRequired() const
  {
    return {numTmpObjectsRequired_, localFunctional_.numTmpObjectsRequired()};
  }

  template< class T, class IntersectionType, class V, class R >
  void assembleLocal(const SpaceInterface< T >& testSpace,
                     const IntersectionType& intersection,
                     Dune::Stuff::LA::VectorInterface< V >& systemVector,
                     std::vector< std::vector< Dune::DynamicVector< R > > >& tmpLocalVectorsContainer,
                     Dune::DynamicVector< size_t >& tmpIndicesContainer) const
  {
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.startTiming("GDT.LocalAssembler.Codim1Vector.assembleLocal");
#endif
    // check
    assert(tmpLocalVectorsContainer.size() >= 2);
    assert(tmpLocalVectorsContainer[0].size() >= numTmpObjectsRequired_);
    assert(tmpLocalVectorsContainer[1].size() >= localFunctional_.numTmpObjectsRequired());
    // get and clear vector
    Dune::DynamicVector< R >& localVector = tmpLocalVectorsContainer[0][0];
    Dune::Stuff::Common::clear(localVector);
    auto& tmpFunctionalVectors = tmpLocalVectorsContainer[1];
    // get entity
    const auto entityPtr = intersection.inside();
    const auto& entity = *entityPtr;
    // apply local functional (results are in localVector)
    localFunctional_.apply(testSpace.baseFunctionSet(entity), intersection, localVector, tmpFunctionalVectors);
    // write local vectors to global
    const size_t size = testSpace.mapper().numDofs(entity);
    assert(tmpIndicesContainer.size() >= size);
    assert(localVector.size() >= size);
    testSpace.mapper().globalIndices(entity, tmpIndicesContainer);
    for (size_t ii = 0; ii < size; ++ii) {
      const size_t globalII = tmpIndicesContainer[ii];
      systemVector.add(globalII, localVector[ii]);
    }
#ifdef DUNE_STUFF_PROFILER_ENABLED
      DSC_PROFILER.stopTiming("GDT.LocalAssembler.Codim1Vector.assembleLocal");
#endif
  } // void assembleLocal(...) const

private:
  const LocalFunctionalType& localFunctional_;
}; // class Codim1Vector


} // namespace LocalAssembler
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMLBER_LOCAL_CODIM1_HH
