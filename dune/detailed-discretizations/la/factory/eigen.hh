#ifndef DUNE_DETAILED_DISCRETIZATIONS_CONTAINER_FACTORY_EIGEN_HH
#define DUNE_DETAILED_DISCRETIZATIONS_CONTAINER_FACTORY_EIGEN_HH

#ifdef HAVE_EIGEN

// dune-common
#include <dune/common/shared_ptr.hh>

// dune-helper-tools
#include <dune/helper-tools/grid/information.hh>

// local
#include "../backend/container/eigen.hh"

namespace Dune {

namespace DetailedDiscretizations {

namespace LA {

namespace Factory {

template <class EntryType>
class Eigen
{
public:
  typedef Dune::DetailedDiscretizations::LA::Backend::Container::Eigen::SparseMatrix<EntryType> SparseMatrixType;

  typedef Dune::DetailedDiscretizations::LA::Backend::Container::Eigen::DenseMatrix<EntryType> DenseMatrixType;

  typedef Dune::DetailedDiscretizations::LA::Backend::Container::Eigen::DenseVector<EntryType> DenseVectorType;

  template <class AnsatzSpaceType, class TestSpaceType>
  static SparseMatrixType createSparseMatrix(const AnsatzSpaceType& ansatzSpace, const TestSpaceType& testSpace)
  {
    // init
    SparseMatrixType matrix(ansatzSpace.map().size(), testSpace.map().size());
    // compute max number of dofs per entity
    const unsigned int maxNumberLocalDofs = std::max(ansatzSpace.map().maxLocalSize(), testSpace.map().maxLocalSize());
    // compute max number of neighbours
    const unsigned int maxNumberNeighbors =
        Dune::HelperTools::Grid::Information::maxNumberOfNeighbors(ansatzSpace.gridPart());
    // reserve
    matrix.reserve((maxNumberNeighbors + 1) * maxNumberLocalDofs);
    // return
    return matrix;
  }

  template <class AnsatzSpaceType, class TestSpaceType>
  static DenseMatrixType createDenseMatrix(const AnsatzSpaceType& ansatzSpace, const TestSpaceType& testSpace)
  {
    // init
    DenseMatrixType matrix(ansatzSpace.map().size(), testSpace.map().size());
    // reserve
    matrix.reserve();
    // return
    return matrix;
  }

  template <class SpaceType>
  static DenseVectorType createDenseVector(const SpaceType& space)
  {
    // init
    DenseVectorType vector(space.map().size());
    // reserve
    vector.reserve();
    // return
    return vector;
  }
}; // class Eigen

} // namespace Factory

} // namespace LA

} // namespace DetailedDiscretizations

} // namespace Dune

#endif // HAVE_EIGEN

#endif // DUNE_DETAILED_DISCRETIZATIONS_CONTAINER_FACTORY_EIGEN_HH
