#ifndef DUNE_DETAILED_DISCRETIZATIONS_CONTAINER_FACTORY_EIGEN_HH
#define DUNE_DETAILED_DISCRETIZATIONS_CONTAINER_FACTORY_EIGEN_HH

#ifdef HAVE_EIGEN

// system
#include <vector>
#include <set>

// dune-common
#include <dune/common/shared_ptr.hh>

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
    std::vector<std::set<unsigned int>> pattern(ansatzSpace.map().size());

    // generate sparsity pattern
    typedef typename AnsatzSpaceType::GridPartType::template Codim<0>::IteratorType IteratorType;
    IteratorType itEnd = ansatzSpace.gridPart().template end<0>();
    for (IteratorType it = ansatzSpace.gridPart().template begin<0>(); it != itEnd; ++it) {
      typename AnsatzSpaceType::GridPartType::template Codim<0>::IteratorType::Entity& entity = *it;
      for (unsigned int i = 0; i < ansatzSpace.baseFunctionSet().local(entity).size(); ++i) {
        for (unsigned int j = 0; j < testSpace.baseFunctionSet().local(entity).size(); ++j) {
          const unsigned int globalI = ansatzSpace.map().toGlobal(entity, i);
          const unsigned int globalJ = testSpace.map().toGlobal(entity, j);
          pattern[globalI].insert(globalJ);
        }
      }
    } // generate sparsity pattern

    // tell pattern to matrix
    for (unsigned int row = 0; row < ansatzSpace.map().size(); ++row) {
      matrix.storage()->startVec(row);
      for (typename std::set<unsigned int>::iterator it = pattern[row].begin(); it != pattern[row].end(); ++it) {
        unsigned int column = *it;
        matrix.storage()->insertBackByOuterInner(row, column);
      }
    } // tell pattern to matrix

    // finalize matrix
    matrix.storage()->finalize();
    matrix.storage()->makeCompressed();

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
