#ifndef DUNE_DETAILED_DISCRETIZATIONS_LA_CONTAINER_FACTORY_EIGEN_HH
#define DUNE_DETAILED_DISCRETIZATIONS_LA_CONTAINER_FACTORY_EIGEN_HH

//#ifdef HAVE_EIGEN

#include <dune/common/shared_ptr.hh>

#include <dune/stuff/la/container/pattern.hh>
#include <dune/stuff/la/container/eigen.hh>

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace LA {
namespace Container {
namespace Factory {

template< class ElementType >
class Eigen
{
public:
  typedef Dune::Stuff::LA::Container::Eigen::SparseMatrix< ElementType > SparseMatrixType;

  typedef Dune::Stuff::LA::Container::Eigen::DenseMatrix< ElementType > DenseMatrixType;

  typedef Dune::Stuff::LA::Container::Eigen::DenseVector< ElementType > DenseVectorType;

  template< class TestSpaceType, class AnsatzSpaceType >
  static Dune::shared_ptr< SparseMatrixType > createSparseMatrix(const TestSpaceType& testSpace,
                                                                 const AnsatzSpaceType& ansatzSpace)
  {
    typedef Dune::Stuff::LA::Container::Pattern::Default PatternType;

    const Dune::shared_ptr< const PatternType > pattern = testSpace.computePattern(ansatzSpace);
    Dune::shared_ptr< SparseMatrixType > sparseMatrix(new SparseMatrixType(testSpace.map().size(),
                                                                           ansatzSpace.map().size(),
                                                                           *pattern));
    return sparseMatrix;
  } // static ... createSparseMatrix(...)

  template< class TestSpaceType, class AnsatzSpaceType >
  static Dune::shared_ptr< DenseMatrixType > createDenseMatrix(const TestSpaceType& testSpace,
                                                               const AnsatzSpaceType& ansatzSpace)
  {
    Dune::shared_ptr< DenseMatrixType > denseMatrix(new DenseMatrixType(testSpace.map().size(),
                                                                        ansatzSpace.map().size()));
    return denseMatrix;
  } // static ... createDenseMatrix(...)

  template< class SpaceType >
  static Dune::shared_ptr< DenseVectorType > createDenseVector(const SpaceType& space)
  {
    Dune::shared_ptr< DenseVectorType > denseVector(new DenseVectorType(space.map().size()));
    return denseVector;
  } // static Dune::shared_ptr< DenseVectorType > createDenseVector(const SpaceType& space)
}; // class Eigen

} // namespace Factory
} // namespace Conatiner
} // namespace LA
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

//#endif // HAVE_EIGEN

#endif // DUNE_DETAILED_DISCRETIZATIONS_LA_CONTAINER_FACTORY_EIGEN_HH
