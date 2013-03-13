#ifndef DUNE_DETAILED_DISCRETIZATIONS_LA_CONTAINER_FACTORY_EIGEN_HH
#define DUNE_DETAILED_DISCRETIZATIONS_LA_CONTAINER_FACTORY_EIGEN_HH

#ifdef HAVE_EIGEN

#include <dune/common/shared_ptr.hh>

#include <dune/stuff/la/container/pattern.hh>
#include <dune/stuff/la/container/eigen.hh>

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace LA {
namespace Container {
namespace Factory {

template <class ElementType>
class Eigen
{
public:
  typedef Dune::Stuff::LA::Container::EigenRowMajorSparseMatrix<ElementType> RowMajorSparseMatrixType;
  typedef Dune::Stuff::LA::Container::EigenDenseMatrix<ElementType> DenseMatrixType;
  typedef Dune::Stuff::LA::Container::EigenDenseVector<ElementType> DenseVectorType;
  typedef Dune::Stuff::LA::Container::SparsityPatternDefault PatternType;

  template <class TestSpaceType, class AnsatzSpaceType>
  static Dune::shared_ptr<RowMajorSparseMatrixType> createRowMajorSparseMatrix(const TestSpaceType& testSpace,
                                                                               const AnsatzSpaceType& ansatzSpace)
  {
    typedef Dune::Stuff::LA::Container::SparsityPatternDefault PatternType;
    const Dune::shared_ptr<const PatternType> pattern = testSpace.computePattern(ansatzSpace);
    return createRowMajorSparseMatrix(testSpace, ansatzSpace, *pattern);
  } // static ... createRowMajorSparseMatrix(...)

  template <class TestSpaceType, class AnsatzSpaceType>
  static Dune::shared_ptr<RowMajorSparseMatrixType> createRowMajorSparseMatrix(const TestSpaceType& testSpace,
                                                                               const AnsatzSpaceType& ansatzSpace,
                                                                               const PatternType& pattern)
  {
    return Dune::make_shared<RowMajorSparseMatrixType>(testSpace.map().size(), ansatzSpace.map().size(), pattern);
  } // static ... createRowMajorSparseMatrix(...)

  template <class TestSpaceType, class AnsatzSpaceType>
  static Dune::shared_ptr<DenseMatrixType> createDenseMatrix(const TestSpaceType& testSpace,
                                                             const AnsatzSpaceType& ansatzSpace)
  {
    return Dune::make_shared<DenseMatrixType>(testSpace.map().size(), ansatzSpace.map().size());
  } // static ... createDenseMatrix(...)

  template <class SpaceType>
  static Dune::shared_ptr<DenseVectorType> createDenseVector(const SpaceType& space)
  {
    return Dune::make_shared<DenseVectorType>(space.map().size());
  } // static Dune::shared_ptr< DenseVectorType > createDenseVector(const SpaceType& space)
}; // class Eigen

} // namespace Factory
} // namespace Conatiner
} // namespace LA
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // HAVE_EIGEN

#endif // DUNE_DETAILED_DISCRETIZATIONS_LA_CONTAINER_FACTORY_EIGEN_HH
