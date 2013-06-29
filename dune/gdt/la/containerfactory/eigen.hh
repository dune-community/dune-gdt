#ifndef DUNE_GDT_LA_CONTAINER_FACTORY_EIGEN_HH
#define DUNE_GDT_LA_CONTAINER_FACTORY_EIGEN_HH

// we need this for the HAVE_EIGEN define!
#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif

#include <dune/stuff/la/container/pattern.hh>
#include <dune/stuff/la/container/eigen.hh>

#include <dune/gdt/space/interface.hh>

namespace Dune {
namespace GDT {


template< class ElementType >
class ContainerFactoryEigen
{
public:
  typedef Dune::Stuff::LA::EigenRowMajorSparseMatrix< ElementType >  RowMajorSparseMatrixType;
  typedef Dune::Stuff::LA::EigenDenseMatrix< ElementType >           DenseMatrixType;
  typedef Dune::Stuff::LA::EigenDenseVector< ElementType >           DenseVectorType;
  typedef Dune::Stuff::LA::SparsityPatternDefault                    PatternType;

  template< class T, class A >
  static RowMajorSparseMatrixType* createRowMajorSparseMatrix(const SpaceInterface< T >& testSpace,
                                                              const SpaceInterface< A >& ansatzSpace)
  {
    const std::shared_ptr< const PatternType > pattern(testSpace.computePattern(ansatzSpace));
    return createRowMajorSparseMatrix(testSpace, ansatzSpace, *pattern);
  } // static ... createRowMajorSparseMatrix(...)

  template< class T, class A >
  static RowMajorSparseMatrixType* createRowMajorSparseMatrix(const SpaceInterface< T >& testSpace,
                                                              const SpaceInterface< A >& ansatzSpace,
                                                              const PatternType& pattern)
  {
    return new RowMajorSparseMatrixType(testSpace.mapper().size(),
                                        ansatzSpace.mapper().size(),
                                        pattern);
  } // static ... createRowMajorSparseMatrix(...)

  template< class T, class A >
  static DenseMatrixType *createDenseMatrix(const SpaceInterface< T >& testSpace,
                                            const SpaceInterface< A >& ansatzSpace)
  {
    return new DenseMatrixType(testSpace.mapper().size(),
                               ansatzSpace.mapper().size());
  } // static ... createDenseMatrix(...)

  template< class S >
  static DenseVectorType* createDenseVector(const SpaceInterface< S >& space)
  {
    return new DenseVectorType(space.mapper().size());
  } // static ... createDenseVector(const SpaceType& space)
}; // class ContainerFactoryEigen


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LA_CONTAINER_FACTORY_EIGEN_HH
