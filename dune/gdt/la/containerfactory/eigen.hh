// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LA_CONTAINER_FACTORY_EIGEN_HH
#define DUNE_GDT_LA_CONTAINER_FACTORY_EIGEN_HH

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
