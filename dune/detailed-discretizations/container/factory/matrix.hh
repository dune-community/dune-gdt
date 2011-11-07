#ifndef DUNE_DETAILED_DISCRETIZATIONS_CONTAINER_FACTORY_MATRIX_HH
#define DUNE_DETAILED_DISCRETIZATIONS_CONTAINER_FACTORY_MATRIX_HH

// dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/fmatrix.hh>

// dune-istl includes
#include <dune/istl/bcrsmatrix.hh>

// dune-detailed-discretizations includes
#include <dune/detailed-discretizations/container/sparsitypattern.hh>

namespace Dune {

namespace DetailedDiscretizations {

//! Contains container classes for storing.
namespace Container {

namespace Matrix {

template <class ContainerImp>
class Factory
{
};

/**
  \brief      Static factory class for matrix classes of type Dune::BCRSMatrix.

  \tparam     BlockType Type to construct a Dune::BCRSMatrix representing the type for
              a block, normally a Dune::FieldMatrix.

  \attention  The sparsity pattern is always created for codim 1 contributions as well! This should be optimized
  \todo       See \attention
 */
template <class BlockType>
class Factory<Dune::BCRSMatrix<BlockType>>
{
public:
  //! Wrapped container type.
  typedef Dune::BCRSMatrix<BlockType> ContainerType;

  //! Return type for create() method.
  typedef std::auto_ptr<ContainerType> AutoPtrType;

  /** @brief Creates a new BCRSMatrix object and returns an auto_ptr pointing
   * to the allocated object.
   *
   * Matrices have size @f$ H \times H @f$ where @f$H@f$ is the number
   * of degrees of freedom in the discrete function space @f$ { \cal X }_H @f$.
   * The matrices' sparsity pattern is determined by the discrete function
   * space's basefunction overlap.
   *
   * @param dfs The discrete function space @f$ { \cal X }_H @f$.
   */
  template <class AnsatzSpaceType, class TestSpaceType>
  static AutoPtrType create(AnsatzSpaceType& ansatzSpace, TestSpaceType& testSpace)
  {
    return AutoPtrType(createPtr(ansatzSpace, testSpace));
  }

  /** @brief Creates a new BCRSMatrix object and returns a pointer to the
   * allocated object.
   *
   * Matrices have size @f$ H \times H @f$ where @f$H@f$ is the number
   * of degrees of freedom in the discrete function space @f$ { \cal X }_H @f$.
   * The matrices' sparsity pattern is determined by the discrete function
   * space's basefunction overlap.
   *
   * @param dfs The discrete function space @f$ { \cal X }_H @f$.
   */
  template <class AnsatzSpaceType, class TestSpaceType>
  static ContainerType* createPtr(AnsatzSpaceType& ansatzSpace, TestSpaceType& testSpace)
  {
    // some types
    typedef typename AnsatzSpaceType::GridPartType GridPartType;

    typedef typename GridPartType::template Codim<0>::IteratorType EntityIteratorType;

    typedef typename EntityIteratorType::Entity EntityType;

    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;

    typedef typename IntersectionIteratorType::Intersection IntersectionType;

    typedef typename IntersectionType::EntityPointer EntityPointerType;

    const unsigned int ansatzSize = ansatzSpace.map().size();
    const unsigned int testSize   = testSpace.map().size();

    typedef Dune::DetailedDiscretizations::Container::SparsityPattern PatternType;

    PatternType sPattern(ansatzSize);

    ContainerType* matrix = new ContainerType(ansatzSize, testSize, ContainerType::random);

    // compute sparsity pattern
    // \todo precompile this in linear subspace
    // \todo use constraints for sparsity pattern
    const EntityIteratorType lastEntity = ansatzSpace.gridPart().template end<0>();
    for (EntityIteratorType entityIterator = ansatzSpace.gridPart().template begin<0>(); entityIterator != lastEntity;
         ++entityIterator) {
      const EntityType& entity = *entityIterator;
      for (unsigned int i = 0; i < ansatzSpace.baseFunctionSet().local(entity).size(); ++i) {
        unsigned int ii = ansatzSpace.map().toGlobal(entity, i);
        for (unsigned int j = 0; j < testSpace.baseFunctionSet().local(entity).size(); ++j) {
          unsigned int jj = testSpace.map().toGlobal(entity, j);
          sPattern.insert(ii, jj);
        }
      }
      // do loop over all intersections
      const IntersectionIteratorType lastIntersection = ansatzSpace.gridPart().iend(entity);
      for (IntersectionIteratorType intIt = ansatzSpace.gridPart().ibegin(entity); intIt != lastIntersection; ++intIt) {
        const IntersectionType& intersection = *intIt;
        // if inner intersection
        if (intersection.neighbor() && !intersection.boundary()) {
          // get neighbouring entity
          const EntityPointerType neighbourPtr = intersection.outside();
          const EntityType& neighbour = *neighbourPtr;
          for (unsigned int i = 0; i < ansatzSpace.baseFunctionSet().local(entity).size(); ++i) {
            unsigned int ii = ansatzSpace.map().toGlobal(entity, i);
            for (unsigned int j = 0; j < testSpace.baseFunctionSet().local(neighbour).size(); ++j) {
              unsigned int jj = testSpace.map().toGlobal(neighbour, j);
              sPattern.insert(ii, jj);
            }
          }
        } // end if inner intersection
      } // done loop over all intersections
    }

    for (unsigned int i = 0; i < sPattern.size(); ++i) {
      matrix->setrowsize(i, sPattern.countNonZeros(i));
    }
    matrix->endrowsizes();

    for (unsigned int i = 0; i < sPattern.size(); ++i) {
      typedef SparsityPattern::NonZeroColIterator ColIterator;
      ColIterator sit = sPattern.begin(i);
      for (; sit != sPattern.end(i); sit++) {
        matrix->addindex(i, *sit);
      }
    }
    matrix->endindices();

    return matrix;
  } // end method createPtr

}; // end class MatrixFactory<BCRSMatrix<T> >

template <class FieldType, int n, int m = n>
class Defaults
{
public:
  typedef Factory<Dune::BCRSMatrix<Dune::FieldMatrix<FieldType, n, n>>> BCRSMatrix;

}; // end class Defaults

} // end namespace Matrix

} // end namespace Container

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_CONTAINER_FACTORY_MATRIX_HH
