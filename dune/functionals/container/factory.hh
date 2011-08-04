#ifndef DUNE_FEM_FUNCTIONALS_CONTAINER_FACTORY_HH
#define DUNE_FEM_FUNCTIONALS_CONTAINER_FACTORY_HH

// dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>

// dune-istl includes
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/bvector.hh>

// local includes
#include "sparsitypattern.hh"

namespace Dune
{

namespace Functionals
{

//! Contains container classes for storing.
namespace Container
{

namespace Matrix
{

///**
// * @brief Interface of static factory class for matrix and vector classes.
// *
// * @tparam ContainerImp Type of the container, for example a vector or matrix container.
// */
//template< class ContainerImp >
//class Factory
//{
//private:
//  class NonImplemented
//  {
//  };

//public:
//  //! Return type for create() method.
//  typedef NonImplemented
//    AutoPtrType;
//  //! Wrapped container type.
//  typedef NonImplemented
//    ContainerType;

//public:
//  /** @brief Creates a new matrix/vector object and returns an auto_ptr
//   * pointing to the allocated object.
//   *
//   * - Matrices have size @f$ H \times H @f$ where @f$H@f$ is the number
//   * of degrees of freedom in the discrete function space @f$ { \cal X }_H @f$.
//   * The matrices' sparsity pattern is determined by the discrete function
//   * space's basefunction overlap.
//   * - Vectors have @f$ H @f$ components
//   *
//   * @param dfs The discrete function space @f$ { \cal X }_H @f$.
//   */
//  template< class DiscFuncSpace >
//  static AutoPtrType create( DiscFuncSpace& dfs )
//  {
//    DUNE_THROW( InvalidStateException, "Factory not implemented for ContainerType!" );
//  }

//  /** @brief Creates a new matrix/vector object and returns a pointer to the
//   * allocated object.
//   *
//   * - Matrices have size @f$ H \times H @f$ where @f$H@f$ is the number
//   * of degrees of freedom in the discrete function space @f$ { \cal X }_H @f$.
//   * The matrices' sparsity pattern is determined by the discrete function
//   * space's basefunction overlap.
//   * - Vectors have @f$ H @f$ components
//   *
//   * @param dfs The discrete function space @f$ { \cal X }_H @f$.
//   */
//  template< class DiscFuncSpace >
//  static ContainerType* createPtr( DiscFuncSpace& dfs )
//  {
//    DUNE_THROW( InvalidStateException, "Factory not implemented for ContainerType!" );
//  }
//}; // end of class Factory

/**
 * @brief Interface of static factory class.
 */
template< class ContainerImp >
class Factory
//  : public Factory< ContainerImp >
{
};

/**
 * @brief Static factory class for matrix classes of type Dune::BCRSMatrix.
 *
 * @tparam T Type to construct a Dune::BCRSMatrix representing the type for
 * a block, normally a Dune::FieldMatrix.
 */
template< class T >
class Factory< Dune::BCRSMatrix< T > >
{
public:

  //! Wrapped container type.
  typedef Dune::BCRSMatrix< T >
    ContainerType;

  //! Return type for create() method.
  typedef std::auto_ptr< ContainerType >
    AutoPtrType;

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
  template< class AnsatzSpaceType, class TestSpaceType >
  static AutoPtrType create( AnsatzSpaceType& ansatzSpace, TestSpaceType& testSpace )
  {
    return AutoPtrType( createPtr( ansatzSpace, testSpace ) );
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
  template< class AnsatzSpaceType, class TestSpaceType >
  static ContainerType* createPtr( AnsatzSpaceType& ansatzSpace, TestSpaceType& testSpace )
  {
    // some types
    typedef typename AnsatzSpaceType::GridPartType
      GridPartType;

    typedef typename GridPartType::template Codim< 0 >::IteratorType
      EntityIteratorType;

    typedef typename EntityIteratorType::Entity
      EntityType;

    const unsigned int ansatzSize = ansatzSpace.map().size();
    const unsigned int testSize = testSpace.map().size();

    typedef Dune::Functionals::Container::SparsityPattern
      PatternType;

    PatternType sPattern( ansatzSize );

    ContainerType* matrix = new ContainerType( ansatzSize,
                                               testSize,
                                               ContainerType::random );

    // compute sparsity pattern
    // \todo precompile this in linear subspace
    // \todo use constraints for sparsity pattern
    const EntityIteratorType lastEntity = ansatzSpace.gridPart().template end< 0 >();
    for(  EntityIteratorType entityIterator = ansatzSpace.gridPart().template begin< 0 >();
          entityIterator != lastEntity;
          ++entityIterator )
    {
      const EntityType& entity = *entityIterator;
      for( unsigned int i = 0; i < ansatzSpace.baseFunctionSet().local( entity ).size(); ++i )
      {
        unsigned int ii = ansatzSpace.map().toGlobal( entity, i );
        for( unsigned int j = 0; j < testSpace.baseFunctionSet().local( entity ).size(); ++j )
        {
          unsigned int jj = testSpace.map().toGlobal( entity, j );
          sPattern.insert( ii, jj );
        }
      }
    }

    for( unsigned int i = 0; i < sPattern.size(); ++i )
    {
      matrix->setrowsize( i, sPattern.countNonZeros( i ) );
    }
    matrix->endrowsizes();

    for( unsigned int i = 0; i < sPattern.size(); ++i )
    {
      typedef SparsityPattern::NonZeroColIterator
        ColIterator;
      ColIterator sit = sPattern.begin( i );
      for( ; sit!=sPattern.end( i ); sit++ )
      {
        matrix->addindex( i, *sit );
      }
    }
    matrix->endindices();

    return matrix;
  } // end method createPtr

}; // end class MatrixFactory<BCRSMatrix<T> >

template< class FieldType, int n, int m = n >
class Defaults
{
public:

  typedef Factory< Dune::BCRSMatrix< Dune::FieldMatrix< FieldType, n, n > > >
    BCRSMatrix;

}; // end class Defaults

} // end namespace Matrix

namespace Vector
{

/**
 * @brief Interface of static factory class.
 */
template< class ContainerImp >
class Factory
//  : public Factory< ContainerImp >
{
};

///**
// * @brief Interface of static factory class for vector classes.
// */
//template< class ContainerImp >
//class VectorFactory
//  : public Factory< ContainerImp >
//{
//};

/**
 * @brief Static factory class for vector classes of type Dune::BlockVector.
 *
 * @tparam T Type to construct a Dune::BlockVector representing the type for
 * a block, normally a Dune::FieldVector.
 */
template< class T >
class Factory< Dune::BlockVector< T > >
{
public:
  //! \copydoc Factory::ContainerType
  typedef Dune::BlockVector< T >
    ContainerType;
  //! \copydoc Factory::AutoPtrType
  typedef Dune::shared_ptr< ContainerType >
    AutoPtrType;

public:
  /** @brief Creates a new vector object and returns an auto_ptr pointing to
   * the allocated object.
   *
   * The vector has @f$ H @f$ components which is the number of degrees of
   * freedom of the given discrete function space @f$ {\cal X}_H @f$.
   *
   * @param dfs The discrete function space @f$ { \cal X }_H @f$.
   */
  template< class DFSType >
  static ContainerType* createPtr( DFSType& dfs )
  {
    const unsigned int numDofs = dfs.map().size();
    ContainerType* bv = new ContainerType( numDofs );
    return bv;
  }

  /** @brief Creates a new Dune::BlockVector object and returns an auto_ptr pointing
   * to the allocated object.
   *
   * Block vectors have size @f$H@f$ where @f$H@f$ is the number
   * of degrees of freedom in the discrete function space @f$ { \cal X }_H @f$.
   *
   * @param dfs The discrete function space @f$ { \cal X }_H @f$.
   */
  template< class DFSType >
  static AutoPtrType create( DFSType& dfs )
  {
    return AutoPtrType( createPtr( dfs ) );
  }
}; // end of VectorFactory<BlockVector<T> >

template< class FieldType, int n >
class Defaults
{
public:

  typedef Factory< Dune::BlockVector< Dune::FieldVector< FieldType, n > > >
    BlockVector;

}; // end class Defaults

} // end namespace Vector

} // end namespace Container

} // end namespace Functionals

} // end namespace Dune


#endif /* end of include guard: DUNE_FEM_FUNCTIONALS_CONTAINER_FACTORY_HH */
