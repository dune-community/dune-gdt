#ifndef DUNE_FEM_FUNCTIONALS_CONTAINER_FACTORY_HH
#define DUNE_FEM_FUNCTIONALS_CONTAINER_FACTORY_HH

#include <dune/common/exceptions.hh>
#include <dune/fem/container/sparsitypattern.hh>

namespace Dune
{
namespace Functionals
{
namespace Container
{

/** @brief interface of static factory class for matrix and vector classes
 */
template< class ContainerImp >
class Factory
{
private:
  class NonImplemented
  {
  };

public:
  //! return type for create() method
  typedef NonImplemented
    AutoPtrType;
  //! wrapped container type
  typedef NonImplemented
    ContainerType;

public:
  /** @brief creates a new matrix/vector object and returns an auto_ptr
   * pointing to the allocated object
   *
   * - Matrices have size @f$ H \times H @f$ where @f$H@f$ is the number
   * of degrees of freedom in the discrete function space @f$ { \cal X }_H @f$.
   * The matrices' sparsity pattern is determined by the discrete function
   * space's basefunction overlap.
   * - Vectors have @f$ H @f$ components
   *
   * @param dfs the discrete function space @f$ { \cal X }_H @f$.
   */
  template< class DiscFuncSpace >
  static AutoPtrType create( DiscFuncSpace& dfs )
  {
    DUNE_THROW( InvalidStateException, "Factory not implemented for ContainerType!" );
  }

  /** @brief creates a new matrix/vector object and returns a pointer to the
   * allocated object
   *
   * - Matrices have size @f$ H \times H @f$ where @f$H@f$ is the number
   * of degrees of freedom in the discrete function space @f$ { \cal X }_H @f$.
   * The matrices' sparsity pattern is determined by the discrete function
   * space's basefunction overlap.
   * - Vectors have @f$ H @f$ components
   *
   * @param dfs the discrete function space @f$ { \cal X }_H @f$.
   */
  template< class DiscFuncSpace >
  static ContainerType* createPtr( DiscFuncSpace& dfs )
  {
    DUNE_THROW( InvalidStateException, "Factory not implemented for ContainerType!" );
  }
}; // end of class Factory

/** @brief interface of static factory class for matrix classes
 */
template< class ContainerImp >
class MatrixFactory
  : public Factory< ContainerImp >
{
};

// specialization for BCRSMatrix
template<class T>
class MatrixFactory< Dune::BCRSMatrix< T > >//Dune::FieldMatrix<double, 1,1> > >
{
public:
//  typedef Dune::FieldMatrix<double, 1,1> T;
  typedef Dune::BCRSMatrix<T>
    ContainerType;
  typedef std::auto_ptr<ContainerType>
    AutoPtrType;

public:


  /** @brief creates a new BCRSMatrix object and returns an auto_ptr pointing
   * to the allocated object
   *
   * Matrices have size @f$ H \times H @f$ where @f$H@f$ is the number
   * of degrees of freedom in the discrete function space @f$ { \cal X }_H @f$.
   * The matrices' sparsity pattern is determined by the discrete function
   * space's basefunction overlap.
   *
   * @param dfs the discrete function space @f$ { \cal X }_H @f$.
   */
  template< class DiscFuncSpace >
  static AutoPtrType create( DiscFuncSpace& dfs )
  {
    return AutoPtrType( createPtr( dfs ) );
  }

  /** @brief creates a new BCRSMatrix object and returns a pointer to the
   * allocated object
   *
   * Matrices have size @f$ H \times H @f$ where @f$H@f$ is the number
   * of degrees of freedom in the discrete function space @f$ { \cal X }_H @f$.
   * The matrices' sparsity pattern is determined by the discrete function
   * space's basefunction overlap.
   *
   * @param dfs the discrete function space @f$ { \cal X }_H @f$.
   */
  template< class DiscFuncSpace >
  static ContainerType* createPtr( DiscFuncSpace& dfs )
  {
    typedef typename DiscFuncSpace::BaseFunctionSetType
      BFS;
    typedef typename DiscFuncSpace::IteratorType
      ItType;
    typedef typename ItType::Entity
      Entity;

    const unsigned int numDofs = dfs.size();
    typedef Dune::Functionals::Container::SparsityPattern
      PatternType;
    PatternType sPattern(numDofs);

    ContainerType* matrix = new ContainerType(numDofs,
                                              numDofs,
                                              ContainerType::random);

    // compute sparsity pattern
    // \todo precompile this in linear subspace
    // \todo use constraints for sparsity pattern
    ItType it = dfs.begin();
    for ( ; it != dfs.end(); ++it )
    {
      const Entity& en = *it;
      const BFS& bfs = dfs.baseFunctionSet( en );

      for (unsigned int i = 0; i < bfs.numBaseFunctions(); ++i)
      {
        unsigned int ii = dfs.mapToGlobal( en, i );
        for (unsigned int j = 0; j < bfs.numBaseFunctions(); ++j)
        {
          unsigned int jj = dfs.mapToGlobal( en, j );
          sPattern.insert( ii, jj );
        }
      }
    }

    for (unsigned int i = 0; i < sPattern.size(); ++i) {
      matrix->setrowsize( i, sPattern.countNonZeros( i ) );
    }
    matrix->endrowsizes();

    for (unsigned int i = 0; i < sPattern.size(); ++i)
    {
      typedef SparsityPattern::NonZeroColIterator
        ColIterator;
      ColIterator sit = sPattern.begin( i );
      for ( ; sit!=sPattern.end( i ); sit++ )
      {
        matrix->addindex(i, *sit);
      }
    }
    matrix->endindices();

    return matrix;
  }

}; // end of MatrixFactory<BCRSMatrix<T> >


/** @brief interface of static factory class for vector classes
 */
template<class ContainerImp>
class VectorFactory
  : public Factory<ContainerImp>
{
};

/** @brief static factory class for vector classes of type Dune::BlockVector
 */
template< class T >
class VectorFactory< Dune::BlockVector< T > >
{
public:
  //! \copydoc Factory::ContainerType
  typedef Dune::BlockVector< T >
    ContainerType;
  //! \copydoc Factory::AutoPtrType;
  typedef std::auto_ptr< ContainerType >
    AutoPtrType;

public:
  /** @brief creates a new vector object and returns an auto_ptr pointing to
   * the allocated object
   *
   * The vector has @f$ H @f$ components which is the number of degrees of
   * freedom of the given discrete function space @f$ {\cal X}_H @f$.
   *
   * @param dfs the discrete function space @f$ { \cal X }_H @f$.
   */
  template< class DFSType >
  static ContainerType* createPtr( DFSType& dfs )
  {
    const unsigned int numDofs = dfs.size();
    ContainerType* bv = new ContainerType(numDofs);
    return bv;
  }

  template< class DFSType >
  static AutoPtrType create( DFSType & dfs)
  {
    return AutoPtrType( createPtr( dfs ) );
  }
}; // end of VectorFactory<BlockVector<T> >

} // end of namespace Container

} // end of namespace Functionals

} // end of namespace Dune


#endif /* end of include guard: DUNE_FEM_FUNCTIONALS_CONTAINER_FACTORY_HH */
