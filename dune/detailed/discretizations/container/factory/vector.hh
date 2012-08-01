#ifndef DUNE_DETAILED_DISCRETIZATIONS_CONTAINER_FACTORY_VECTOR_HH
#define DUNE_DETAILED_DISCRETIZATIONS_CONTAINER_FACTORY_VECTOR_HH

// dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/fvector.hh>

// dune-istl includes
#include <dune/istl/bvector.hh>

// dune-detailed-discretizations includes
#include <dune/detailed/discretizations/container/sparsitypattern.hh>

namespace Dune
{

namespace Detailed {

namespace Discretizations
{

//! Contains container classes for storing.
namespace Container
{

namespace Vector
{

/**
 * @brief Interface of static factory class.
 */
template< class ContainerImp >
class Factory
{
};

/**
 * @brief Static factory class for vector classes of type Dune::BlockVector.
 *
 * @tparam T Type to construct a Dune::BlockVector representing the type for
 * a block, normally a Dune::FieldVector.
 */
template< class BlockType >
class Factory< Dune::BlockVector< BlockType > >
{
public:
  //! \copydoc Factory::ContainerType
  typedef Dune::BlockVector< BlockType >
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
}; // end of VectorFactory

template< class FieldType, int n >
class Defaults
{
public:

  typedef Factory< Dune::BlockVector< Dune::FieldVector< FieldType, n > > >
    BlockVector;

}; // end class Defaults

} // end namespace Vector

} // end namespace Container

} // namespace Discretizations

} // namespace Detailed

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_CONTAINER_FACTORY_VECTOR_HH
