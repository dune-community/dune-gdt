#ifndef DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM1_VECTOR_HH
#define DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM1_VECTOR_HH

// std includes
#include <vector>

namespace Dune
{

namespace Detailed {

namespace Discretizations
{

namespace Assembler
{

namespace Local
{

namespace Codim1
{

/**
  \todo Add neumann boundary treatment
  \todo When adding neumann: think of numTmpObjectsRequired()!
  \todo Add penalty parameter
  **/
template< class LocalFunctionalImp >
class Vector
{
public:

  typedef LocalFunctionalImp
    LocalFunctionalType;

  typedef Vector< LocalFunctionalType >
    ThisType;

  typedef typename LocalFunctionalType::RangeFieldType
    RangeFieldType;

  //! constructor
  Vector( const LocalFunctionalType localFunctional )
    : localFunctional_( localFunctional )
  {
  }

private:
  //! copy constructor
  Vector( const ThisType& other )
    : localFunctional_( other.localFunctional() )
  {
  }

public:
  const LocalFunctionalType& localFunctional() const
  {
    return localFunctional_;
  }

  /**
    \todo Add neumann treatment here!
    **/
  std::vector< unsigned int > numTmpObjectsRequired() const
  {
    std::vector< unsigned int > ret( 2, 0 );
    // we require 1 tmp vector in this local assembler
    ret[0] = 1;
    // the functional itself requires that much local matrices
    ret[1] = localFunctional_.numTmpObjectsRequired();
    return ret;
  }

  template< class TestSpaceType,
            class EntityType,
            class SystemVectorType,
            class LocalVectorType >
  void assembleLocal( const TestSpaceType& testSpace,
                      const EntityType& entity,
                      SystemVectorType& systemVector,
                      std::vector< std::vector< LocalVectorType > >& tmpLocalVectorsContainer ) const
  {
    // get the local basefunction set
    typedef typename TestSpaceType::BaseFunctionSetType::LocalBaseFunctionSetType
      LocalTesBaseFunctionSetType;

    const LocalTesBaseFunctionSetType localTestBaseFunctionSet = testSpace.baseFunctionSet().local( entity );

    // check tmp local vectors
    assert( tmpLocalVectorsContainer.size() > 1 );
    std::vector< LocalVectorType >& tmpLocalVectors = tmpLocalVectorsContainer[0];
    if( tmpLocalVectors.size() < 1 )
    {
      tmpLocalVectors.resize( 1, LocalVectorType( testSpace.map().maxLocalSize(), RangeFieldType( 0.0 ) ) );
    }

    // some types
    typedef typename TestSpaceType::GridViewType
      GridViewType;

    typedef typename GridViewType::IntersectionIterator
      IntersectionIteratorType;

    typedef typename IntersectionIteratorType::Intersection
      IntersectionType;

    typedef typename IntersectionType::EntityPointer
      EntityPointerType;

    const GridViewType& gridView = testSpace.gridView();

    const IntersectionIteratorType lastIntersection = gridView.iend( entity );

    // do loop over all intersections
    for( IntersectionIteratorType intIt = gridView.ibegin( entity ); intIt != lastIntersection; ++intIt )
    {
      const IntersectionType& intersection = *intIt;

      if( !intersection.neighbor() && intersection.boundary() ) // if boundary intersection
      {
//        const unsigned int boundaryId = intersection.boundaryId();

//        // if dirichlet boundary intersection
//        if( boundaryId == 2 )
//        {
          localFunctional_.applyLocal(  localTestBaseFunctionSet,
                                        intersection,
                                        tmpLocalVectors[0],
                                        tmpLocalVectorsContainer[1] );

          // write local vector to global
          addToVector( testSpace, entity, tmpLocalVectors[0], systemVector );

//        } // end if dirichlet boundary intersection
//        else if( boundaryId == 3 ) // if neumann boundary intersection
//        {
//        } // end if neumann boundary intersection
      }// end if boundary intersection
    } // done loop over all intersections
  } // end method assembleLocal

private:

  //! assignment operator
  ThisType& operator=( const ThisType& );

  template< class TestSpaceType,
            class EntityType,
            class LocalVectorType,
            class SystemVectorType >
  void addToVector( const TestSpaceType& testSpace,
                    const EntityType& entity,
                    const LocalVectorType& localVector,
                    SystemVectorType& systemVector ) const
  {
    for( unsigned int j = 0; j < testSpace.baseFunctionSet().local( entity ).size(); ++j )
    {
      const unsigned int globalJ = testSpace.map().toGlobal( entity, j );

      systemVector[globalJ] += localVector[j];
    }
  } // end method addToVector

  const LocalFunctionalType localFunctional_;
}; // end class Vector

} // end namespace Codim1

} // end namespace Local

} // end namespace Assembler

} // namespace Discretizations

} // namespace Detailed

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_CODIM1_VECTOR_HH
