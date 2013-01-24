/**
  \file vector.hh
  \todo doc
  **/

#ifndef DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_COMBINED_VECTOR_HH
#define DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_COMBINED_VECTOR_HH

//std includes
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

namespace Combined
{

/**
  \brief  Combined local vector assembler.
  \todo   doc
  **/
template< class FirstImp, class SecondImp >
class Vector
{
public:
  typedef FirstImp
    FirstType;

  typedef SecondImp
    SecondType;

  typedef Vector< FirstType, SecondType >
    ThisType;

  /**
    \brief  Constructor.
    \todo   doc
    **/
  Vector( const FirstType& first, const SecondType& second )
    : first_( first ),
      second_( second )
  {
  }

private:
  /**
    \brief  Copy constructor.
    \todo   doc
    **/
  Vector( const ThisType& other )
    : first_( other.first() ),
      second_( other.second() )
  {
  }

public:
  /**
    \brief  Returns first local vector assembler.
    \todo   doc
    **/
  const FirstType& first() const
  {
    return first_;
  }

  /**
    \brief  Returns second local vector assembler.
    \todo   doc
    **/
  const SecondType& second() const
  {
    return second_;
  }

  /**
    \brief  Returns the number of max tmp storage, which are required.
    \todo   doc
    **/
  std::vector< unsigned int > numTmpObjectsRequired() const
  {
    const auto& numTmpObjectsRequiredFirst = first_.numTmpObjectsRequired();
    const auto& numTmpObjectsRequiredSecond = second_.numTmpObjectsRequired();
    return { std::max( numTmpObjectsRequiredFirst[0], numTmpObjectsRequiredSecond[0] ),
             std::max( numTmpObjectsRequiredFirst[1], numTmpObjectsRequiredSecond[1] )};
  }

  /**
    \brief  Main method.
    \todo   doc
    **/
  template< class TestSpaceType,
            class EntityType,
            class SystemVectorType,
            class LocalVectorType >
  void assembleLocal( const TestSpaceType& testSpace,
                      const EntityType& entity,
                      SystemVectorType& systemVector,
                      std::vector< std::vector< LocalVectorType > >& tmpLocalVectorsContainer ) const
  {
    first_.assembleLocal(   testSpace,
                            entity,
                            systemVector,
                            tmpLocalVectorsContainer );
    second_.assembleLocal(  testSpace,
                            entity,
                            systemVector,
                            tmpLocalVectorsContainer );
  } // end method assembleLocal

private:
  /**
    \brief  Assignment operator.
    \todo   doc
    **/
  ThisType& operator=( const ThisType& );

  const FirstType& first_;
  const SecondType& second_;
}; // end class Vector

} // end namespace Combined

} // end namespace Local

} // end namespace Assembler

} // namespace Discretizations

} // namespace Detailed

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_ASSEMLBER_LOCAL_COMBINED_VECTOR_HH
