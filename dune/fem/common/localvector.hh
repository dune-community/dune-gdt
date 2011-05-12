#ifndef DUNE_FEM_FUNCTIONALS_COMMON_DOFVECTOR_HH
#define DUNE_FEM_FUNCTIONALS_COMMON_DOFVECTOR_HH

namespace Dune
{

namespace Functionals
{

namespace Common
{

/**
  * \brief  This class represents a local DoF vector.
  *         It is based upon std::vector and should be replaced by something clever in the future!
  *
  * \tparam ElementType
  *         Type of one element, usually double or RangeFieldType.
  **/
template< class ElementType >
class LocalVector
{
public:

  /**
    * \brief      Initializes an empty vector, according to the given size.
    *
    * \param[in]  size
    *             Size, the vector should have, usually the number of local DoFs.
    **/
  LocalVector( const unsigned int size )
  {
    // resize
    storage_.resize( size, 0.0 );
  }

  /**
    * \brief      Initializes a DoF vector and sets its entries to the
    *             corresponding entries of the given localFunction.
    *
    * \tparam     LocalFunctionType
    *             Type of the given local function, usually something that fulfills the Dune::LocalFunction interface.
    *
    * \param[in]  localFunction
    *             Local function, whose DoFs are going to be assigned to the vector. Has to provide the method
    *             numDofs() and read access via an operator[].
    **/
  template< class LocalFunctionType >
  LocalVector( const LocalFunctionType& localFunction )
  {
    // resize
    storage_.resize( localFunction.numDofs(), 0.0 );

    // copy entries
    for( int i = 0; i < localFunction.numDofs(); ++i )
    {
      storage_[i] = localFunction[i];
    }
  }

  /**
    * \brief      Returns the size.
    *
    * \param[out] const unsigned int
    *             Number of entries.
    */
  const unsigned int size() const
  {
    return storage_.size();
  }

  /**
    * \brief      Random read and write access.
    *
    * \param[in]  i
    *             Number of the element.
    *
    * \param[out] ElementType&
    *             Reference to the element at position i, writable.
    **/
  ElementType& operator[]( const unsigned int i )
  {
    return storage_[i];
  }

  /**
    * \brief      Random read access.
    * \param[in]  i
    *             Number of the element.
    *
    * \param[out] const ElementType
    *             Reference to the element at position i, readable only.
    **/
  const ElementType operator[]( const unsigned int i ) const
  {
    return storage_[i];
  }

  /**
    * \brief      Scalar product of two local vectors of same type.
    *
    * \param[in]  other
    *             The LocalVector, that will be multiplied with from the right.
    *
    * \param[out] ElementType
    *             Result of the scalar product.
    **/
  ElementType operator*( const LocalVector< ElementType >& other ) const
  {
    assert( storage_.size() == other.size() );
    ElementType result = 0.0;

    for(  unsigned ii = 0;
          ii < storage_.size();
          ++ii )
    {
      result += storage_[ii] * other[ii];
    }

    return result;
  }

private:

  std::vector< ElementType > storage_;

}; // end of class LocalDoFVector

} // end of namespace Common

} // end of namespace Functionals

} // end of namespace Dune

#endif // end DUNE_FEM_FUNCTIONALS_COMMON_DOFVECTOR_HH
