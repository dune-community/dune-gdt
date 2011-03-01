#ifndef DUNE_FEM_FUNCTIONALS_DOFVECTOR_HH
#define DUNE_FEM_FUNCTIONALS_DOFVECTOR_HH

namespace Dune
{

namespace Functionals
{

/**
  * \brief      This class represents a local DoF vector.
  *
  *             It is based upon std::vector and should be replaced by something clever in the future!
  *
  * \todo       Doc me, please!
  **/
template< class ElementType >
class LocalDoFVector
{
public:
  /**
    * \brief    Initializes an empty vector, according to the given size.
    **/
  LocalDoFVector( const unsigned int size )
    : size_( size )
  {
    // resize
    storage_.resize( size, 0.0 );
  }

  /**
    * \brief    Initializes a DoF vector and sets its entries to the
    *           corresponding entries of the given localFunction.
    **/
  template< class LocalFunctionType >
  LocalDoFVector( const LocalFunctionType& localFunction )
    : size_( localFunction.numDofs() )
  {
    // resize
    storage_.resize( localFunction.numDofs() );

    // copy entries
    for(  int ii = 0;
          ii < localFunction.numDofs();
          ++ii )
    {
      storage_[ii] = localFunction[ii];
    }
  }

  /**
    * \brief    Returns the size.
    */
  const unsigned int size() const
  {
    return size_;
  }

  /**
    * \brief    Random read and write access.
    **/
  ElementType& operator[]( const unsigned int ii )
  {
    return storage_[ii];
  }

  /**
    * \brief    Random read access.
    **/
  const ElementType operator[]( const unsigned int ii ) const
  {
    return storage_[ii];
  }

  /**
    * \brief    Scalar product of two local DoF vectors of same type.
    **/
  ElementType operator*( const LocalDoFVector< ElementType >& other ) const
  {
    assert( size_ == other.size() );
    ElementType result = 0.0;

    for(  unsigned ii = 0;
          ii < size_;
          ++ii )
    {
      result += storage_[ii] * other[ii];
    }

    return result;
  }

private:
  std::vector< ElementType > storage_;
  const unsigned int size_;

}; // end class LocalDoFVector

} // end namespace Functionals

} // end namespace Dune

#endif // end DUNE_FEM_FUNCTIONALS_DOFVECTOR_HH
