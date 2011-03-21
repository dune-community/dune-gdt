#ifndef LOCALMATRIX_5K9BIZ34
#define LOCALMATRIX_5K9BIZ34



namespace Dune
{
namespace Fem
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
class LocalMatrix
{
public:
  /**
    * \brief    Initializes an empty vector, according to the given size.
    **/
  LocalMatrix( const unsigned int n, const unsigned int m )
    : n_( n ), m_(m)
  {
    // resize
    storage_.resize( n_ * m_, 0.0 );
  }

  /**
    * \brief    Returns the storage size (n*m).
    */
  const unsigned int size() const
  {
    return n_*m_;
  }

  ElementType & get( const unsigned int i, const unsigned int j )
  {
    return storage_[i*m_+j];
  }

  void set( const unsigned int i, const unsigned int j, ElementType val )
  {
    storage_[i*m_+j] = val;
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

private:
  std::vector< ElementType > storage_;
  const unsigned int n_;
  const unsigned int m_;

}; // end class LocalDoFVector

} // end namespace Functionals

} // end namespace Fem

} // end namespace Dune

#endif /* end of include guard: LOCALMATRIX_5K9BIZ34 */
