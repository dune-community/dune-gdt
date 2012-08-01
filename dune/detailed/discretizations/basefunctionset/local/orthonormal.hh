#ifndef DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_LOCAL_ORTHONORMAL_HH
#define DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_LOCAL_ORTHONORMAL_HH

namespace Dune
{

namespace Detailed {

namespace Discretizations
{

namespace BaseFunctionSet
{

namespace Local
{

template< class BaseFunctionSetImp >
class Orthonormal
{
public:

  typedef BaseFunctionSetImp
    BaseFunctionSetType;

  typedef typename BaseFunctionSetType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType
    FunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::GridElementType
    GridElementType;

  enum{ polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder };

  typedef Orthonormal< BaseFunctionSetType >
    ThisType;

  typedef typename FunctionSpaceType::DomainFieldType
    DomainFieldType;

  typedef typename FunctionSpaceType::DomainType
    DomainType;

  typedef typename FunctionSpaceType::RangeFieldType
    RangeFieldType;

  typedef typename FunctionSpaceType::RangeType
    RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType
    JacobianRangeType;

  typedef typename FunctionSpaceType::HessianRangeType
    HessianRangeType;

private:
  typedef typename BaseFunctionSetType::HostLocalBaseFunctionSetType
    HostLocalBaseFunctionSetType;

public:
  //! constructor
  Orthonormal( const BaseFunctionSetType& baseFunctionSet, const GridElementType& gridElement )
    : baseFunctionSet_( baseFunctionSet ),
      gridElement_( gridElement ),
      hostLocalBaseFunctionSet_( baseFunctionSet_.hostLocalBaseFunctionSet_ ),
      size_( hostLocalBaseFunctionSet_.size() ),
      order_( hostLocalBaseFunctionSet_.order() )
  {
    hostLocalBaseFunctionSet_.size();
  }

  //! copy constructor
  Orthonormal( const ThisType& other )
    : baseFunctionSet_( other.baseFunctionSet() ),
      gridElement_( other.gridElement() ),
      hostLocalBaseFunctionSet_( other.hostLocalBaseFunctionSet_ ),
      size_( other.size() ),
      order_( other.order() )
  {
  }

  const BaseFunctionSetType& baseFunctionSet() const
  {
    return baseFunctionSet_;
  }

  const GridElementType& gridElement() const
  {
    return gridElement_;
  }

  unsigned int size() const
  {
    return size_;
  }

  int order() const
  {
    return order_;
  }

  void evaluate( const DomainType& x, std::vector< RangeType >& ret) const
  {
    assert( ret.size() == size_ );
    hostLocalBaseFunctionSet_.evaluateFunction( x, ret );
  }

  void jacobian( const DomainType& x, std::vector< JacobianRangeType >& ret ) const
  {
    assert( ret.size() == size_ );

    // some types we will need
    typedef typename GridElementType::Geometry
      ElementGeometryType;
    typedef typename ElementGeometryType::Jacobian
      JacobianInverseTransposedType;
    typedef typename JacobianRangeType::row_type
      JacobianRowType;

    // geometry and jacobian inverse transposed
    const ElementGeometryType& elementGeometry = gridElement_.geometry();
    const JacobianInverseTransposedType& jacobianInverseTransposed = elementGeometry.jacobianInverseTransposed( x );

    hostLocalBaseFunctionSet_.evaluateJacobian( x, ret );

    // evaluate
    JacobianRowType tmp( 0.0 );
    for( unsigned int i = 0; i < size_; ++i )
    {
      // transpose for each dim of range
      const unsigned int dimRange = DiscreteFunctionSpaceType::dimRange;
      for( unsigned int row = 0; row < dimRange; ++row )
      {
        // transpose
        tmp = ret[i][row];
        jacobianInverseTransposed.mv( tmp, ret[i][row] );
      }
    }
  }

private:

  //! assignment operator
  ThisType& operator=( const ThisType& );

  const BaseFunctionSetType& baseFunctionSet_;
  const GridElementType& gridElement_;
  const HostLocalBaseFunctionSetType& hostLocalBaseFunctionSet_;
  const unsigned int size_;
  const int order_;

}; // end class Orthonormal

} // end namespace Local

} // end namespace Common

} // namespace Discretizations

} // namespace Detailed

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_LOCAL_ORTHONORMAL_HH
