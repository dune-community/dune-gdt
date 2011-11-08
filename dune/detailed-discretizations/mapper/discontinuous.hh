#ifndef DUNE_DETAILED_DISCRETIZATIONS_MAPPER_DISCONTINUOUS_HH
#define DUNE_DETAILED_DISCRETIZATIONS_MAPPER_DISCONTINUOUS_HH

namespace Dune
{

namespace DetailedDiscretizations
{

namespace Mapper
{

/**
  \attention assumes that we have same number of local basefunctions on each element and only one element type in the
             grid view
  **/
template< class DiscreteFunctionSpaceImp >
class Discontinuous
{
public:

  typedef DiscreteFunctionSpaceImp
    DiscreteFunctionSpaceType;

  typedef Discontinuous< DiscreteFunctionSpaceType >
    ThisType;

  typedef typename DiscreteFunctionSpaceType::GridViewType::IndexSet
    IndexSetType;

  enum{ polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder };

  Discontinuous( const DiscreteFunctionSpaceType& discreteFunctionSpace )
    : discreteFunctionSpace_( discreteFunctionSpace ),
      indexSet_( discreteFunctionSpace_.gridView().indexSet() ),
      maxLocalSize_( discreteFunctionSpace.baseFunctionSet().local( *(discreteFunctionSpace.gridElementBegin()) ).size() ),
      size_( maxLocalSize()*indexSet_.size(0) )
  {
  }

  ~Discontinuous()
  {
  }

  const DiscreteFunctionSpaceType& discreteFunctionSpace()
  {
    return discreteFunctionSpace_;
  }

  template< class EntityType >
  unsigned int toGlobal( const EntityType& entity, const unsigned int localDofNumber ) const
  {
//    std::cout << indexSet_.index(entity) << ", " << localDofNumber << ": " << indexSet_.index(entity)*maxLocalSize() + localDofNumber << std::endl;
    return indexSet_.index(entity)*maxLocalSize() + localDofNumber;
  }

  unsigned int size() const
  {
    return size_;
  }

  unsigned int maxLocalSize() const
  {
    return maxLocalSize_;
  }

private:
  //! copy constructor
  Discontinuous( const ThisType& );

  //! assignment operator
  ThisType& operator=( const ThisType& );

  const DiscreteFunctionSpaceType& discreteFunctionSpace_;
  const IndexSetType& indexSet_;
  const unsigned int maxLocalSize_;
  const unsigned int size_;

}; // end class Discontinuous

} // end namespace Mapper

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_MAPPER_DISCONTINUOUS_HH
