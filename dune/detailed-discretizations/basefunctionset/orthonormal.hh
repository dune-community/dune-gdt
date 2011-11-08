#ifndef DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_ORTHONORMAL_HH
#define DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_ORTHONORMAL_HH

// dune-localfunctions includes
#include <dune/localfunctions/orthonormal.hh>

// dune-detailed-discretizations includes
#include <dune/detailed-discretizations/basefunctionset/local/orthonormal.hh>

namespace Dune {

namespace DetailedDiscretizations {

namespace BaseFunctionSet {

template <class DiscreteFunctionSpaceImp>
class Orthonormal
{
public:
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  typedef Orthonormal<DiscreteFunctionSpaceType> ThisType;

  typedef typename DiscreteFunctionSpaceType::GridElementType GridElementType;

  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::GridViewType GridViewType;

  enum
  {
    polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder
  };

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

  typedef Dune::DetailedDiscretizations::BaseFunctionSet::Local::Orthonormal<ThisType> LocalBaseFunctionSetType;

private:
  typedef Dune::OrthonormalLocalFiniteElement<GridViewType::dimension, DomainType, RangeType> LocalFiniteElementType;

  typedef typename LocalFiniteElementType::Traits::LocalBasisType HostLocalBaseFunctionSetType;

public:
  Orthonormal(const DiscreteFunctionSpaceType& space)
    : space_(space)
    , localFiniteElement_(space_.gridElementBegin()->type().id(), polynomialOrder)
    , hostLocalBaseFunctionSet_(localFiniteElement_.localBasis())
  {
  } // end constructor

  ~Orthonormal()
  {
  }

  const DiscreteFunctionSpaceType& space() const
  {
    return space_;
  }

  LocalBaseFunctionSetType local(const GridElementType& element) const
  {
    return LocalBaseFunctionSetType(*this, element);
  }

private:
  //! copy constructor
  Orthonormal(const ThisType& other);

  //! assignment operator
  ThisType& operator=(const ThisType&);

  friend class Dune::DetailedDiscretizations::BaseFunctionSet::Local::Orthonormal<ThisType>;

  const DiscreteFunctionSpaceType& space_;
  const LocalFiniteElementType localFiniteElement_;
  const HostLocalBaseFunctionSetType& hostLocalBaseFunctionSet_;

}; // end class Orthonormal

} // end namespace BaseFunctionSet

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_ORTHONORMAL_HH
