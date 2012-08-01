#ifndef DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH
#define DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH

// dune-detailed-discretizations includes
//#include <dune/detailed/discretizations/basefunctionset/local/lagrange.hh>
#include <dune/detailed/discretizations/constraints/dirichlet.hh>

namespace Dune
{

namespace Detailed
{

namespace Discretizations {

namespace DiscreteFunctionSpace
{

namespace Sub
{

namespace Linear
{

template< class SuperSpaceImp >
class Dirichlet
{
public:

  typedef SuperSpaceImp
    SuperSpaceType;

  typedef Dirichlet< SuperSpaceType >
    ThisType;

  typedef Dune::Detailed::Discretizations::Constraints::DirichletZero< SuperSpaceType >
    ConstraintsType;

  typedef typename SuperSpaceType::FunctionSpaceType
    FunctionSpaceType;

  typedef typename SuperSpaceType::GridPartType
    GridPartType;

  typedef typename SuperSpaceType::GridViewType
    GridViewType;

  enum{ polynomialOrder = SuperSpaceType::polynomialOrder };

  typedef typename SuperSpaceType::MapperType
    MapperType;

  typedef typename SuperSpaceType::BaseFunctionSetType
    BaseFunctionSetType;

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

  static const unsigned int dimDomain = SuperSpaceType::dimDomain;

  static const unsigned int dimRange = SuperSpaceType::dimRange;

  typedef typename SuperSpaceType::PatternType PatternType;

  /**
      @name Convenience
      @{
   **/
  typedef typename SuperSpaceType::IteratorType
    IteratorType;

  typedef typename SuperSpaceType::EntityType
    EntityType;
  /**
      @}
   **/

  Dirichlet( const SuperSpaceType& superSpace )
    : superSpace_( superSpace ),
      constraints_( superSpace_ )
  {
  }

  const SuperSpaceType& superSpace() const
  {
    return superSpace_;
  }

  const ConstraintsType& constraints() const
  {
    return constraints_;
  }

  const GridPartType& gridPart() const
  {
    return superSpace_.gridPart();
  }

  const GridViewType& gridView() const
  {
    return superSpace_.gridView();
  }

  const MapperType& map() const
  {
    return superSpace_.map();
  }

  const BaseFunctionSetType& baseFunctionSet() const
  {
    return superSpace_.baseFunctionSet();
  }

  int order() const
  {
    return superSpace_.order();
  }

  bool continuous() const
  {
    return superSpace_.continuous();
  }

  /**
      @name Convenience methods
      @{
   **/
  IteratorType begin() const
  {
    return superSpace_.gridPart().template begin< 0 >();
  }

  IteratorType end() const
  {
    return superSpace_.gridPart().template end< 0 >();
  }
  /**
      @}
   **/

  template< class OtherDiscreteFunctionSpaceType>
  PatternType computePattern(const OtherDiscreteFunctionSpaceType& other) const
  {
    return superSpace_.computePattern(other);
  }

  PatternType computePattern() const
  {
    return superSpace_.computePattern();
  }

private:

  //! copy constructor
  Dirichlet( const ThisType& );

  //! assignment operator
  ThisType& operator=( const ThisType& );

  const SuperSpaceType& superSpace_;
  const ConstraintsType constraints_;

}; // end class Dirichlet

} // end namespace Linear

} // end namespace Sub

} // end namespace DiscreteFunctionSpace

} // namespace Discretizations

} // end namespace Discretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH
