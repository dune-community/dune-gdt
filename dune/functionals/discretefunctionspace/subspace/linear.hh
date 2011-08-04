#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH

// dune-functionals includes
#include <dune/functionals/basefunctionset/local/lagrange.hh>
#include <dune/functionals/constraints/dirichlet.hh>

namespace Dune
{

namespace Functionals
{

namespace DiscreteFunctionSpace
{

namespace Subspace
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

  typedef Dune::Functionals::Constraints::DirichletZero< SuperSpaceType >
    ConstraintsType;

  typedef typename SuperSpaceType::GridPartType
    GridPartType;

  typedef typename SuperSpaceType::FunctionSpaceType
    FunctionSpaceType;

  enum{ polynomialOrder = SuperSpaceType::polynomialOrder };

  typedef typename SuperSpaceType::MapperType
    MapperType;

  typedef typename FunctionSpaceType::DomainType
    DomainType;

  typedef typename FunctionSpaceType::DomainFieldType
    DomainFieldType;

  typedef typename FunctionSpaceType::RangeFieldType
    RangeFieldType;

  typedef typename FunctionSpaceType::RangeType
    RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType
    JacobianRangeType;

  typedef typename FunctionSpaceType::HessianRangeType
    HessianRangeType;

  typedef typename SuperSpaceType::LocalBaseFunctionSetType
    LocalBaseFunctionSetType;

  static const unsigned int dimDomain = SuperSpaceType::dimDomain;

  static const unsigned int dimRange = SuperSpaceType::dimRange;

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

  const MapperType& map() const
  {
    return superSpace_.map();
  }

  int order() const
  {
    return superSpace_.order();
  }

  template< class EntityType >
  const LocalBaseFunctionSetType localBaseFunctionSet( const EntityType& entity ) const
  {
    return superSpace_.localBaseFunctionSet( entity );
  }

  bool continuous() const
  {
    return superSpace_.continuous();
  }

private:
  const SuperSpaceType& superSpace_;
  const ConstraintsType constraints_;

}; // end class Dirichlet

} // end namespace Linear

} // end namespace Subspace

} // end namespace DiscreteFunctionSpace

} // end namespace Functionals

} // end namespace Dune


#endif /* end of include guard: DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH */
