#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH

// dune-functionals includes
#include <dune/functionals/common/localbasefunction.hh>
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
class DirichletZero
{
public:

  typedef SuperSpaceImp
    SuperSpaceType;

  typedef DirichletZero< SuperSpaceType >
    ThisType;

  typedef Dune::Functionals::Constraints::DirichletZero< SuperSpaceType >
    ConstraintsType;

  typedef typename SuperSpaceType::GridPartType
    GridPartType;

  typedef typename SuperSpaceType::FunctionSpaceType
    FunctionSpaceType;

  typedef typename SuperSpaceType::EntityType
    EntityType;

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

  typedef Dune::Functionals::Common::LocalBaseFunctionSet< ThisType >
    LocalBaseFunctionSetType;

  /**
    \defgroup dune-fem related
    \{
    **/
  typedef typename SuperSpaceType::BaseFunctionSetType
    BaseFunctionSetType;

  typedef typename SuperSpaceType::IteratorType
    IteratorType;
  /**
    \}
    **/

  DirichletZero( const SuperSpaceType& superSpace )
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

  const unsigned int size() const
  {
    return superSpace_.size();
  }

  const int numMaxLocalDoFs() const
  {
    return superSpace_.numMaxLocalDoFs();
  }

  const int order() const
  {
    return superSpace_.order();
  }

  const LocalBaseFunctionSetType localBaseFunctionSet( const EntityType& entity ) const
  {
    return LocalBaseFunctionSetType( *this, entity );
  }

  /**
    \defgroup dune-fem related
    \{
    **/
  IteratorType begin() const
  {
    return superSpace_.begin();
  }

  const IteratorType end() const
  {
    return superSpace_.end();
  }

  const BaseFunctionSetType baseFunctionSet( const EntityType& entity ) const
  {
    return superSpace_.baseFunctionSet( entity );
  }

  const int mapToGlobal( const EntityType& entity, const int localDof) const
  {
    return superSpace_.mapToGlobal( entity, localDof);
  }
  /**
    \}
    **/

private:
  const SuperSpaceType& superSpace_;
  const ConstraintsType constraints_;

}; // end of class DirichletZero

} // end namespace Linear

} // end of namespace Subspace

} // end namespace DiscreteFunctionSpace

} // end of namespace Functionals

} // end of namespace Dune


#endif /* end of include guard: DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH */
