#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH

// dune-functionals includes
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

  typedef typename SuperSpaceType::LocalBaseFunctionProviderType
    LocalBaseFunctionProviderType;

  typedef Dune::Functionals::Constraints::DirichletZero< SuperSpaceType >
    DirichletZeroConstraintsType;

  typedef typename SuperSpaceType::GridPartType
    GridPartType;

  typedef typename SuperSpaceType::FunctionSpaceType
    FunctionSpaceType;

  typedef typename SuperSpaceType::RangeFieldType
    RangeFieldType;

  typedef typename SuperSpaceType::DomainType
    DomainType;

  typedef typename SuperSpaceType::EntityType
    EntityType;

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

  const DirichletZeroConstraintsType& constraints() const
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

  const int numMaxDoFs() const
  {
    return superSpace_.numMaxDoFs();
  }

  const int order() const
  {
    return superSpace_.order();
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
  const DirichletZeroConstraintsType constraints_;

}; // end of class DirichletZero

} // end namespace Linear

} // end of namespace Subspace

} // end namespace DiscreteFunctionSpace

} // end of namespace Functionals

} // end of namespace Dune


#endif /* end of include guard: DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH */
