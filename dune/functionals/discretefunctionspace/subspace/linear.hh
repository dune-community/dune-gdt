#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_LINEAR_HH

// dune-functionals includes
#include <dune/functionals/constraints/dirichlet.hh>

namespace Dune {

namespace Functionals {

namespace DiscreteFunctionSpace {

namespace Subspace {

namespace Linear {

template <class SuperSpaceImp>
class DirichletZero
{
public:
  typedef SuperSpaceImp SuperSpaceType;

  typedef Dune::Functionals::Constraints::DirichletZero<SuperSpaceType> DirichletZeroConstraintsType;

  typedef typename SuperSpaceType::GridPartType GridPartType;

  typedef typename SuperSpaceType::FunctionSpaceType FunctionSpaceType;

  DirichletZero(const SuperSpaceType& superSpace)
    : superSpace_(superSpace)
    , constraints_(superSpace_)
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
