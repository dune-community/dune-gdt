#ifndef DUNE_FEM_FUNCTIONALS_SUBSPACE_SUBSPACES_HH
#define DUNE_FEM_FUNCTIONALS_SUBSPACE_SUBSPACES_HH


namespace Dune
{

namespace Functionals
{

namespace Subspace
{

// \todo should prepare sparsity patterns and such things!
template< class DiscFuncSpace, class Constraints >
class Linear
{
public:
  typedef DiscFuncSpace
    DiscreteFunctionSpaceType;

  typedef Constraints
    ConstraintsType;
public:
  Linear( DiscreteFunctionSpaceType& space,
          ConstraintsType& constraints )
    : space_(space),
      constraints_(constraints)
  {
  }
private:
  DiscreteFunctionSpaceType &space_;
  ConstraintsType           &constraints_;
}; // end of class Linear

template< class LinearSubspace, class OffsetFunction >
class Affine
{
public:
  typedef LinearSubspace
    LinearSubSpaceType;

  typedef OffsetFunction
    OffsetFunctionType;
public:
  Affine( LinearSubspace& linear,
          OffsetFunctionType& offset )
    : linear_(linear),
      offset_(offset)
  {
  }
private:
  LinearSubspace     &linear_;
  OffsetFunctionType &offset_;
}; // end of class Affine

} // end of namespace Subspace

} // end of namespace Functionals

} // end of namespace Dune


#endif /* end of include guard: DUNE_FEM_FUNCTIONALS_SUBSPACE_SUBSPACES_HH */
