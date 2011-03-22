#ifndef SUBSPACES_WBM7R4K5
#define SUBSPACES_WBM7R4K5


namespace Dune {
namespace Fem {
namespace Functional {
namespace Subspace {

// \todo should prepare sparsity patterns and such things!
template <class DiscFuncSpace, class Constraints>
class Linear
{
public:
  typedef DiscFuncSpace DiscreteFunctionSpaceType;

  typedef Constraints ConstraintsType;

public:
  Linear(DiscreteFunctionSpaceType& space, ConstraintsType& constraints)
    : space_(space)
    , constraints_(constraints)
  {
  }

private:
  DiscreteFunctionSpaceType& space_;
  ConstraintsType& constraints_;
}; // end of class Linear

template <class LinearSubspace, class OffsetFunction>
class Affine
{
public:
  typedef LinearSubspace LinearSubSpaceType;

  typedef OffsetFunction OffsetFunctionType;

public:
  Affine(LinearSubspace& linear, OffsetFunctionType& offset)
    : linear_(linear)
    , offset_(offset)
  {
  }

private:
  LinearSubspace& linear_;
  OffsetFunctionType& offset_;
}; // end of class Affine

} // end of namespace Constraints
} // end of namespace Functional
} // end of namespace Fem
} // end of namespace Dune


#endif /* end of include guard: SUBSPACES_WBM7R4K5 */
