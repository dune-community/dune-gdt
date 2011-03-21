#ifndef DIRICHLET_SHETXCOC
#define DIRICHLET_SHETXCOC

namespace Dune {
namespace Fem {
namespace Functional {
namespace Constraints {

template <class DiscFuncSpace>
class Dirichlet
{
public:
  typedef DiscFuncSpace DiscreteFunctionSpace;
  typedef typename DiscreteFunctionSpace::GridPartType GridPartType;

public:
  Dirichlet(DiscFuncSpace& space)
    : space_(space)
    , gridPart_(space.gridPart())
  {
  }

private:
  DiscreteFunctionSpace& space_;
  GridPartType& gridPart_;
}; // end class Dirichlet

} // end of namespace Constraints
} // end of namespace Functional
} // end of namespace Fem
} // end of namespace Dune

#endif /* end of include guard: DIRICHLET_SHETXCOC */
