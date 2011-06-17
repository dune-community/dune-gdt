#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_AFFINE_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_AFFINE_HH

//dune-fem includes
#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>

namespace Dune
{

namespace Functionals
{

namespace DiscreteFunctionSpace
{

namespace Subspace
{

namespace Affine
{

template< class DirichletZeroSpaceImp >
class Dirichlet
{
public:

  typedef DirichletZeroSpaceImp
    DirichletZeroSpaceType;

  typedef Dune::AdaptiveDiscreteFunction< typename DirichletZeroSpaceType::SuperSpaceType >
    AffineShiftType;

  Dirichlet(  const DirichletZeroSpaceType& dirichletZeroSpace/*,
              const std::string xExpression,
              const std::string yExpression = "",
              const std::string zExpression = ""*/)
    : dirichletZeroSpace_( dirichletZeroSpace )/*,
      affineShift_( xExpression, yExpression, zExpression )*/
  {
  }

private:

  const DirichletZeroSpaceType& dirichletZeroSpace_;
//  const AffineShiftType affineShift_;
}; // end class Dirichlet

} // end of namespace Affine

} // end namespace Subspace

} // end namespace DiscreteFunctionSpace

} // end of namespace Functionals

} // end of namespace Dune

#endif // DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_AFFINE_HH
