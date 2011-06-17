#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_AFFINE_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_AFFINE_HH

// dune-fem includes
#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>

// dune-fem-tools includes
#include <dune/fem-tools/function/runtimefunction.hh>
#include <dune/fem-tools/space/projection.hh>

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

  typedef typename DirichletZeroSpaceType::SuperSpaceType
    SuperSpaceType;

  typedef typename DirichletZeroSpaceType::FunctionSpaceType
    FunctionSpaceType;

  typedef Dune::FemTools::RuntimeFunction< FunctionSpaceType >
    RuntimeFunctionType;

  typedef Dune::AdaptiveDiscreteFunction< typename SuperSpaceType::HostSpaceType >
    AffineShiftType;

  Dirichlet( const DirichletZeroSpaceType& dirichletZeroSpace, const std::string expression = "[0.0;0.0;0.0]" )
    : dirichletZeroSpace_( dirichletZeroSpace ),
      runtimeFunction_( expression ),
      affineShift_( "affineShift", dirichletZeroSpace_.superSpace().hostSpace() )
  {
    Dune::FemTools::Space::BetterL2Projection::project( runtimeFunction_, affineShift_ );
  }

  const DirichletZeroSpaceType& baseSpace() const
  {
    return dirichletZeroSpace_;
  }

  const SuperSpaceType& superSpace() const
  {
    return dirichletZeroSpace_.superSpace();
  }

private:

  const DirichletZeroSpaceType& dirichletZeroSpace_;
  const RuntimeFunctionType runtimeFunction_;
  AffineShiftType affineShift_;
}; // end class Dirichlet

} // end of namespace Affine

} // end namespace Subspace

} // end namespace DiscreteFunctionSpace

} // end of namespace Functionals

} // end of namespace Dune

#endif // DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_AFFINE_HH
