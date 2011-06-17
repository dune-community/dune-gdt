#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_AFFINE_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_AFFINE_HH

// dune-fem includes
#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>

// dune-fem-tools includes
#include <dune/fem-tools/function/runtimefunction.hh>
#include <dune/fem-tools/space/projection.hh>

namespace Dune {

namespace Functionals {

namespace DiscreteFunctionSpace {

namespace Subspace {

namespace Affine {

template <class DirichletZeroSpaceImp>
class Dirichlet
{
public:
  typedef DirichletZeroSpaceImp DirichletZeroSpaceType;

  typedef typename DirichletZeroSpaceType::SuperSpaceType SuperSpaceType;

  typedef typename DirichletZeroSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef Dune::FemTools::Function::Runtime<FunctionSpaceType> RuntimeFunctionType;

  typedef Dune::AdaptiveDiscreteFunction<typename SuperSpaceType::HostSpaceType> AffineShiftType;

  typedef typename DirichletZeroSpaceType::EntityType EntityType;

  typedef typename DirichletZeroSpaceType::RangeFieldType RangeFieldType;

  typedef typename DirichletZeroSpaceType::DomainType DomainType;

  /**
    \defgroup dune-fem related
    \{
    **/
  typedef typename SuperSpaceType::BaseFunctionSetType BaseFunctionSetType;

  typedef typename SuperSpaceType::IteratorType IteratorType;
  /**
    \}
    **/

  Dirichlet(const DirichletZeroSpaceType& dirichletZeroSpace, const std::string expression = "[0.0;0.0;0.0]")
    : dirichletZeroSpace_(dirichletZeroSpace)
    , runtimeFunction_(expression)
    , affineShift_("affineShift", dirichletZeroSpace_.superSpace().hostSpace())
  {
    Dune::FemTools::Projection::Dirichlet::project(runtimeFunction_, affineShift_);
  }

  const DirichletZeroSpaceType& baseSpace() const
  {
    return dirichletZeroSpace_;
  }

  const SuperSpaceType& superSpace() const
  {
    return dirichletZeroSpace_.superSpace();
  }

  const AffineShiftType& affineShift() const
  {
    return affineShift_;
  }

  const int numMaxDoFs() const
  {
    return dirichletZeroSpace_.numMaxDoFs();
  }

  const int order() const
  {
    return dirichletZeroSpace_.order();
  }

  /**
    \defgroup dune-fem related
    \{
    **/
  IteratorType begin() const
  {
    return dirichletZeroSpace_.begin();
  }

  const IteratorType end() const
  {
    return dirichletZeroSpace_.end();
  }

  const BaseFunctionSetType baseFunctionSet(const EntityType& entity) const
  {
    return dirichletZeroSpace_.baseFunctionSet(entity);
  }

  const int mapToGlobal(const EntityType& entity, const int localDof) const
  {
    return dirichletZeroSpace_.mapToGlobal(entity, localDof);
  }
  /**
    \}
    **/

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
