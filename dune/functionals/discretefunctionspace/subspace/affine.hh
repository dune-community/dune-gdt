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
  typedef DirichletZeroSpaceImp BaseSpaceType;

  typedef typename BaseSpaceType::SuperSpaceType SuperSpaceType;

  typedef typename BaseSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef Dune::FemTools::Function::Runtime<FunctionSpaceType> RuntimeFunctionType;

  typedef Dune::AdaptiveDiscreteFunction<typename SuperSpaceType::HostSpaceType> AffineShiftType;

  typedef typename BaseSpaceType::ConstraintsType ConstraintsType;

  typedef typename BaseSpaceType::EntityType EntityType;

  typedef typename BaseSpaceType::RangeFieldType RangeFieldType;

  typedef typename BaseSpaceType::DomainType DomainType;

  /**
    \defgroup dune-fem related
    \{
    **/
  typedef typename SuperSpaceType::BaseFunctionSetType BaseFunctionSetType;

  typedef typename SuperSpaceType::IteratorType IteratorType;
  /**
    \}
    **/

  Dirichlet(const BaseSpaceType& baseSpace, const std::string expression = "[0.0;0.0;0.0]")
    : baseSpace_(baseSpace)
    , runtimeFunction_(expression)
    , affineShift_("affineShift", baseSpace.superSpace().hostSpace())
  {
    Dune::FemTools::Projection::Dirichlet::project(runtimeFunction_, affineShift_);
  }

  const BaseSpaceType& baseSpace() const
  {
    return baseSpace_;
  }

  const SuperSpaceType& superSpace() const
  {
    return baseSpace_.superSpace();
  }

  const AffineShiftType& affineShift() const
  {
    return affineShift_;
  }

  const int numMaxDoFs() const
  {
    return baseSpace_.numMaxDoFs();
  }

  const int order() const
  {
    return baseSpace_.order();
  }

  const ConstraintsType& constraints() const
  {
    return baseSpace_.constraints();
  }

  /**
    \defgroup dune-fem related
    \{
    **/
  IteratorType begin() const
  {
    return baseSpace_.begin();
  }

  const IteratorType end() const
  {
    return baseSpace_.end();
  }

  const BaseFunctionSetType baseFunctionSet(const EntityType& entity) const
  {
    return baseSpace_.baseFunctionSet(entity);
  }

  const int mapToGlobal(const EntityType& entity, const int localDof) const
  {
    return baseSpace_.mapToGlobal(entity, localDof);
  }
  /**
    \}
    **/

private:
  const BaseSpaceType& baseSpace_;
  const RuntimeFunctionType runtimeFunction_;
  AffineShiftType affineShift_;
}; // end class Dirichlet

} // end of namespace Affine

} // end namespace Subspace

} // end namespace DiscreteFunctionSpace

} // end of namespace Functionals

} // end of namespace Dune

#endif // DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_AFFINE_HH
