#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_AFFINE_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_AFFINE_HH

// dune-fem includes
#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>

// dune-fucntionals includes
#include <dune/functionals/basefunctionset/local/lagrange.hh>
#include <dune/functionals/discretefunction/continuous.hh>

// dune-fem-tools includes
#include <dune/fem-tools/function/functiontools.hh>
#include <dune/fem-tools/function/runtimefunction.hh>
#include <dune/fem-tools/space/projection.hh>

namespace Dune {

namespace Functionals {

namespace DiscreteFunctionSpace {

namespace Subspace {

namespace Affine {

template <class BaseSpaceImp>
class Dirichlet
{
public:
  typedef BaseSpaceImp BaseSpaceType;

  typedef Dirichlet<BaseSpaceType> ThisType;

  typedef typename BaseSpaceType::SuperSpaceType SuperSpaceType;

  typedef typename BaseSpaceType::GridPartType GridPartType;

  typedef typename BaseSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef Dune::FemTools::Function::Runtime<FunctionSpaceType> RuntimeFunctionType;

  //  typedef Dune::AdaptiveDiscreteFunction< typename SuperSpaceType::HostSpaceType >
  typedef Dune::Functionals::DiscreteFunction::Continuous::BlockVector<SuperSpaceType> AffineShiftType;

  typedef typename BaseSpaceType::ConstraintsType ConstraintsType;

  typedef typename BaseSpaceType::EntityType EntityType;

  typedef typename BaseSpaceType::DomainType DomainType;

  typedef typename BaseSpaceType::DomainFieldType DomainFieldType;

  typedef typename BaseSpaceType::RangeFieldType RangeFieldType;

  typedef typename BaseSpaceType::RangeType RangeType;

  typedef typename BaseSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename BaseSpaceType::HessianRangeType HessianRangeType;

  typedef Dune::Functionals::Common::LocalBaseFunctionSet<ThisType> LocalBaseFunctionSetType;

  /**
    \defgroup dune-fem related
    \{
    **/
  typedef typename SuperSpaceType::BaseFunctionSetType BaseFunctionSetType;

  typedef typename SuperSpaceType::IteratorType IteratorType;
  /**
    \}
    **/

  static const unsigned int dimDomain = BaseSpaceType::dimDomain;

  static const unsigned int dimRange = BaseSpaceType::dimRange;

  Dirichlet(const BaseSpaceType& baseSpace, const std::string expression = "[0.0;0.0;0.0]")
    : baseSpace_(baseSpace)
    , runtimeFunction_(expression)
    ,
    //      affineShift_( "affineShift", baseSpace.superSpace().hostSpace() )
    affineShift_(baseSpace.superSpace(), "affineShift", runtimeFunction_, "dirichlet")
  {
    //    Dune::FemTools::Projection::Dirichlet::project( runtimeFunction_, affineShift_ );
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

  const LocalBaseFunctionSetType localBaseFunctionSet(const EntityType& entity) const
  {
    return LocalBaseFunctionSetType(*this, entity);
  }

  const int numMaxLocalDoFs() const
  {
    return baseSpace_.numMaxLocalDoFs();
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
  const AffineShiftType affineShift_;
}; // end class Dirichlet

} // end of namespace Affine

} // end namespace Subspace

} // end namespace DiscreteFunctionSpace

} // end of namespace Functionals

} // end of namespace Dune

#endif // DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_SUBSPACE_AFFINE_HH
