#ifndef DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_SUBSPACE_AFFINE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_SUBSPACE_AFFINE_HH

// dune-helper-tools includes
#include <dune/helper-tools/function/runtime.hh>

// dune-detailed-discretizations includes
#include <dune/detailed-discretizations/basefunctionset/local/lagrange.hh>
#include <dune/detailed-discretizations/discretefunction/continuous.hh>

namespace Dune {

namespace DetailedDiscretizations {

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

  typedef typename BaseSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef typename BaseSpaceType::GridPartType GridPartType;

  enum
  {
    polynomialOrder = BaseSpaceType::polynomialOrder
  };

private:
  typedef Dune::HelperTools::Function::Runtime<FunctionSpaceType> RuntimeFunctionType;

public:
  typedef Dune::DetailedDiscretizations::DiscreteFunction::Continuous::BlockVector<SuperSpaceType> AffineShiftType;

  typedef typename BaseSpaceType::ConstraintsType ConstraintsType;

  typedef typename BaseSpaceType::BaseFunctionSetType BaseFunctionSetType;

  typedef typename BaseSpaceType::DomainType DomainType;

  typedef typename BaseSpaceType::DomainFieldType DomainFieldType;

  typedef typename BaseSpaceType::RangeFieldType RangeFieldType;

  typedef typename BaseSpaceType::RangeType RangeType;

  typedef typename BaseSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename BaseSpaceType::HessianRangeType HessianRangeType;

  typedef typename BaseSpaceType::MapperType MapperType;

  static const unsigned int dimDomain = BaseSpaceType::dimDomain;

  static const unsigned int dimRange = BaseSpaceType::dimRange;

  /**
      @name Convenience Types
      @{
   **/
  typedef typename SuperSpaceType::IteratorType IteratorType;

  typedef typename SuperSpaceType::EntityType EntityType;
  /**
      @}
   **/

  /**
    \brief  constructor, taking an analytical expression for the boundary data, does a dirichlet projection to generate
            the affine shift
    **/
  Dirichlet(const BaseSpaceType& baseSpace, const std::string expression = "[0.0;0.0;0.0]")
    : baseSpace_(baseSpace)
    , affineShift_(baseSpace.superSpace(), "affineShift", RuntimeFunctionType(expression), "dirichlet")
  {
  }

  /**
    \brief  constructor, taking a discrete function as affine shift
    **/
  Dirichlet(const BaseSpaceType& baseSpace, const AffineShiftType& affineShift)
    : baseSpace_(baseSpace)
    , affineShift_(affineShift)
  {
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

  const GridPartType& gridPart() const
  {
    return baseSpace_.gridPart();
  }

  const BaseFunctionSetType& baseFunctionSet() const
  {
    return baseSpace_.baseFunctionSet();
  }

  int order() const
  {
    return baseSpace_.order();
  }

  const ConstraintsType& constraints() const
  {
    return baseSpace_.constraints();
  }

  const MapperType& map() const
  {
    return baseSpace_.map();
  }

  bool continuous() const
  {
    return baseSpace_.continuous();
  }

  /**
      @name Convenience methods
      @{
   **/
  IteratorType begin() const
  {
    return baseSpace_.gridPart().template begin<0>();
  }

  IteratorType end() const
  {
    return baseSpace_.gridPart().template end<0>();
  }
  /**
      @}
   **/

private:
  //! copy constructor
  Dirichlet(const ThisType&);

  //! assignment operator
  ThisType& operator=(const ThisType&);

  const BaseSpaceType& baseSpace_;
  const AffineShiftType affineShift_;
}; // end class Dirichlet

} // end namespace Affine

} // end namespace Subspace

} // end namespace DiscreteFunctionSpace

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_SUBSPACE_AFFINE_HH
