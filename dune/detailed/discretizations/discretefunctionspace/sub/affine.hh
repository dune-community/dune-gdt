#ifndef DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_SUBSPACE_AFFINE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_SUBSPACE_AFFINE_HH

// dune-common
#include <dune/common/shared_ptr.hh>

// dune-stuff
#include <dune/stuff/function/expression.hh>

// dune-detailed-discretizations
#include <dune/detailed/discretizations/discretefunction/default.hh>

namespace Dune {

namespace Detailed {

namespace Discretizations {

namespace DiscreteFunctionSpace {

namespace Sub {

namespace Affine {

template <class BaseSpaceImp, class AffineShiftImp>
class Dirichlet
{
public:
  typedef BaseSpaceImp BaseSpaceType;

  typedef AffineShiftImp AffineShiftType;

  typedef Dirichlet<BaseSpaceType, AffineShiftType> ThisType;

  typedef typename BaseSpaceType::SuperSpaceType SuperSpaceType;

  typedef typename BaseSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef typename BaseSpaceType::GridPartType GridPartType;

  typedef typename BaseSpaceType::GridViewType GridViewType;

  enum
  {
    polynomialOrder = BaseSpaceType::polynomialOrder
  };

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

  typedef typename BaseSpaceType::PatternType PatternType;

  Dirichlet(const BaseSpaceType& baseSpace, const Dune::shared_ptr<const AffineShiftType> affineShift)
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

  const Dune::shared_ptr<const AffineShiftType> affineShift() const
  {
    return affineShift_;
  }

  const GridPartType& gridPart() const
  {
    return baseSpace_.gridPart();
  }

  const GridViewType& gridView() const
  {
    return baseSpace_.gridView();
  }

  const BaseFunctionSetType& baseFunctionSet() const
  {
    return baseSpace_.baseFunctionSet();
  }

  unsigned int order() const
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

  template <class LocalGridPartType, class OtherDiscreteFunctionSpaceType>
  Dune::shared_ptr<const PatternType> computeLocalPattern(const LocalGridPartType& localGridPart,
                                                          const OtherDiscreteFunctionSpaceType& other) const
  {
    return baseSpace_.computeLocalPattern(localGridPart, other);
  }

  template <class LocalGridPartType>
  Dune::shared_ptr<const PatternType> computeLocalPattern(const LocalGridPartType& localGridPart) const
  {
    return baseSpace_.computeLocalPattern(localGridPart);
  }

  template <class CouplingGridPartType, class OutsideDiscreteFunctionSpaceType>
  Dune::shared_ptr<const PatternType> computeCouplingPattern(const CouplingGridPartType& couplingGridPart,
                                                             const OutsideDiscreteFunctionSpaceType& outerSpace) const
  {
    return baseSpace_.computeCouplingPattern(couplingGridPart, outerSpace);
  }

  template <class OtherDiscreteFunctionSpaceType>
  Dune::shared_ptr<const PatternType> computePattern(const OtherDiscreteFunctionSpaceType& other) const
  {
    return baseSpace_.computePattern(other);
  }

  Dune::shared_ptr<const PatternType> computePattern() const
  {
    return baseSpace_.computePattern();
  }

private:
  Dirichlet(const ThisType&);
  ThisType& operator=(const ThisType&);

  const BaseSpaceType& baseSpace_;
  const Dune::shared_ptr<const AffineShiftType> affineShift_;
}; // end class Dirichlet

} // end namespace Affine

} // end namespace Sub

} // end namespace DiscreteFunctionSpace

} // namespace Discretizations

} // end namespace Detailed

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTIONSPACE_SUBSPACE_AFFINE_HH
