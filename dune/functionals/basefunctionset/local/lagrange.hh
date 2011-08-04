#ifndef DUNE_FUNCTIONALS_BASEFUNCTIONSET_LOCAL_LAGRANGE_HH
#define DUNE_FUNCTIONALS_BASEFUNCTIONSET_LOCAL_LAGRANGE_HH

// dune-fem includes
//#include <dune/fem/space/lagrangespace/lagrangespace.hh>

// dune-functionals includes
//#include <dune/functionals/discretefunctionspace/continuous/lagrangefemadapter.hh>

// dune-fem-tools includes
#include <dune/fem-tools/common/printing.hh>
#include <dune/fem-tools/common/string.hh>
#include <dune/fem-tools/grid/entity.hh>

namespace Dune {

namespace Functionals {

namespace BaseFunctionSet {

namespace Local {

template <class BaseFunctionSetImp>
class Lagrange
{
public:
  typedef BaseFunctionSetImp BaseFunctionSetType;

  typedef typename BaseFunctionSetType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

  typedef typename DiscreteFunctionSpaceType::EntityType EntityType;

  enum
  {
    polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder
  };

  typedef Lagrange<BaseFunctionSetType> ThisType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

  Lagrange(const BaseFunctionSetType& baseFunctionSet, const EntityType& entity)
    : baseFunctionSet_(baseFunctionSet)
    , entity_(entity)
  {
    // get the host basefunctioset
    typedef typename BaseFunctionSetType::BaseFunctionSetType HostBaseFunctionSetType;
    const HostBaseFunctionSetType tmpBaseFunctionSet = baseFunctionSet_.baseFunctionSet(entity_);

    size_ = tmpBaseFunctionSet.numBaseFunctions();
    // this is still fishy
    order_ = baseFunctionSet_.space().order();
  }

  const BaseFunctionSetType& baseFunctionSet() const
  {
    return baseFunctionSet_;
  }

  const EntityType& entity() const
  {
    return entity_;
  }

  unsigned int size() const
  {
    return size_;
  }

  int order() const
  {
    return order_;
  }

  void evaluate(const DomainType& x, std::vector<RangeType>& ret) const
  {
    // get the host basefunctioset
    typedef typename BaseFunctionSetType::BaseFunctionSetType HostBaseFunctionSetType;
    const HostBaseFunctionSetType baseFunctionSet = baseFunctionSet_.baseFunctionSet(entity_);

    // and evaluate
    for (unsigned int i = 0; i < size_; ++i) {
      baseFunctionSet.evaluate(i, x, ret[i]);
    }
  }

  void jacobian(const DomainType& x, std::vector<JacobianRangeType>& ret) const
  {
    // some types we will need
    typedef typename EntityType::Geometry EntityGeometryType;
    typedef typename EntityGeometryType::Jacobian JacobianInverseTransposedType;
    typedef typename JacobianRangeType::row_type JacobianRowType;

    // get the host basefunctioset
    typedef typename BaseFunctionSetType::BaseFunctionSetType HostBaseFunctionSetType;
    const HostBaseFunctionSetType baseFunctionSet = baseFunctionSet_.baseFunctionSet(entity_);

    // geometry and jacobian inverse transposed
    const EntityGeometryType& entityGeometry                       = entity_.geometry();
    const JacobianInverseTransposedType& jacobianInverseTransposed = entityGeometry.jacobianInverseTransposed(x);

    // some tmp storage
    JacobianRangeType jacobianUntransposed(0.0);
    JacobianRangeType jacobianTransposed(0.0);

    // evaluate
    for (unsigned int i = 0; i < size_; ++i) {
      // get untransposed jacobian
      baseFunctionSet.jacobian(i, x, jacobianUntransposed);

      // transpose for each dim of range
      const unsigned int dimRange = DiscreteFunctionSpaceType::dimRange;
      for (unsigned int row = 0; row < dimRange; ++row) {
        // transpose
        jacobianInverseTransposed.mv(jacobianUntransposed[row], jacobianTransposed[row]);
        ret[i][row] = jacobianTransposed[row];
      }
    }
  }

private:
  //  friend class Dune::Functionals::DiscreteFunctionSpace::Continuous::LagrangeFemAdapter< DiscreteFunctionSpaceType
  //  >;

  const BaseFunctionSetType& baseFunctionSet_;
  const EntityType& entity_;
  unsigned int size_;
  int order_;

}; // end class Lagrange

} // end namespace Local

} // end namespace Common

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_BASEFUNCTIONSET_LOCAL_LAGRANGE_HH
