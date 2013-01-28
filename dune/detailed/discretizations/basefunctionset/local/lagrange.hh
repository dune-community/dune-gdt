#ifndef DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_LOCAL_LAGRANGE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_LOCAL_LAGRANGE_HH

// system
#include <vector>

// dune-common
#include <dune/common/shared_ptr.hh>
#include <dune/common/timer.hh>

namespace Dune {

namespace Detailed {

namespace Discretizations {

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

  typedef typename GridPartType::template Codim<0>::EntityType EntityType;

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

  //! constructor
  Lagrange(const BaseFunctionSetType& baseFunctionSet, const EntityType& entity)
    : baseFunctionSet_(baseFunctionSet)
    , entity_(entity)
    , size_(0)
    , order_(0)
  {
    // get the host basefunctioset
    typedef typename BaseFunctionSetType::BaseFunctionSetType HostBaseFunctionSetType;
    const HostBaseFunctionSetType tmpBaseFunctionSet = baseFunctionSet_.baseFunctionSet(entity_);

    size_ = tmpBaseFunctionSet.size();
    // this is still fishy, i.e. p-adaptivity
    order_ = baseFunctionSet_.space().order();
  }

  //! copy constructor
  Lagrange(const ThisType& other)
    : baseFunctionSet_(other.baseFunctionSet())
    , entity_(other.entity())
    , size_(other.size())
    , order_(other.order())
  {
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
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
      baseFunctionSet.evaluate(i, x, ret[i]);
#pragma GCC diagnostic pop
    }
  }

  void jacobian(const DomainType& x, std::vector<JacobianRangeType>& ret) const
  {
    //    typedef typename BaseFunctionSetType::JacobianCacheMapType JacobianCacheMapType;
    //    const JacobianCacheMapType& constJacobianCacheMap = baseFunctionSet_.jacobianCacheMap();
    //    JacobianCacheMapType& jacobianCacheMap = const_cast< JacobianCacheMapType& >(constJacobianCacheMap);
    //    if (jacobianCacheMap.find(x) != jacobianCacheMap.end()) {
    //        ret = jacobianCacheMap[x];
    //    } else {
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
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
      baseFunctionSet.jacobian(i, x, jacobianUntransposed);
#pragma GCC diagnostic pop

      // transpose for each dim of range
      const unsigned int dimRange = DiscreteFunctionSpaceType::dimRange;
      for (unsigned int row = 0; row < dimRange; ++row) {
        // transpose
        jacobianInverseTransposed.mv(jacobianUntransposed[row], jacobianTransposed[row]);
        ret[i][row] = jacobianTransposed[row];
      }
    }
    //      jacobianCacheMap[x] = ret;
    //    }
  }

private:
  //! assignment operator
  ThisType& operator=(const ThisType&);

  const BaseFunctionSetType& baseFunctionSet_;
  const EntityType& entity_;
  unsigned int size_;
  int order_;

}; // end class Lagrange

} // end namespace Local

} // end namespace Common

} // namespace Discretizations

} // namespace Detailed

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_LOCAL_LAGRANGE_HH
