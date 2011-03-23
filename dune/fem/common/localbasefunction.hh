#ifndef DUNE_FEM_FUNCTIONALS_COMMON_LOCALBASEFUNCTIONSET_HH
#define DUNE_FEM_FUNCTIONALS_COMMON_LOCALBASEFUNCTIONSET_HH

// dune-fem includes
#include <dune/fem/function/localfunction/localfunction.hh>


namespace Dune {

namespace Functionals {

namespace Common {

template <class DiscreteFunctionSpaceImp>
class LocalBaseFunctionProvider
{
public:
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  class LocalBaseFunction
  {

  public:
    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;

    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;

    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef typename DiscreteFunctionSpaceType::HessianRangeType HessianRangeType;

    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;

    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

    typedef typename DiscreteFunctionSpaceType::IteratorType EntityIteratorType;

    typedef typename EntityIteratorType::Entity EntityType;

    typedef typename EntityType::Geometry EntityGeometryType;

    static const int dimDomain = DiscreteFunctionSpaceType::dimDomain;

    static const int dimRange = DiscreteFunctionSpaceType::dimRange;

    LocalBaseFunction(const EntityType& entity, const BaseFunctionSetType& baseFunctionSet, const int localDoFNumber)
      : entity_(entity)
      , baseFunctionSet_(baseFunctionSet)
      , localDoFNumber_(localDoFNumber)
    {
    }

    ~LocalBaseFunction()
    {
    }

    template <class PointType>
    void evaluate(const PointType& x, RangeType& ret) const
    {
      baseFunctionSet_.evaluate(localDoFNumber_, x, ret);
    }

    template <class PointType>
    void jacobian(const PointType& x, JacobianRangeType& ret) const
    {
      // some types we will need
      typedef typename EntityGeometryType::Jacobian JacobianInverseTransposedType;

      typedef typename JacobianRangeType::row_type JacobianRowType;

      // geometry and jacobian inverse transposed
      const EntityGeometryType& entityGeometry                       = entity_.geometry();
      const JacobianInverseTransposedType& jacobianInverseTransposed = entityGeometry.jacobianInverseTransposed(x);

      // get untransposed jacobian
      JacobianRangeType jacobianUntransposed(0.0);
      baseFunctionSet_.jacobian(localDoFNumber_, x, jacobianUntransposed);

      // do for each dim of Range
      for (unsigned int row = 0; row < ret.N(); ++row) {
        // transpose
        JacobianRowType jacobian(0.0);
        jacobianInverseTransposed.mv(jacobianUntransposed[0], jacobian);

        // return
        ret[row] = jacobian;
      }
    }

    const EntityType& entity() const
    {
      return entity_;
    }

  private:
    const EntityType& entity_;
    const BaseFunctionSetType baseFunctionSet_;
    const int localDoFNumber_;

  }; // end class LocalBaseFunction

  typedef LocalBaseFunction LocalBaseFunctionType;

  typedef typename LocalBaseFunctionType::EntityType EntityType;

  LocalBaseFunctionProvider(const DiscreteFunctionSpaceType& discreteFunctionSpace)
    : discreteFunctionSpace_(discreteFunctionSpace)
  {
  }

  ~LocalBaseFunctionProvider()
  {
  }

  const LocalBaseFunctionType provide(const EntityType& entity, const int localDoFNumber) const
  {
    return LocalBaseFunctionType(entity, discreteFunctionSpace_.baseFunctionSet(entity), localDoFNumber);
  }

private:
  const DiscreteFunctionSpaceType& discreteFunctionSpace_;

}; // end class LocalBaseFunctionProvider

} // end namespace Common

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FEM_FUNCTIONALS_COMMON_LOCALBASEFUNCTIONSET_HH
