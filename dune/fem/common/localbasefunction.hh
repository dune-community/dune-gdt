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

  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

  typedef typename DiscreteFunctionSpaceType::IteratorType EntityIteratorType;

  typedef typename EntityIteratorType::Entity EntityType;

  typedef typename EntityType::Geometry EntityGeometryType;

private:
  class LocalBaseFunction
  {

  public:
    typedef LocalBaseFunction LocalBaseFunctionType;

    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;

    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;

    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef typename DiscreteFunctionSpaceType::HessianRangeType HessianRangeType;

    typedef typename EntityGeometryType::LocalCoordinate LocalCoordinateType;

    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;

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
      baseFunctionSet_.jacobian(localDoFNumber_, x, ret);
    }

    template <int diffOrder, class PointType>
    void evaluate(const FieldVector<int, diffOrder>& diffVariable, const PointType& x, RangeType& ret) const
    {
      baseFunctionSet_.evaluate(localDoFNumber_, diffVariable, x, ret);
    }

  private:
    const EntityType& entity_;
    const BaseFunctionSetType& baseFunctionSet_;
    const int localDoFNumber_;

  }; // end class LocalBaseFunction

public:
  typedef LocalBaseFunction LocalBaseFunctionType;

  LocalBaseFunctionProvider(const DiscreteFunctionSpaceType& discreteFunctionSpace)
    : discreteFunctionSpace_(discreteFunctionSpace)
  {
  }

  ~LocalBaseFunctionProvider()
  {
  }

  const LocalBaseFunctionType provide(const EntityType& entity, const int localDoFNumber) const
  {

    return LocalBaseFunction(entity, discreteFunctionSpace_.baseFunctionSet(entity), localDoFNumber);
  }

private:
  const DiscreteFunctionSpaceType& discreteFunctionSpace_;

}; // end class LocalBaseFunctionProvider

} // end namespace Common

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FEM_FUNCTIONALS_COMMON_LOCALBASEFUNCTIONSET_HH
