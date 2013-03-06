#ifndef DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTION_DEFAULT_HH
#define DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTION_DEFAULT_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/grid/io/file/vtk/function.hh>

#include "local.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace DiscreteFunction {

template< class DiscreteFunctionSpaceImp, class VectorImp >
class DefaultConst;

template< class DiscreteFunctionSpaceImp, class VectorImp >
class Default
  : public Dune::VTKFunction< typename DiscreteFunctionSpaceImp::GridViewType >
{
private:
  typedef Dune::VTKFunction< typename DiscreteFunctionSpaceImp::GridViewType > BaseType;

public:
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  typedef VectorImp VectorType;

  typedef Default< DiscreteFunctionSpaceType, VectorType > ThisType;

  typedef DefaultConst< DiscreteFunctionSpaceType, VectorType > ConstType;

  typedef typename DiscreteFunctionSpaceType::GridViewType::template Codim< 0 >::Iterator::Entity EntityType;

  static const int polOrder = DiscreteFunctionSpaceType::polynomialOrder;

  typedef Dune::Detailed::Discretizations::DiscreteFunction::LocalConst< ThisType > ConstLocalFunctionType;

  typedef Dune::Detailed::Discretizations::DiscreteFunction::Local< ThisType > LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  static const int dimDomain = DiscreteFunctionSpaceType::dimDomain;

  static const int dimRange = DiscreteFunctionSpaceType::dimRange;

  typedef typename VectorType::size_type size_type;

  Default(const DiscreteFunctionSpaceType& _space,
          Dune::shared_ptr< VectorType > _vector,
          const std::string _name = "discrete_function")
    : BaseType(),
      space_(_space)
    , vector_(_vector)
    , name_(_name)
  {
    assert(vector_->size() == space_.map().size() && "Given vector has wrong size!");
  }

  Default(const DiscreteFunctionSpaceType& _space,
          const std::string _name = "discrete_function")
    : BaseType(),
      space_(_space)
    , vector_(new VectorType(space_.map().size()))
    , name_(_name)
  {
    assert(vector_->size() == space_.map().size() && "Given vector has wrong size!");
  }

  Dune::shared_ptr< const ConstType > createConst() const
  {
    Dune::shared_ptr< const ConstType > ret(new ConstType(*this));
    return ret;
  }

private:
  Default(const ThisType&);
  ThisType& operator=(const ThisType&);

public:
  const DiscreteFunctionSpaceType& space() const
  {
    return space_;
  }

  virtual std::string name() const
  {
    return name_;
  }

  ConstLocalFunctionType localFunction(const EntityType& entity) const
  {
    return ConstLocalFunctionType(*this, entity);
  }

  LocalFunctionType localFunction(const EntityType& entity)
  {
    return LocalFunctionType(*this, entity);
  }

  Dune::shared_ptr< VectorType > vector()
  {
    return vector_;
  }

  const Dune::shared_ptr< const VectorType > vector() const
  {
    return vector_;
  }

  /**
      @name Convenience methods
      @{
   **/
  size_type size() const
  {
    return vector_->size();
  }
  /**
      @}
   **/

  void clear()
  {
    const RangeFieldType zero(0);
    for (size_type ii = 0; ii < vector_->size(); ++ii)
      vector_->set(ii, zero);
  }

  /**
      @name Methods to comply with the Dune::VTKFunction interface
      @{
   **/
  virtual int ncomps() const
  {
    return dimRange;
  }

  virtual double evaluate(int component, const EntityType& entity, const DomainType& x) const
  {
    assert(vector_->size() == space_.map().size() && "Given vector has wrong size!");
    RangeType ret(0.0);
    localFunction(entity).evaluate(x, ret);
    return ret[component];
  }
  /**
      @}
     **/

private:
  const DiscreteFunctionSpaceType& space_;
  Dune::shared_ptr< VectorType > vector_;
  const std::string name_;
}; // class Defaul


template< class DiscreteFunctionSpaceImp, class VectorImp >
class DefaultConst
  : public Dune::VTKFunction< typename DiscreteFunctionSpaceImp::GridViewType >
{
private:
  typedef Dune::VTKFunction< typename DiscreteFunctionSpaceImp::GridViewType > BaseType;

public:
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  typedef VectorImp VectorType;

  typedef DefaultConst< DiscreteFunctionSpaceType, VectorType > ThisType;

  typedef Default< DiscreteFunctionSpaceType, VectorType > NonConstType;

  typedef typename DiscreteFunctionSpaceType::GridViewType::template Codim< 0 >::Iterator::Entity EntityType;

  static const int polOrder = DiscreteFunctionSpaceType::polynomialOrder;

  typedef Dune::Detailed::Discretizations::DiscreteFunction::LocalConst< ThisType > ConstLocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  static const int dimDomain = DiscreteFunctionSpaceType::dimDomain;

  static const int dimRange = DiscreteFunctionSpaceType::dimRange;

  typedef typename VectorType::size_type size_type;

  DefaultConst(const DiscreteFunctionSpaceType& _space,
               const Dune::shared_ptr< const VectorType > _vector,
               const std::string _name = "discrete_function")
    : BaseType()
    , space_(_space)
    , name_(_name)
    , vector_(_vector)
  {
    assert(vector_->size() == space_.map().size() && "Given vector has wrong size!");
  }

  DefaultConst(const NonConstType& _nonConstOther)
    : BaseType()
    , space_(_nonConstOther.space())
    , name_(_nonConstOther.name())
    , vector_(_nonConstOther.vector())
  {}

private:
  DefaultConst(const ThisType&);
  ThisType& operator=(const ThisType&);

public:
  const DiscreteFunctionSpaceType& space() const
  {
    return space_;
  }

  virtual std::string name() const
  {
    return name_;
  }

  ConstLocalFunctionType localFunction(const EntityType& entity) const
  {
    return ConstLocalFunctionType(*this, entity);
  }

  const Dune::shared_ptr< const VectorType > vector() const
  {
    return vector_;
  }

  /**
      @name Convenience methods
      @{
   **/
  size_type size() const
  {
    return vector_->size();
  }
  /**
      @}
   **/

  /**
      @name Methods to comply with the Dune::VTKFunction interface
      @{
   **/
  virtual int ncomps() const
  {
    return dimRange;
  }

  virtual double evaluate(int component, const EntityType& entity, const DomainType& x) const
  {
    assert(vector_->size() == space_.map().size() && "Given vector has wrong size!");
    RangeType ret(0.0);
    localFunction(entity).evaluate(x, ret);
    return ret[component];
  }
  /**
      @}
     **/

private:
  const DiscreteFunctionSpaceType& space_;
  const std::string name_;
  const Dune::shared_ptr< const VectorType > vector_;
}; // class DefaulConst

} // namespace DiscreteFunction
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTION_DEFAULT_HH
