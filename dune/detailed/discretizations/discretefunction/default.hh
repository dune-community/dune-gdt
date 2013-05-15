#ifndef DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTION_DEFAULT_HH
#define DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTION_DEFAULT_HH

#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/io/file/vtk/function.hh>

#include <dune/stuff/la/container/interface.hh>

#include <dune/detailed/discretizations/space/interface.hh>

#include "local.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {


// template< class DiscreteFunctionSpaceImp, class VectorImp >
// class DefaultConst;


// template< class DiscreteFunctionSpaceImp, class VectorImp >
// class Default
//  : public Dune::VTKFunction< typename DiscreteFunctionSpaceImp::GridViewType >
//{
// private:
//  typedef Dune::VTKFunction< typename DiscreteFunctionSpaceImp::GridViewType > BaseType;

// public:
//  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

//  typedef VectorImp VectorType;

//  typedef Default< DiscreteFunctionSpaceType, VectorType > ThisType;

//  typedef DefaultConst< DiscreteFunctionSpaceType, VectorType > ConstType;

//  typedef typename DiscreteFunctionSpaceType::GridViewType::template Codim< 0 >::Iterator::Entity EntityType;

//  static const int polOrder = DiscreteFunctionSpaceType::polynomialOrder;

//  typedef Dune::Detailed::Discretizations::DiscreteFunction::LocalConst< ThisType > ConstLocalFunctionType;

//  typedef Dune::Detailed::Discretizations::DiscreteFunction::Local< ThisType > LocalFunctionType;

//  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

//  typedef typename FunctionSpaceType::DomainType DomainType;

//  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

//  typedef typename FunctionSpaceType::RangeType RangeType;

//  static const int dimDomain = DiscreteFunctionSpaceType::dimDomain;

//  static const int dimRange = DiscreteFunctionSpaceType::dimRange;

//  typedef typename VectorType::size_type size_type;

//  Default(const DiscreteFunctionSpaceType& _space,
//          std::shared_ptr< VectorType > _vector,
//          const std::string _name = "discrete_function")
//    : BaseType(),
//      space_(_space)
//    , vector_(_vector)
//    , name_(_name)
//  {
//    assert(vector_->size() == space_.map().size() && "Given vector has wrong size!");
//  }

//  Default(const DiscreteFunctionSpaceType& _space,
//          const std::string _name = "discrete_function")
//    : BaseType(),
//      space_(_space)
//    , vector_(new VectorType(space_.map().size()))
//    , name_(_name)
//  {
//    assert(vector_->size() == space_.map().size() && "Given vector has wrong size!");
//  }

//  std::shared_ptr< const ConstType > createConst() const
//  {
//    std::shared_ptr< const ConstType > ret(new ConstType(*this));
//    return ret;
//  }

// private:
//  Default(const ThisType&);
//  ThisType& operator=(const ThisType&);

// public:
//  const DiscreteFunctionSpaceType& space() const
//  {
//    return space_;
//  }

//  virtual std::string name() const
//  {
//    return name_;
//  }

//  ConstLocalFunctionType localFunction(const EntityType& entity) const
//  {
//    return ConstLocalFunctionType(*this, entity);
//  }

//  LocalFunctionType localFunction(const EntityType& entity)
//  {
//    return LocalFunctionType(*this, entity);
//  }

//  std::shared_ptr< VectorType > vector()
//  {
//    return vector_;
//  }

//  const std::shared_ptr< const VectorType > vector() const
//  {
//    return vector_;
//  }

//  /**
//      @name Convenience methods
//      @{
//   **/
//  size_type size() const
//  {
//    return vector_->size();
//  }
//  /**
//      @}
//   **/

//  void clear()
//  {
//    const RangeFieldType zero(0);
//    for (size_type ii = 0; ii < vector_->size(); ++ii)
//      vector_->set(ii, zero);
//  }

//  /**
//      @name Methods to comply with the Dune::VTKFunction interface
//      @{
//   **/
//  virtual int ncomps() const
//  {
//    return dimRange;
//  }

//  virtual double evaluate(int component, const EntityType& entity, const DomainType& x) const
//  {
//    assert(vector_->size() == space_.map().size() && "Given vector has wrong size!");
//    RangeType ret(0.0);
//    localFunction(entity).evaluate(x, ret);
//    return ret[component];
//  }
//  /**
//      @}
//     **/

// private:
//  const DiscreteFunctionSpaceType& space_;
//  std::shared_ptr< VectorType > vector_;
//  const std::string name_;
//}; // class Defaul


template <class SpaceImp, class VectorImp>
class DiscreteFunctionDefaultConst : public Dune::VTKFunction<typename SpaceImp::GridPartType::GridViewType>,
                                     public Dune::Stuff::LocalizableFunction
{
public:
  typedef typename SpaceInterface<typename SpaceImp::Traits>::derived_type SpaceType;
  typedef typename Dune::Stuff::LA::Container::VectorInterface<typename VectorImp::Traits>::derived_type VectorType;

  typedef typename SpaceType::EntityType EntityType;

  typedef DiscreteFunctionLocalConst<SpaceType, VectorType> ConstLocalFunctionType;
  typedef typename ConstLocalFunctionType::DomainType DomainType;
  typedef typename ConstLocalFunctionType::RangeType RangeType;

  DiscreteFunctionDefaultConst(const SpaceType& sp, const std::shared_ptr<const VectorType> vec,
                               const std::string nm = "discrete_function")
    : /*BaseType()
    , */ space_(sp)
    , name_(nm)
    , vector_(vec)
  {
    assert(vector_->size() == space_.mapper().size() && "Given vector has wrong size!");
  }

  const SpaceType& space() const
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

  const std::shared_ptr<const VectorType> vector() const
  {
    return vector_;
  }

  /** \defgroup vtk ´´Methods to comply with the Dune::VTKFunction interface.'' */
  /* @{ */
  virtual int ncomps() const
  {
    return SpaceType::dimRangeRows;
  }

  virtual double evaluate(int component, const EntityType& entity, const DomainType& x) const
  {
    RangeType ret(0);
    localFunction(entity).evaluate(x, ret);
    return ret[component];
  }
  /* @} */

private:
  const SpaceType& space_;
  const std::string name_;
  const std::shared_ptr<const VectorType> vector_;
}; // class DiscreteFunctionDefaultConst


template <class SpaceImp, class VectorImp>
class DiscreteFunctionDefault : public DiscreteFunctionDefaultConst<SpaceImp, VectorImp>
{
  typedef DiscreteFunctionDefaultConst<SpaceImp, VectorImp> BaseType;

public:
  typedef typename SpaceInterface<typename SpaceImp::Traits>::derived_type SpaceType;
  typedef typename Dune::Stuff::LA::Container::VectorInterface<typename VectorImp::Traits>::derived_type VectorType;

  typedef typename SpaceType::EntityType EntityType;

  typedef DiscreteFunctionLocal<SpaceType, VectorType> LocalFunctionType;
  typedef typename LocalFunctionType::DomainType DomainType;

  DiscreteFunctionDefault(const SpaceType& sp, std::shared_ptr<VectorType> vec,
                          const std::string nm = "discrete_function")
    : BaseType(sp, vec, nm)
    , nonConstVector_(vec)
  {
  }

  //  DiscreteFunctionDefault(const SpaceType& sp,
  //                          const std::string nm = "discrete_function")
  //    : nonConstVector_(std::make_shared< VectorType >(sp.mapper().size()))
  //    , BaseType(sp, nonConstVector_, nm)
  //  {}

  LocalFunctionType localFunction(const EntityType& entity)
  {
    return LocalFunctionType(*this, entity);
  }

  std::shared_ptr<VectorType> vector()
  {
    return nonConstVector_;
  }

  //  virtual std::string name() const
  //  {
  //    return BaseType::name();
  //  }

  //  /** \defgroup vtk ´´Methods to comply with the Dune::VTKFunction interface.'' */
  //  /* @{ */
  //  virtual int ncomps() const
  //  {
  //    return BaseType::ncomps();
  //  }

  //  virtual double evaluate(int component, const EntityType& entity, const DomainType& x) const
  //  {
  //    BaseType::evaluate(component, entity, x);
  //  }
  //  /* @} */

private:
  std::shared_ptr<VectorType> nonConstVector_;
}; // class DiscreteFunctionDefault


} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTION_DEFAULT_HH
