#ifndef DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTION_CONTINUOUS_HH
#define DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTION_CONTINUOUS_HH

// dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

// dune-grid
#include <dune/grid/io/file/vtk/function.hh>

// local includes
#include "local.hh"

namespace Dune
{

namespace Detailed {

namespace Discretizations
{

namespace DiscreteFunction
{

template< class DiscreteFunctionSpaceImp, class VectorBackendImp >
class Default
  : public Dune::VTKFunction< typename DiscreteFunctionSpaceImp::GridViewType >
{
public:
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  typedef VectorBackendImp StorageType;

  typedef Default< DiscreteFunctionSpaceType, StorageType >
    ThisType;

private:
  typedef Dune::VTKFunction< typename DiscreteFunctionSpaceImp::GridViewType > BaseType;

public:
  typedef typename DiscreteFunctionSpaceType::GridViewType::template Codim< 0 >::Iterator::Entity EntityType;

  enum{ polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder };

  typedef Dune::Detailed::Discretizations::DiscreteFunction::Local< ThisType >
    LocalFunctionType;

  typedef Dune::Detailed::Discretizations::DiscreteFunction::LocalConst< ThisType >
    ConstLocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType
    FunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::DomainFieldType
    DomainFieldType;

  typedef typename DiscreteFunctionSpaceType::DomainType
    DomainType;

  typedef typename DiscreteFunctionSpaceType::RangeFieldType
    RangeFieldType;

  typedef typename DiscreteFunctionSpaceType::RangeType
    RangeType;

  typedef typename DiscreteFunctionSpaceType::JacobianRangeType
    JacobianRangeType;

  static const int dimDomain = DiscreteFunctionSpaceType::dimDomain;

  static const int dimRange = DiscreteFunctionSpaceType::dimRange;

  Default(const DiscreteFunctionSpaceType& space,
          const std::string name = "default" )
    : BaseType(),
      space_(space),
      storage_(space_.map().size()),
      name_(name)
  {
//    clear();
  }

  Default(const DiscreteFunctionSpaceType& space,
          StorageType storage,
          const std::string name = "default")
    : BaseType(),
      space_(space),
      storage_(storage),
      name_(name)
  {
  }

//  template< class FunctionType >
//  BlockVector(  const DiscreteFunctionSpaceType& discreteFunctionSpace,
//                const std::string name,
//                const FunctionType function,
//                const std::string projectionType )
//    : space_( discreteFunctionSpace ),
//      storage_( space_.map().size() ),
//      name_( name )
//  {
//    if( projectionType.compare( "dirichlet" ) == 0 )
//    {
//      Dune::HelperTools::Projection::Dirichlet::project( function, *this );
//    }
//    else
//    {
//      throw Dune::NotImplemented();
//    }
//  }

private:
  //! copy constructor
  Default(const ThisType& other)
    : BaseType(),
      space_(other.space()),
      storage_(space_.map().size()),
      name_("copyOF" + other.name())
  {
    for (unsigned int i = 0; i < storage_.size(); ++i)
    {
      storage_[i] = storage[i];
    }
  }

  //! assignment operator
  ThisType& operator=( const ThisType& other )
  {
    if(this != other)
    {
      assert(other.space().map().size() == this->space().map().size());
      for (unsigned int i = 0; i < storage_.size(); ++i)
      {
        storage_[i] = storage[i];
      }
    }
    return *this;
  }

public:
  const DiscreteFunctionSpaceType& space() const
  {
    return space_;
  }

  virtual std::string name() const
  {
    return name_;
  }

  void setName( const std::string& newName = "" )
  {
    name_ = newName;
  }

//  void clear()
//  {
//    storage_ = 0.0;
//  }

  const StorageType& storage() const
  {
    return storage_;
  }

  StorageType& storage()
  {
    return storage_;
  }

  RangeFieldType& operator[] ( const unsigned int globalDofNumber )
  {
    return storage_[globalDofNumber];
  }

  const RangeFieldType& operator[] ( const unsigned int globalDofNumber ) const
  {
    return storage_[globalDofNumber];
  }

//  template< class EntityType >
  LocalFunctionType localFunction( const EntityType& entity )
  {
    return LocalFunctionType( (*this), entity );
  }

//  template< class EntityType >
  ConstLocalFunctionType localFunction( const EntityType& entity ) const
  {
    return ConstLocalFunctionType( (*this), entity );
  }

//  /**
//    \attention  This is not correct for order 0
//    \todo       fix me
//    **/
//  bool continuous() const
//  {
//    return true;
//  }

  /**
      @name Convenience methods
      @{
   **/

  /**
    \attention  someone should think about this at some point (i.e. h-adaptivity)
    **/
  int order() const
  {
    return space_.order();
  }

  unsigned int size() const
  {
    return space_.map().size();
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

  virtual RangeFieldType evaluate(int component, const EntityType& entity, const DomainType& x) const
  {
    RangeType ret(0.0);
    localFunction(entity).evaluate(x, ret);
    return ret[component];
  }

  /**
      @}
     **/
private:

  const DiscreteFunctionSpaceType& space_;
  StorageType storage_;
  std::string name_;

}; // end class Default

} // end namespace DiscreteFunction

} // namespace Discretizations

} // namespace Detailed

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTION_CONTINUOUS_HH
