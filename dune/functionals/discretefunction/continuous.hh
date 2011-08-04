#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTION_CONTINUOUS_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTION_CONTINUOUS_HH

// dune-common includes
#include <dune/common/fvector.hh>

// dune-istl includes
#include <dune/istl/bvector.hh>

// local includes
#include "local.hh"

namespace Dune
{

namespace Functionals
{

namespace DiscreteFunction
{

namespace Continuous
{

template< class ContinuousDiscreteFunctionSpaceImp >
class BlockVector
{
public:

  typedef ContinuousDiscreteFunctionSpaceImp
    DiscreteFunctionSpaceType;

  typedef BlockVector< DiscreteFunctionSpaceType >
    ThisType;

  typedef Dune::Functionals::DiscreteFunction::Local< ThisType >
    LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::EntityType
    EntityType;

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

  typedef Dune::BlockVector< Dune::FieldVector< RangeFieldType, 1 > >
    StorageType;

  BlockVector( const DiscreteFunctionSpaceType& discreteFunctionSpace, const std::string name = "continuousBlockVectorFunction" )
    : space_( discreteFunctionSpace ),
      storage_( discreteFunctionSpace.size() ),
      name_( name )
  {

  }

  const DiscreteFunctionSpaceType& space() const
  {
    return space_;
  }

  const std::string name() const
  {
    return name_;
  }

  void rename( const std::string& newName = "" )
  {
    name_ = newName;
  }

  void clear()
  {
    std::cout << "DiscreteFunction::Continuous::BlockVector::clear() does not do anything!" << std::endl;
  }

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
    return storage_[globalDofNumber][0];
  }

  const RangeFieldType& operator[] ( const unsigned int globalDofNumber ) const
  {
    return storage_[globalDofNumber][0];
  }

  LocalFunctionType localFunction( const EntityType& entity )
  {
    return LocalFunctionType( (*this), entity );
  }

  const LocalFunctionType localFunction( const EntityType& entity ) const
  {
    return LocalFunctionType( (*this), entity );
  }

private:

  const DiscreteFunctionSpaceType& space_;
  StorageType storage_;
  std::string name_;

}; // end class BlockVector

} // end namespace Continuous

} // end namespace DiscreteFunction

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_DISCRETEFUNCTION_CONTINUOUS_HH
