#ifndef DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTION_DISCONTINUOUS_HH
#define DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTION_DISCONTINUOUS_HH

// dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

// dune-istl includes
#include <dune/istl/bvector.hh>

// dune-helper-tools includes
#include <dune/helper-tools/discretefunctionspace/projection/dirichlet.hh>

// local includes
#include "local.hh"

namespace Dune {

namespace DetailedDiscretizations {

namespace DiscreteFunction {

namespace Discontinuous {

template <class ContinuousDiscreteFunctionSpaceImp>
class BlockVector
{
public:
  typedef ContinuousDiscreteFunctionSpaceImp DiscreteFunctionSpaceType;

  enum
  {
    polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder
  };

  typedef BlockVector<DiscreteFunctionSpaceType> ThisType;

  typedef Dune::DetailedDiscretizations::DiscreteFunction::Local<ThisType> LocalFunctionType;

  typedef Dune::DetailedDiscretizations::DiscreteFunction::LocalConst<ThisType> ConstLocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;

  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;

  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

  typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

  static const int dimDomain = DiscreteFunctionSpaceType::dimDomain;

  static const int dimRange = DiscreteFunctionSpaceType::dimRange;

  typedef Dune::BlockVector<Dune::FieldVector<RangeFieldType, dimRange>> StorageType;

  BlockVector(const DiscreteFunctionSpaceType& discreteFunctionSpace,
              const std::string name = "continuousBlockVectorFunction")
    : space_(discreteFunctionSpace)
    , storage_(discreteFunctionSpace.map().size())
    , name_(name)
  {
    clear();
  }

  BlockVector(const DiscreteFunctionSpaceType& discreteFunctionSpace, const StorageType& storage,
              const std::string name = "continuousBlockVectorFunction")
    : space_(discreteFunctionSpace)
    , storage_(space_.map().size())
    , name_(name)
  {
    assert(storage.size() == storage_.size());
    for (unsigned int i = 0; i < storage_.size(); ++i) {
      storage_[i] = storage[i];
    }
  }

  template <class FunctionType>
  BlockVector(const DiscreteFunctionSpaceType& discreteFunctionSpace, const std::string name,
              const FunctionType function, const std::string projectionType)
    : space_(discreteFunctionSpace)
    , storage_(space_.map().size())
    , name_(name)
  {
    if (projectionType.compare("dirichlet") == 0) {
      Dune::HelperTools::Projection::Dirichlet::project(function, *this);
    } else {
      throw Dune::NotImplemented();
    }
  }

  //! copy constructor
  BlockVector(const ThisType& other)
    : space_(other.space())
    , storage_(space_.map().size())
    , name_("copyOF" + other.name())
  {
    for (unsigned int i = 0; i < storage_.size(); ++i) {
      operator[](i) = other[i];
    }
  }

private:
  //! assignment operator
  ThisType& operator=(const ThisType& other)
  {
    if (this != other) {
      assert(other.space().map().size() == this->space().map().size());
      for (unsigned int i = 0; i < storage_.size(); ++i) {
        operator[](i) = other[i];
      }
    }
    return *this;
  }

public:
  const DiscreteFunctionSpaceType& space() const
  {
    return space_;
  }

  const std::string name() const
  {
    return name_;
  }

  void setName(const std::string& newName = "")
  {
    name_ = newName;
  }

  void clear()
  {
    storage_ = 0.0;
  }

  const StorageType& storage() const
  {
    return storage_;
  }

  StorageType& storage()
  {
    return storage_;
  }

  RangeFieldType& operator[](const unsigned int globalDofNumber)
  {
    return storage_[globalDofNumber][0];
  }

  const RangeFieldType& operator[](const unsigned int globalDofNumber) const
  {
    return storage_[globalDofNumber][0];
  }

  template <class EntityType>
  LocalFunctionType localFunction(const EntityType& entity)
  {
    return LocalFunctionType((*this), entity);
  }

  template <class EntityType>
  ConstLocalFunctionType localFunction(const EntityType& entity) const
  {
    return ConstLocalFunctionType((*this), entity);
  }

  /**
    \attention  This is not correct for order 0
    \todo       fix me
    **/
  bool continuous() const
  {
    return true;
  }

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

private:
  const DiscreteFunctionSpaceType& space_;
  StorageType storage_;
  std::string name_;

}; // end class BlockVector

} // end namespace Discontinuous

} // end namespace DiscreteFunction

} // end namespace DetailedDiscretizations

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_DISCRETEFUNCTION_DISCONTINUOUS_HH
