// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Kirsten Weber

#ifndef DUNE_GDT_DISCRETEFUNCTION_DEFAULT_HH
#define DUNE_GDT_DISCRETEFUNCTION_DEFAULT_HH

#include <memory>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/stuff/common/disable_warnings.hh>
# include <dune/common/fvector.hh>

# include <dune/grid/io/file/vtk.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/la/container/interfaces.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/common/memory.hh>

#include <dune/gdt/spaces/interface.hh>

#include "local.hh"

namespace Dune {
namespace GDT {


template< class SpaceImp, class VectorImp >
class ConstDiscreteFunction
  : public Stuff::LocalizableFunctionInterface< typename SpaceImp::EntityType,
                                                typename SpaceImp::DomainFieldType, SpaceImp::dimDomain,
                                                typename SpaceImp::RangeFieldType, SpaceImp::dimRange, SpaceImp::dimRangeCols >
{
  static_assert(std::is_base_of< SpaceInterface< typename SpaceImp::Traits >, SpaceImp >::value,
                "SpaceImp has to be derived from SpaceInterface!");
  static_assert(std::is_base_of< Dune::Stuff::LA::VectorInterface< typename VectorImp::Traits >, VectorImp >::value,
                "VectorImp has to be derived from Stuff::LA::VectorInterface!");
  static_assert(std::is_same< typename SpaceImp::RangeFieldType, typename VectorImp::ScalarType >::value,
                "Types do not match!");
  typedef Stuff::LocalizableFunctionInterface
      < typename SpaceImp::EntityType, typename SpaceImp::DomainFieldType, SpaceImp::dimDomain,
        typename SpaceImp::RangeFieldType, SpaceImp::dimRange, SpaceImp::dimRangeCols >
    BaseType;
  typedef ConstDiscreteFunction< SpaceImp, VectorImp > ThisType;
public:
  typedef SpaceImp                              SpaceType;
  typedef VectorImp                             VectorType;
  typedef typename BaseType::EntityType         EntityType;
  typedef typename BaseType::LocalfunctionType  LocalfunctionType;

  static const unsigned int                   dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainFieldType  DomainFieldType;
  typedef typename BaseType::DomainType       DomainType;

  static const unsigned int                     dimRange = BaseType::dimRange;
  static const unsigned int                     dimRangeCols = BaseType::dimRangeCols;
  typedef typename BaseType::RangeFieldType     RangeFieldType;
  typedef typename BaseType::RangeType          RangeType;
  typedef typename BaseType::JacobianRangeType  JacobianRangeType;

  typedef ConstLocalDiscreteFunction< SpaceType, VectorType > ConstLocalDiscreteFunctionType;

  ConstDiscreteFunction(const SpaceType& sp,
                        const VectorType& vec,
                        const std::string nm = "dune.gdt.constdiscretefunction")
    : space_(sp)
    , vector_(vec)
    , name_(nm)
  {
    assert(vector_.size() == space_.mapper().size() && "Given vector has wrong size!");
  }

  ConstDiscreteFunction(const ThisType& other)
    : space_(other.space_)
    , vector_(other.vector_)
    , name_(other.name_)
  {}

  ConstDiscreteFunction(ThisType&& source)
    : space_(source.space_)
    , vector_(source.vector_)
    , name_(std::move(source.name_))
  {}

  virtual ~ConstDiscreteFunction() {}

  ThisType& operator=(const ThisType& other) = delete;

  virtual ThisType* copy() const DS_OVERRIDE
  {
    return new ThisType(*this);
  }

  virtual std::string name() const DS_OVERRIDE
  {
    return name_;
  }

  const SpaceType& space() const
  {
    return space_;
  }

  const VectorType& vector() const
  {
    return vector_;
  }

  ConstLocalDiscreteFunctionType local_discrete_function(const EntityType& entity) const
  {
    assert(space_.grid_view()->indexSet().contains(entity));
    return ConstLocalDiscreteFunctionType(space_, vector_, entity);
  }

  virtual std::unique_ptr< LocalfunctionType > local_function(const EntityType& entity) const DS_OVERRIDE
  {
    return std::unique_ptr< ConstLocalDiscreteFunctionType >
        (new ConstLocalDiscreteFunctionType(local_discrete_function(entity)));
  }

  void visualize(const std::string filename,
                 const bool subsampling = (SpaceType::polOrder > 1),
                 VTK::OutputType vtk_output_type = VTK::appendedraw) const
  {
    BaseType::template visualize< typename SpaceType::GridViewType >(*(space().grid_view()),
                                                                     filename,
                                                                     subsampling,
                                                                     vtk_output_type);
  }

  bool dofs_valid() const
  {
    return vector().valid();
  }

protected:
  const SpaceType& space_;
private:
  const VectorType& vector_;
  const std::string name_;
}; // class ConstDiscreteFunction


namespace internal {


template< class VectorType >
class RefProvider
{
public:
  virtual ~RefProvider() {}

  virtual VectorType& ref() = 0;
};


template< class VectorType >
class RefProviderByRef
  : public RefProvider< VectorType >
{
public:
  RefProviderByRef(VectorType& vec)
    : vector_(vec)
  {}

  virtual ~RefProviderByRef() {}

  virtual VectorType& ref() DS_OVERRIDE DS_FINAL
  {
    return vector_;
  }

private:
  VectorType& vector_;
};


template< class VectorType >
class RefProviderByPtr
  : public RefProvider< VectorType >
{
public:
  RefProviderByPtr(const size_t sz)
    : vector_(new VectorType(sz))
  {}

  virtual ~RefProviderByPtr() {}

  virtual VectorType& ref() DS_OVERRIDE DS_FINAL
  {
    return *vector_;
  }

private:
  std::unique_ptr< VectorType > vector_;
};


template< class VectorType >
class VectorProvider
{
public:
  VectorProvider(VectorType& ref)
    : ref_provider_(Stuff::Common::make_unique< RefProviderByRef< VectorType > >(ref))
  {}

  VectorProvider(const size_t sz)
    : ref_provider_(Stuff::Common::make_unique< RefProviderByPtr< VectorType > >(sz))
  {}

  VectorProvider(const VectorProvider< VectorType >& other)
    : ref_provider_(Stuff::Common::make_unique< RefProviderByPtr< VectorType > >(other.vector_ref().size()))
  {
    vector_ref() = other.vector_ref();
  }

  VectorProvider(VectorProvider< VectorType >&& source)
    : ref_provider_(std::move(source.ref_provider_))
  {}

  VectorProvider< VectorType >& operator=(const VectorProvider< VectorType >& other)
  {
    if (this != &other) {
      vector_ref() = other.vector_ref();
    }
    return *this;
  }

  VectorType& vector_ref()
  {
    return ref_provider_->ref();
  }

  const VectorType& vector_ref() const
  {
    return ref_provider_->ref();
  }

private:
  std::unique_ptr< RefProvider< VectorType > > ref_provider_;
};


} // namespace internal


template< class SpaceImp, class VectorImp >
class DiscreteFunction
  : internal::VectorProvider< VectorImp >
  , public ConstDiscreteFunction< SpaceImp, VectorImp >
{
  typedef internal::VectorProvider< VectorImp > VectorProviderBaseType;
  typedef ConstDiscreteFunction< SpaceImp, VectorImp >  BaseType;
  typedef DiscreteFunction< SpaceImp, VectorImp >       ThisType;
public:
  typedef typename BaseType::SpaceType          SpaceType;
  typedef typename BaseType::VectorType         VectorType;
  typedef typename BaseType::EntityType         EntityType;
  typedef typename BaseType::LocalfunctionType  LocalfunctionType;

  typedef LocalDiscreteFunction< SpaceType, VectorType > LocalDiscreteFunctionType;

  DiscreteFunction(const SpaceType& sp,
                   VectorType& vec,
                   const std::string nm = "dune.gdt.discretefunction")
    : VectorProviderBaseType(vec)
    , BaseType(sp, VectorProviderBaseType::vector_ref(), nm)
  {}

  DiscreteFunction(const SpaceType& sp,
                   const std::string nm = "dune.gdt.discretefunction")
    : VectorProviderBaseType(sp.mapper().size())
    , BaseType(sp, VectorProviderBaseType::vector_ref(), nm)
  {}

  DiscreteFunction(const ThisType& other)
    : VectorProviderBaseType(other)
    , BaseType(other.space(), VectorProviderBaseType::vector_ref(), other.name())
  {}

  DiscreteFunction(ThisType&& source)
    : VectorProviderBaseType(source)
    , BaseType(source.space(), VectorProviderBaseType::vector_ref(), source.name())
  {}

  ~DiscreteFunction() {}

  ThisType& operator=(const ThisType& other) = delete;

  virtual ThisType* copy() const DS_OVERRIDE
  {
    return new ThisType(*this);
  }

  using BaseType::vector;

  VectorType& vector()
  {
    return this->vector_ref();
  }

  using BaseType::local_discrete_function;

  LocalDiscreteFunctionType local_discrete_function(const EntityType& entity)
  {
    assert(space_.grid_view()->indexSet().contains(entity));
    return LocalDiscreteFunctionType(space_, this->vector_ref(), entity);
  }

private:
  using BaseType::space_;
}; // class DiscreteFunction


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETEFUNCTION_DEFAULT_HH
