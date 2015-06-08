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
#include <dune/common/fvector.hh>

#include <dune/grid/io/file/vtk.hh>

#include <dune/stuff/la/container/interfaces.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/common/memory.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/productinterface.hh>

#include "local.hh"

namespace Dune {
namespace GDT {

//forward
template< class SpaceImp, class VectorImp >
class ConstDiscreteFunction;


namespace internal {


template< size_t ii, class SpaceImp, class VectorImp, class Traits, bool is_factor_space = false >
struct visualize_helper
{
  static void visualize_factor(const std::string filename,
                               const bool subsampling,
                               const VTK::OutputType vtk_output_type,
                               const DS::PerThreadValue< SpaceImp >& space,
                               const VectorImp& vector)
  {
    assert( ii == 0 );
    ConstDiscreteFunction< SpaceImp, VectorImp > discrete_function(space, vector);
    discrete_function.visualize(filename, subsampling, vtk_output_type);
    std::cout << std::is_base_of< typename Dune::GDT::ProductSpaceInterface< Traits >, SpaceImp >::value;
  }
};

template< size_t ii, class SpaceImp, class VectorImp, class Traits >
struct visualize_helper< ii, SpaceImp, VectorImp, Traits, true  > {
  static void visualize_factor(const std::string filename,
                               const bool subsampling,
                               const VTK::OutputType vtk_output_type,
                               const DS::PerThreadValue< SpaceImp >& space,
                               const VectorImp& vector)
  {
    assert(ii < SpaceImp::num_factors);
    const auto& factor_space = space->template factor< ii >();
    VectorImp factor_vector(factor_space.mapper().size());
    const auto it_end = space->grid_view().template end< 0 >();
    for (auto it = space->grid_view().template begin< 0 >(); it != it_end; ++it) {
      const auto& entity = *it;
      for (size_t jj = 0; jj < factor_space.mapper().numDofs(entity); ++jj)
        factor_vector[factor_space.mapper().mapToGlobal(entity, jj)] = vector[space->factor_mapper().mapToGlobal(ii, entity, jj)];
    }
    ConstDiscreteFunction< typename SpaceImp::FactorSpaceType, VectorImp > factor_discrete_function(factor_space, factor_vector);
    factor_discrete_function.visualize(filename, subsampling, vtk_output_type);
  }
};


} // namespace internal


template< class SpaceImp, class VectorImp >
class ConstDiscreteFunction
  : public Stuff::LocalizableFunctionInterface< typename SpaceImp::EntityType,
                                                typename SpaceImp::DomainFieldType, SpaceImp::dimDomain,
                                                typename SpaceImp::RangeFieldType, SpaceImp::dimRange, SpaceImp::dimRangeCols >
{
  static_assert(is_space< SpaceImp >::value, "SpaceImp has to be derived from SpaceInterface!");
  static_assert(Stuff::LA::is_vector< VectorImp >::value,
                "VectorImp has to be derived from Stuff::LA::VectorInterface!");
  static_assert(std::is_same< typename SpaceImp::RangeFieldType, typename VectorImp::ScalarType >::value,
                "Types do not match!");
  typedef Stuff::LocalizableFunctionInterface
      < typename SpaceImp::EntityType, typename SpaceImp::DomainFieldType, SpaceImp::dimDomain,
        typename SpaceImp::RangeFieldType, SpaceImp::dimRange, SpaceImp::dimRangeCols >
                                                       BaseType;
  typedef ConstDiscreteFunction< SpaceImp, VectorImp > ThisType;
public:
  typedef SpaceImp                             SpaceType;
  typedef typename SpaceImp::Traits            SpaceTraits;
  typedef VectorImp                            VectorType;
  typedef typename BaseType::EntityType        EntityType;
  typedef typename BaseType::LocalfunctionType LocalfunctionType;

  typedef ConstLocalDiscreteFunction< SpaceType, VectorType > ConstLocalDiscreteFunctionType;

  ConstDiscreteFunction(const SpaceType& sp, const VectorType& vec, const std::string nm = "gdt.constdiscretefunction")
    : space_(sp)
    , vector_(vec)
    , name_(nm)
  {
    assert(vector_.size() == space_->mapper().size() && "Given vector has wrong size!");
  }

  ConstDiscreteFunction(const ThisType& other)
    : space_(other.space())
    , vector_(other.vector_)
    , name_(other.name_)
  {}

  ConstDiscreteFunction(ThisType&& source)
    : space_(source.space())
    , vector_(source.vector_)
    , name_(source.name_)
  {}

  virtual ~ConstDiscreteFunction() {}

  ThisType& operator=(const ThisType& other) = delete;

  virtual std::string name() const override
  {
    return name_;
  }

  const SpaceType& space() const
  {
    return *space_;
  }

  const VectorType& vector() const
  {
    return vector_;
  }

  std::unique_ptr< ConstLocalDiscreteFunctionType > local_discrete_function(const EntityType& entity) const
  {
    assert(space_->grid_view().indexSet().contains(entity));
    return DSC::make_unique< ConstLocalDiscreteFunctionType >(*space_, vector_, entity);
  }

  virtual std::unique_ptr< LocalfunctionType > local_function(const EntityType& entity) const override
  {
    return local_discrete_function(entity);
  }

  void visualize(const std::string filename,
                 const bool subsampling = (SpaceType::polOrder > 1),
                 const VTK::OutputType vtk_output_type = VTK::appendedraw) const
  {
    BaseType::template visualize< typename SpaceType::GridViewType >(space().grid_view(),
                                                                     filename,
                                                                     subsampling,
                                                                     vtk_output_type);
  } // ... visualize(...)

  template< size_t ii >
  void visualize_factor(const std::string filename,
                        const bool subsampling = (SpaceType::polOrder > 1),
                        const VTK::OutputType vtk_output_type = VTK::appendedraw) const
  {
    internal::visualize_helper
        < ii,
          SpaceType,
          VectorType,
          SpaceTraits,
          std::is_base_of< typename Dune::GDT::ProductSpaceInterface< SpaceTraits >, SpaceType >::value >
                              ::visualize_factor(filename, subsampling, vtk_output_type, space_, vector_);
  } // ... visualize_factor< ii >(...)

  bool dofs_valid() const
  {
    return vector().valid();
  }

protected:
  const DS::PerThreadValue< SpaceType > space_;
private:
  const VectorType& vector_;
  const std::string name_;
}; // class ConstDiscreteFunction


template< class SpaceImp, class VectorImp >
class DiscreteFunction
  : Stuff::Common::StorageProvider< VectorImp >
  , public ConstDiscreteFunction< SpaceImp, VectorImp >
{
  typedef Stuff::Common::StorageProvider< VectorImp >  VectorProviderBaseType;
  typedef ConstDiscreteFunction< SpaceImp, VectorImp > BaseType;
  typedef DiscreteFunction< SpaceImp, VectorImp >      ThisType;
public:
  typedef typename BaseType::SpaceType         SpaceType;
  typedef typename BaseType::VectorType        VectorType;
  typedef typename BaseType::EntityType        EntityType;
  typedef typename BaseType::LocalfunctionType LocalfunctionType;

  typedef LocalDiscreteFunction< SpaceType, VectorType > LocalDiscreteFunctionType;

  DiscreteFunction(const SpaceType& sp, VectorType& vec, const std::string nm = "gdt.discretefunction")
    : VectorProviderBaseType(vec)
    , BaseType(sp, VectorProviderBaseType::storage_access(), nm)
  {}

  DiscreteFunction(const SpaceType& sp, VectorType&& vec, const std::string nm = "gdt.discretefunction")
    : VectorProviderBaseType(vec)
    , BaseType(sp, VectorProviderBaseType::storage_access(), nm)
  {}

  DiscreteFunction(const SpaceType& sp, const std::string nm = "gdt.discretefunction")
    : VectorProviderBaseType(new VectorType(sp.mapper().size()))
    , BaseType(sp, VectorProviderBaseType::storage_access(), nm)
  {}

  // manual copy ctor needed bc. of the storage provider
  DiscreteFunction(const ThisType& other)
    : VectorProviderBaseType(new VectorType(other.vector()))
    , BaseType(other.space(), VectorProviderBaseType::storage_access(), other.name())
  {}

  // manual move ctor needed bc. of the storage provider
  DiscreteFunction(ThisType&& source)
    : VectorProviderBaseType(new VectorType(source.vector()))
    , BaseType(source.space(), VectorProviderBaseType::storage_access(), source.name())
  {}

  virtual ~DiscreteFunction() {}

  ThisType& operator=(const ThisType& other) = delete;

  using BaseType::vector;

  VectorType& vector()
  {
    return this->storage_access();
  }

  using BaseType::local_discrete_function;

  std::unique_ptr< LocalDiscreteFunctionType > local_discrete_function(const EntityType& entity)
  {
    assert(space_->grid_view().indexSet().contains(entity));
    return DSC::make_unique< LocalDiscreteFunctionType >(*space_, this->storage_access(), entity);
  }

private:
  using BaseType::space_;
}; // class DiscreteFunction


template< class SpaceType, class VectorType >
ConstDiscreteFunction< SpaceType,
                       VectorType > make_const_discrete_function(const SpaceType& space,
                                                                 const VectorType& vector,
                                                                 const std::string nm = "gdt.constdiscretefunction")
{
  return ConstDiscreteFunction< SpaceType, VectorType >(space, vector, nm);
}


template< class SpaceType, class VectorType >
typename std::enable_if< is_space< SpaceType >::value && Stuff::LA::is_vector< VectorType >::value,
                         DiscreteFunction< SpaceType, VectorType > >::type
make_discrete_function(const SpaceType& space, VectorType& vector, const std::string nm = "gdt.discretefunction")
{
  return DiscreteFunction< SpaceType, VectorType >(space, vector, nm);
}


/**
 * This can be used like \code
auto discrete_function = make_discrete_function< VectorType >(space);
\endcode
 */
template< class VectorType, class SpaceType >
DiscreteFunction< SpaceType, VectorType > make_discrete_function(const SpaceType& space,
                                                                 const std::string nm = "gdt.discretefunction")
{
  return DiscreteFunction< SpaceType, VectorType >(space, nm);
}


namespace internal {


template< class D >
struct is_const_discrete_function_helper
{
  DSC_has_typedef_initialize_once(SpaceType)
  DSC_has_typedef_initialize_once(VectorType)

  static const bool is_candidate = DSC_has_typedef(SpaceType)< D >::value && DSC_has_typedef(SpaceType)< D >::value;
};


} // namespace internal


template< class D, bool candidate = internal::is_const_discrete_function_helper< D >::is_candidate >
struct is_const_discrete_function
  : public std::is_base_of< ConstDiscreteFunction< typename D::SpaceType, typename D::VectorType >, D >
{};


template< class D >
struct is_const_discrete_function< D, false >
  : public std::false_type
{};


template< class D, bool candidate = internal::is_const_discrete_function_helper< D >::is_candidate >
struct is_discrete_function
  : public std::is_base_of< DiscreteFunction< typename D::SpaceType, typename D::VectorType >, D >
{};


template< class D >
struct is_discrete_function< D, false >
  : public std::false_type
{};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETEFUNCTION_DEFAULT_HH
