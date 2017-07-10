// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2011 - 2017)
//   Kirsten Weber   (2013)
//   Rene Milk       (2014 - 2017)
//   Tobias Leibner  (2014, 2016)

#ifndef DUNE_GDT_DISCRETEFUNCTION_DEFAULT_HH
#define DUNE_GDT_DISCRETEFUNCTION_DEFAULT_HH

#include <memory>
#include <vector>
#include <type_traits>

#include <dune/common/fvector.hh>

#include <dune/grid/io/file/vtk.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/tuple.hh>

#include <dune/xt/la/container.hh>

#include <dune/xt/functions/interfaces.hh>

#include <dune/gdt/local/discretefunction.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/fv/default.hh>

namespace Dune {
namespace GDT {

// forward
template <class SpaceImp, class VectorImp>
class ConstDiscreteFunction;


namespace internal {


template <size_t ii, bool is_product_space = false>
struct visualize_helper
{
  template <class DiscreteFunctionType>
  static void visualize_factor(const std::string filename_prefix,
                               const std::string filename_suffix,
                               const bool subsampling,
                               const VTK::OutputType vtk_output_type,
                               const DiscreteFunctionType& discrete_function)
  {
    static_assert(ii == 0, "Space is not a product space, so there is no factor other than 0.");
    discrete_function.visualize(filename_prefix + "_factor_" + Dune::XT::Common::to_string(ii) + "_" + filename_suffix,
                                subsampling,
                                vtk_output_type);
  }
};

template <size_t ii>
struct visualize_helper<ii, true>
{
  template <class DiscreteFunctionType>
  static void visualize_factor(const std::string filename_prefix,
                               const std::string filename_suffix,
                               const bool subsampling,
                               const VTK::OutputType vtk_output_type,
                               const DiscreteFunctionType& discrete_function)
  {
    static_assert(ii < DiscreteFunctionType::SpaceType::num_factors, "This factor does not exist.");
    const auto& space = discrete_function.space();
    const auto& factor_space = space.template factor<ii>();
    typename DiscreteFunctionType::VectorType factor_vector(factor_space.mapper().size());
    const auto it_end = space.grid_layer().template end<0>();
    for (auto it = space.grid_layer().template begin<0>(); it != it_end; ++it) {
      const auto& entity = *it;
      for (size_t jj = 0; jj < factor_space.mapper().numDofs(entity); ++jj)
        factor_vector.set_entry(factor_space.mapper().mapToGlobal(entity, jj),
                                discrete_function.vector().get_entry(space.mapper().mapToGlobal(ii, entity, jj)));
    }
    ConstDiscreteFunction<
        typename XT::Common::tuple_element<ii, typename DiscreteFunctionType::SpaceType::SpaceTupleType>::type,
        typename DiscreteFunctionType::VectorType>
        factor_discrete_function(factor_space, factor_vector);
    factor_discrete_function.visualize(filename_prefix + "_factor_" + Dune::XT::Common::to_string(ii) + "_"
                                           + filename_suffix,
                                       subsampling,
                                       vtk_output_type);
  }
};

// discrete_exact_solution_ time for loop to visualize all factors of a product space, see
// http://stackoverflow.com/a/11081785
template <size_t index, size_t N>
struct static_for_loop
{
  template <class DiscreteFunctionType>
  static void visualize(const std::string filename_prefix,
                        const std::string filename_suffix,
                        const bool subsampling,
                        const VTK::OutputType vtk_output_type,
                        const DiscreteFunctionType& discrete_function)
  {
    static_for_loop<index, N / 2>::visualize(
        filename_prefix, filename_suffix, subsampling, vtk_output_type, discrete_function);
    static_for_loop<index + N / 2, N - N / 2>::visualize(
        filename_prefix, filename_suffix, subsampling, vtk_output_type, discrete_function);
  }
};

// specialization of static for loop to end the loop
template <size_t index>
struct static_for_loop<index, 1>
{
  template <class DiscreteFunctionType>
  static void visualize(const std::string filename_prefix,
                        const std::string filename_suffix,
                        const bool subsampling,
                        const VTK::OutputType vtk_output_type,
                        const DiscreteFunctionType& discrete_function)
  {
    visualize_helper<index, true>::visualize_factor(
        filename_prefix, filename_suffix, subsampling, vtk_output_type, discrete_function);
  }
};

template <size_t index>
struct static_for_loop<index, 0>
{
  template <class DiscreteFunctionType>
  static void visualize(const std::string /*filename_prefix*/,
                        const std::string /*filename_suffix*/,
                        const bool /*subsampling*/,
                        const VTK::OutputType /*vtk_output_type*/,
                        const DiscreteFunctionType& /*discrete_function*/)
  {
  }
};


} // namespace internal


template <class SpaceImp, class VectorImp>
class ConstDiscreteFunction : public XT::Functions::LocalizableFunctionInterface<typename SpaceImp::EntityType,
                                                                                 typename SpaceImp::DomainFieldType,
                                                                                 SpaceImp::dimDomain,
                                                                                 typename SpaceImp::RangeFieldType,
                                                                                 SpaceImp::dimRange,
                                                                                 SpaceImp::dimRangeCols>
{
  static_assert(is_space<SpaceImp>::value, "SpaceImp has to be derived from SpaceInterface!");
  static_assert(XT::LA::is_vector<VectorImp>::value, "VectorImp has to be derived from XT::LA::VectorInterface!");
  static_assert(std::is_same<typename SpaceImp::RangeFieldType, typename VectorImp::ScalarType>::value,
                "Types do not match!");
  typedef XT::Functions::LocalizableFunctionInterface<typename SpaceImp::EntityType,
                                                      typename SpaceImp::DomainFieldType,
                                                      SpaceImp::dimDomain,
                                                      typename SpaceImp::RangeFieldType,
                                                      SpaceImp::dimRange,
                                                      SpaceImp::dimRangeCols>
      BaseType;
  typedef XT::Common::ConstStorageProvider<VectorImp> VectorStorageProvider;
  typedef ConstDiscreteFunction<SpaceImp, VectorImp> ThisType;

public:
  typedef SpaceImp SpaceType;
  typedef typename SpaceType::Traits SpaceTraits;
  typedef VectorImp VectorType;
  using typename BaseType::EntityType;
  using typename BaseType::LocalfunctionType;

  typedef ConstLocalDiscreteFunction<SpaceType, VectorType> ConstLocalDiscreteFunctionType;

  ConstDiscreteFunction(const SpaceType& sp, const VectorType& vec, const std::string nm = "gdt.constdiscretefunction")
    : space_(new XT::Common::PerThreadValue<SpaceType>(sp))
    , vector_(new VectorStorageProvider(vec))
    , name_(nm)
  {
    if (vector().size() != (*space_)->mapper().size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "space.mapper().size(): " << (*space_)->mapper().size() << "\n   "
                                           << "vector.size(): "
                                           << vector_->access().size());
  }

  ConstDiscreteFunction(const ThisType& other)
    : space_(new XT::Common::PerThreadValue<SpaceType>(other.space()))
    , vector_(new VectorStorageProvider(other.vector()))
    , name_(other.name_)
  {
  }

  ConstDiscreteFunction(ThisType&& source)
    : space_(std::move(source.space_))
    , vector_(std::move(source.vector_))
    , name_(source.name_)
  {
  }

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = default;

  std::string name() const override final
  {
    return name_;
  }

  const SpaceType& space() const
  {
    return *(*space_);
  }

  const VectorType& vector() const
  {
    return vector_->access();
  }

  std::unique_ptr<ConstLocalDiscreteFunctionType> local_discrete_function(const EntityType& entity) const
  {
    assert((*space_)->grid_layer().indexSet().contains(entity));
    return Dune::XT::Common::make_unique<ConstLocalDiscreteFunctionType>(*(*space_), vector_->access(), entity);
  }

  std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const override final
  {
    return local_discrete_function(entity);
  }

  using BaseType::visualize;

  /**
   * \brief Visualizes the function using Dune::XT::Functions::LocalizableFunctionInterface::visualize on the grid layer
   *        associated with the space.
   * \sa    Dune::XT::Functions::LocalizableFunctionInterface::visualize
   * \note  Subsampling is enabled by default for functions of order greater than one.
   */
  void visualize(const std::string filename,
                 const bool subsampling = (SpaceType::polOrder > 1),
                 const VTK::OutputType vtk_output_type = VTK::appendedraw) const
  {
    visualize(filename, "", subsampling, vtk_output_type);
  }

  void visualize(const std::string filename_prefix,
                 const std::string filename_suffix,
                 const bool subsampling = (SpaceType::polOrder > 1),
                 const VTK::OutputType vtk_output_type = VTK::appendedraw) const
  {
    redirect_visualize<>()(*this, filename_prefix, filename_suffix, subsampling, vtk_output_type);
  }

  template <size_t ii>
  void visualize_factor(const std::string filename_prefix,
                        const std::string filename_suffix = "",
                        const bool subsampling = (SpaceType::polOrder > 1),
                        const VTK::OutputType vtk_output_type = VTK::appendedraw) const
  {
    internal::visualize_helper<ii, is_product_space<SpaceType>::value>::visualize_factor(
        filename_prefix, filename_suffix, subsampling, vtk_output_type, *this);
  }

  bool dofs_valid() const
  {
    return vector().valid();
  }

protected:
  template <bool product_space = is_product_space<SpaceType>::value, class Anything = void>
  struct redirect_visualize
  {
    void operator()(const ThisType& self,
                    const std::string filename_prefix,
                    const std::string filename_suffix,
                    const bool subsampling,
                    const VTK::OutputType vtk_output_type)
    {
      static_cast<const BaseType&>(self).template visualize<typename SpaceType::GridLayerType>(
          self.space().grid_layer(), filename_prefix + filename_suffix, subsampling, vtk_output_type);
    }
  };

  template <class Anything>
  struct redirect_visualize<true, Anything>
  {
    void operator()(const ThisType& self,
                    const std::string filename_prefix,
                    const std::string filename_suffix,
                    const bool subsampling,
                    const VTK::OutputType vtk_output_type)
    {
      internal::static_for_loop<0,
                                ProductSpaceInterface<SpaceTraits,
                                                      SpaceType::dimDomain,
                                                      SpaceType::dimRange,
                                                      SpaceType::dimRangeCols>::num_factors>::visualize(filename_prefix,
                                                                                                        filename_suffix,
                                                                                                        subsampling,
                                                                                                        vtk_output_type,
                                                                                                        self);
    }
  };

  std::unique_ptr<XT::Common::PerThreadValue<SpaceType>> space_;

private:
  std::unique_ptr<VectorStorageProvider> vector_;
  std::string name_;
}; // class ConstDiscreteFunction


template <class SpaceImp, class VectorImp = typename XT::LA::Container<typename SpaceImp::RangeFieldType>::VectorType>
class DiscreteFunction : XT::Common::StorageProvider<VectorImp>, public ConstDiscreteFunction<SpaceImp, VectorImp>
{
  typedef XT::Common::StorageProvider<VectorImp> VectorProviderBaseType;
  typedef ConstDiscreteFunction<SpaceImp, VectorImp> BaseType;
  typedef DiscreteFunction<SpaceImp, VectorImp> ThisType;

public:
  typedef typename BaseType::SpaceType SpaceType;
  typedef typename BaseType::VectorType VectorType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::LocalfunctionType LocalfunctionType;

  typedef LocalDiscreteFunction<SpaceType, VectorType> LocalDiscreteFunctionType;

  DiscreteFunction(const SpaceType& sp, VectorType& vec, const std::string nm = "gdt.discretefunction")
    : VectorProviderBaseType(vec)
    , BaseType(sp, VectorProviderBaseType::access(), nm)
  {
  }

  DiscreteFunction(const SpaceType& sp, VectorType&& vec, const std::string nm = "gdt.discretefunction")
    : VectorProviderBaseType(std::move(vec))
    , BaseType(sp, VectorProviderBaseType::access(), nm)
  {
  }

  DiscreteFunction(const SpaceType& sp, const std::string nm = "gdt.discretefunction")
    : VectorProviderBaseType(new VectorType(sp.mapper().size(), 0, 2 * XT::Common::threadManager().max_threads()))
    , BaseType(sp, VectorProviderBaseType::access(), nm)
  {
  }

  DiscreteFunction(const ThisType& other)
    : VectorProviderBaseType(VectorType(other.vector()))
    , BaseType(other.space(), VectorProviderBaseType::access(), other.name())
  {
  }

  DiscreteFunction(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = default;

  using BaseType::vector;

  VectorType& vector()
  {
    return this->access();
  }

  using BaseType::local_discrete_function;

  std::unique_ptr<LocalDiscreteFunctionType> local_discrete_function(const EntityType& entity)
  {
    assert((*space_)->grid_layer().indexSet().contains(entity));
    return Dune::XT::Common::make_unique<LocalDiscreteFunctionType>(*(*space_), this->access(), entity);
  }

private:
  using BaseType::space_;
}; // class DiscreteFunction


template <class SpaceType, class VectorType>
ConstDiscreteFunction<SpaceType, VectorType> make_const_discrete_function(
    const SpaceType& space, const VectorType& vector, const std::string nm = "gdt.constdiscretefunction")
{
  return ConstDiscreteFunction<SpaceType, VectorType>(space, vector, nm);
}


template <class SpaceType, class VectorType>
typename std::enable_if<is_space<SpaceType>::value && XT::LA::is_vector<VectorType>::value,
                        DiscreteFunction<SpaceType, VectorType>>::type
make_discrete_function(const SpaceType& space, VectorType& vector, const std::string nm = "gdt.discretefunction")
{
  return DiscreteFunction<SpaceType, VectorType>(space, vector, nm);
}

template <class SpaceType, class VectorType>
typename std::enable_if<is_space<SpaceType>::value && XT::LA::is_vector<VectorType>::value,
                        DiscreteFunction<SpaceType, VectorType>>::type
make_discrete_function(const SpaceType& space, VectorType&& vector, const std::string nm = "gdt.discretefunction")
{
  return DiscreteFunction<SpaceType, VectorType>(space, std::move(vector), nm);
}


/**
 * This can be used like \code
auto discrete_function = make_discrete_function<VectorType>(space);
\endcode
 */
template <class VectorType, class SpaceType>
DiscreteFunction<SpaceType, VectorType> make_discrete_function(const SpaceType& space,
                                                               const std::string nm = "gdt.discretefunction")
{
  return DiscreteFunction<SpaceType, VectorType>(space, nm);
}


namespace internal {


template <class D>
struct is_const_discrete_function_helper
{
  DXTC_has_typedef_initialize_once(SpaceType);
  DXTC_has_typedef_initialize_once(VectorType);

  static const bool is_candidate = DXTC_has_typedef(SpaceType)<D>::value && DXTC_has_typedef(VectorType)<D>::value;
};


} // namespace internal


template <class D, bool candidate = internal::is_const_discrete_function_helper<D>::is_candidate>
struct is_const_discrete_function
    : public std::is_base_of<ConstDiscreteFunction<typename D::SpaceType, typename D::VectorType>, D>
{
};


template <class D>
struct is_const_discrete_function<D, false> : public std::false_type
{
};


template <class D, bool candidate = internal::is_const_discrete_function_helper<D>::is_candidate>
struct is_discrete_function : public std::is_base_of<DiscreteFunction<typename D::SpaceType, typename D::VectorType>, D>
{
};


template <class D>
struct is_discrete_function<D, false> : public std::false_type
{
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETEFUNCTION_DEFAULT_HH
