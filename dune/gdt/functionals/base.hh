// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tim Keil        (2017)

#ifndef DUNE_GDT_FUNCTIONALS_BASE_HH
#define DUNE_GDT_FUNCTIONALS_BASE_HH

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/grid/walker/apply-on.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/la/container/vector-interface.hh>

#include <dune/gdt/assembler/wrapper.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/functionals/interfaces.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/type_traits.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forward, required for the traits
template <class V, class S, class GL = typename S::GridLayerType, class F = typename V::RealType>
class VectorFunctionalBase;


namespace internal {


template <class VectorImp, class SpaceImp, class GridLayerImp, class FieldImp>
class VectorFunctionalBaseTraits
{
  static_assert(XT::LA::is_vector<VectorImp>::value, "VectorType has to be derived from XT::LA::vectorInterface!");
  static_assert(is_space<SpaceImp>::value, "SpaceType has to be derived from SpaceInterface!");
  static_assert(std::is_same<XT::Grid::extract_entity_t<typename SpaceImp::GridLayerType>,
                             XT::Grid::extract_entity_t<GridLayerImp>>::value,
                "SpaceType and GridLayerType have to match!");

public:
  typedef VectorFunctionalBase<VectorImp, SpaceImp, GridLayerImp, FieldImp> derived_type;
  typedef FieldImp FieldType;
};


} // namespace internal


/**
 * \note Does a const_cast in apply(), not sure yet if this is fine.
 */
template <class VectorImp, class SpaceImp, class GridLayerImp, class FieldImp>
class VectorFunctionalBase
    : public FunctionalInterface<internal::VectorFunctionalBaseTraits<VectorImp, SpaceImp, GridLayerImp, FieldImp>>,
      public SystemAssembler<SpaceImp, GridLayerImp>
{
  typedef FunctionalInterface<internal::VectorFunctionalBaseTraits<VectorImp, SpaceImp, GridLayerImp, FieldImp>>
      BaseFunctionalType;
  typedef SystemAssembler<SpaceImp, GridLayerImp> BaseAssemblerType;
  typedef VectorFunctionalBase<VectorImp, SpaceImp, GridLayerImp, FieldImp> ThisType;

public:
  typedef internal::VectorFunctionalBaseTraits<VectorImp, SpaceImp, GridLayerImp, FieldImp> Traits;
  typedef typename BaseAssemblerType::AnsatzSpaceType SpaceType;
  typedef typename SpaceType::BaseFunctionSetType TestBaseType;
  using typename BaseAssemblerType::GridLayerType;
  typedef VectorImp VectorType;
  using typename BaseAssemblerType::IntersectionType;
  using typename BaseFunctionalType::FieldType;
  using typename BaseFunctionalType::derived_type;

  template <class... Args>
  explicit VectorFunctionalBase(VectorType& vec, Args&&... args)
    : BaseAssemblerType(std::forward<Args>(args)...)
    , vector_(vec)
  {
    if (vector_.access().size() != this->test_space().mapper().size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "vector.size(): " << vector_.access().size() << "\n"
                                   << "space().mapper().size(): "
                                   << this->space().mapper().size());
  } // VectorFunctionalBase(...)

  /// \todo Guard against copy and move ctor (Args = ThisType)!
  template <class... Args>
  explicit VectorFunctionalBase(Args&&... args)
    : BaseAssemblerType(std::forward<Args>(args)...)
    , vector_(new VectorType(this->test_space().mapper().size(), 0.0))
  {
  }

  /// \sa SystemAssembler
  VectorFunctionalBase(const ThisType& other) = delete;
  VectorFunctionalBase(ThisType&& source) = delete;
  VectorFunctionalBase(ThisType& other) = delete; // <- b.c. of the too perfect forwarding ctor

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  const VectorType& vector() const
  {
    return vector_.access();
  }

  VectorType& vector()
  {
    return vector_.access();
  }

  const SpaceType& space() const
  {
    return this->ansatz_space();
  }

  using BaseAssemblerType::append;

  ThisType& append(
      const LocalVolumeFunctionalInterface<TestBaseType, FieldType>& local_volume_functional,
      const XT::Grid::ApplyOn::WhichEntity<GridLayerType>* where = new XT::Grid::ApplyOn::AllEntities<GridLayerType>())
  {
    this->append(local_volume_functional, vector_.access(), where);
    return *this;
  }

  ThisType& append(const LocalFaceFunctionalInterface<TestBaseType, IntersectionType, FieldType>& local_face_functional,
                   const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where =
                       new XT::Grid::ApplyOn::AllIntersections<GridLayerType>())
  {
    this->append(local_face_functional, vector_.access(), where);
    return *this;
  }

  template <class S>
  FieldType apply(const XT::LA::VectorInterface<S>& source) const
  {
    const_cast<ThisType&>(*this).assemble();
    return vector().dot(source.as_imp());
  }

  template <class S>
  FieldType apply(const ConstDiscreteFunction<SpaceType, S>& source) const
  {
    return apply(source.vector());
  }

private:
  Dune::XT::Common::StorageProvider<VectorType> vector_;
}; // class VectorFunctionalBase


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_FUNCTIONALS_BASE_HH
