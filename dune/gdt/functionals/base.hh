// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#ifndef DUNE_GDT_FUNCTIONALS_BASE_HH
#define DUNE_GDT_FUNCTIONALS_BASE_HH

#include <dune/xt/common/exceptions.hh>
#include <dune/stuff/grid/walker/apply-on.hh>
#include <dune/stuff/la/container/vector-interface.hh>

#include <dune/gdt/assembler/wrapper.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/functionals/interfaces.hh>
#include <dune/gdt/spaces/interface.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forward, required for the traits
template <class V, class S, class GV = typename S::GridViewType, class F = typename V::RealType>
class VectorFunctionalBase;


namespace internal {


template <class VectorImp, class SpaceImp, class GridViewImp, class FieldImp>
class VectorFunctionalBaseTraits
{
  static_assert(Stuff::LA::is_vector<VectorImp>::value,
                "VectorType has to be derived from Stuff::LA::vectorInterface!");
  static_assert(is_space<SpaceImp>::value, "SpaceType has to be derived from SpaceInterface!");
  static_assert(std::is_same<typename SpaceImp::GridViewType::template Codim<0>::Entity,
                             typename GridViewImp::template Codim<0>::Entity>::value,
                "SpaceType and GridViewType have to match!");

public:
  typedef VectorFunctionalBase<VectorImp, SpaceImp, GridViewImp, FieldImp> derived_type;
  typedef FieldImp FieldType;
};


} // namespace internal


/**
 * \note Does a const_cast in apply(), not sure yet if this is fine.
 */
template <class VectorImp, class SpaceImp, class GridViewImp, class FieldImp>
class VectorFunctionalBase
    : public FunctionalInterface<internal::VectorFunctionalBaseTraits<VectorImp, SpaceImp, GridViewImp, FieldImp>>,
      public SystemAssembler<SpaceImp, GridViewImp>
{
  typedef FunctionalInterface<internal::VectorFunctionalBaseTraits<VectorImp, SpaceImp, GridViewImp, FieldImp>>
      BaseFunctionalType;
  typedef SystemAssembler<SpaceImp, GridViewImp> BaseAssemblerType;
  typedef VectorFunctionalBase<VectorImp, SpaceImp, GridViewImp, FieldImp> ThisType;

public:
  typedef internal::VectorFunctionalBaseTraits<VectorImp, SpaceImp, GridViewImp, FieldImp> Traits;
  typedef typename BaseAssemblerType::AnsatzSpaceType SpaceType;
  using typename BaseAssemblerType::GridViewType;
  typedef VectorImp VectorType;
  using typename BaseFunctionalType::FieldType;
  using typename BaseFunctionalType::derived_type;

public:
  template <class... Args>
  explicit VectorFunctionalBase(VectorType& vec, Args&&... args)
    : BaseAssemblerType(std::forward<Args>(args)...)
    , vector_(vec)
  {
    if (vector_.access().size() != this->test_space().mapper().size())
      DUNE_THROW(Stuff::Exceptions::shapes_do_not_match,
                 "vector.size(): " << vector_.access().size() << "\n"
                                   << "space().mapper().size(): "
                                   << this->space().mapper().size());
  } // VectorFunctionalBase(...)

  template <class... Args>
  explicit VectorFunctionalBase(Args&&... args)
    : BaseAssemblerType(std::forward<Args>(args)...)
    , vector_(new VectorType(this->test_space().mapper().size(), 0.0))
  {
  }

  VectorFunctionalBase(ThisType&& source) = default;

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

  using BaseAssemblerType::add;

  template <class F>
  void add(const LocalVolumeFunctionalInterface<F>& local_volume_functional,
           const DSG::ApplyOn::WhichEntity<GridViewType>* where = new DSG::ApplyOn::AllEntities<GridViewType>())
  {
    typedef internal::LocalVolumeFunctionalWrapper<ThisType,
                                                   typename LocalVolumeFunctionalInterface<F>::derived_type,
                                                   VectorType>
        WrapperType;
    this->codim0_functors_.emplace_back(
        new WrapperType(this->test_space_, where, local_volume_functional.as_imp(), vector_.access()));
  }

  template <class F>
  void
  add(const LocalFaceFunctionalInterface<F>& local_face_functional,
      const DSG::ApplyOn::WhichIntersection<GridViewType>* where = new DSG::ApplyOn::AllIntersections<GridViewType>())
  {
    typedef internal::LocalFaceFunctionalWrapper<ThisType,
                                                 typename LocalFaceFunctionalInterface<F>::derived_type,
                                                 VectorType>
        WrapperType;
    this->codim1_functors_.emplace_back(
        new WrapperType(this->test_space_, where, local_face_functional.as_imp(), vector_.access()));
  }

  template <class S>
  FieldType apply(const Stuff::LA::VectorInterface<S>& source) const
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
