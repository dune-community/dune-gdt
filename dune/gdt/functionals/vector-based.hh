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

#ifndef DUNE_GDT_FUNCTIONALS_VECTOR_BASED_HH
#define DUNE_GDT_FUNCTIONALS_VECTOR_BASED_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/la/container/vector-interface.hh>
#include <dune/xt/grid/functors/interfaces.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/filters.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/assembler/functional-assemblers.hh>
#include <dune/gdt/local/functionals/interfaces.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * \brief Base class for linear functionals which are given by an assembled vector.
 *
 * \note See FunctionalInterface for a description of the template arguments.
 *
 * \sa FunctionalInterface
 * \sa VectorBasedFunctional
 */
template <class V, class GV, size_t r = 1, size_t rC = 1, class F = double>
class ConstVectorBasedFunctional : public FunctionalInterface<V, GV, r, rC, F>
{
  // All other types are checked elsewhere.
  static_assert(XT::LA::is_vector<V>::value, "");

  using ThisType = ConstVectorBasedFunctional;
  using FunctionalBaseType = FunctionalInterface<V, GV, r, rC, F>;

public:
  using DofFieldType = typename V::ScalarType;

  using typename FunctionalBaseType::FieldType;
  using typename FunctionalBaseType::SourceSpaceType;
  using typename FunctionalBaseType::SourceType;
  using typename FunctionalBaseType::SourceVectorType;

  ConstVectorBasedFunctional(const SourceSpaceType& source_spc, const SourceVectorType& vec)
    : source_space_(source_spc)
    , vector_(vec)
  {
    if (vector_.size() != source_space_.mapper().size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "vector_.size() = " << vector_.size() << "\n"
                                     << "source_space_.mapper().size() = " << source_space_.mapper().size());
  }

  ConstVectorBasedFunctional(const ThisType&) = default;
  ConstVectorBasedFunctional(ThisType&&) = default;

  bool linear() const override
  {
    return true;
  }

  const SourceSpaceType& source_space() const override
  {
    return source_space_;
  }

  const SourceVectorType& vector() const
  {
    return vector_;
  }

  FieldType apply(const SourceType& source) const override
  {
    try {
      return vector_.dot(source.dofs().vector());
    } catch (const XT::Common::Exceptions::shapes_do_not_match& ee) {
      DUNE_THROW(Exceptions::functional_error,
                 "when applying vector to source dofs!\n\nThis was the original error: " << ee.what());
    }
  } // ... apply(...)

private:
  const SourceSpaceType& source_space_;
  const SourceVectorType& vector_;
}; // class ConstVectorBasedFunctional


template <class GV, size_t r, size_t rC, class F, class V>
ConstVectorBasedFunctional<typename XT::LA::VectorInterface<V>::derived_type, GV, r, rC, F>
make_vector_functional(const SpaceInterface<GV, r, rC, F>& space, const XT::LA::VectorInterface<V>& vector)
{
  return ConstVectorBasedFunctional<typename XT::LA::VectorInterface<V>::derived_type, GV, r, rC, F>(space,
                                                                                                     vector.as_imp());
}


/**
 * \brief Base class for linear functionals which are assembled into a vector.
 *
 * Similar to the GlobalAssembler, we derive from the XT::Grid::Walker and povide custom append() methods to allow to
 * add local element and intersection functionals. In contrast to the GlobalAssembler we already hold the target vector
 * (or create one of appropriate size), into which we want to assemble. The functional is assembled by walking over the
 * given assembly_gid_view (which defaults to the one fom the given space). This allows to assemble a functional only on
 * a smaller grid view than the one given from the space (similar functionality could be acchieved by appending this
 * functional to another walker and by providing an appropriate filter).
 *
 * \note One could achieve similar functionality by deriving from GlobalAssembler directly, which would slightly
 *       simplify the implementation of the append methods. However, we do not want to expose the other append methods
 *       of GlobalAssembler here (it should not be possible to append a local opeator or two-form to a functional), but
 *       want to expose the ones from the XT::Grid::Walker (it should be possible to append other element or
 *       intersection functors).
 *
 * \note For convenience, we use the same type of vector to assemble this functional into, that we use to model the
 *       source (see the documentation in FunctionalInterface), although other choices are in general conceivable.
 *
 * \note See ConstVectorBasedFunctional and FunctionalInterface for a description of the template arguments.
 *
 * \sa ConstVectorBasedFunctional
 * \sa FunctionalInterface
 * \sa XT::Grid::Walker
 * \sa GlobalAssembler
 */
template <class V, class GV, size_t r = 1, size_t rC = 1, class F = double, class AssemblyGridView = GV>
class VectorBasedFunctional
  : XT::Common::StorageProvider<V>
  , public ConstVectorBasedFunctional<V, GV, r, rC, F>
  , public XT::Grid::Walker<AssemblyGridView>
{
  // All other types are checked elsewhere.
  static_assert(std::is_same<XT::Grid::extract_entity_t<GV>, XT::Grid::extract_entity_t<AssemblyGridView>>::value,
                "We cannot handle different element types!");

  using ThisType = VectorBasedFunctional;
  using VectorStorage = XT::Common::StorageProvider<V>;
  using FunctionalBaseType = ConstVectorBasedFunctional<V, GV, r, rC, F>;
  using WalkerBaseType = XT::Grid::Walker<AssemblyGridView>;

public:
  using AssemblyGridViewType = AssemblyGridView;
  using DofFieldType = typename V::ScalarType;

  using typename FunctionalBaseType::FieldType;
  using typename FunctionalBaseType::SourceSpaceType;
  using typename FunctionalBaseType::SourceType;
  using typename FunctionalBaseType::SourceVectorType;

  using typename WalkerBaseType::ElementType;
  using ElementFilterType = XT::Grid::ElementFilter<AssemblyGridViewType>;
  using IntersectionFilterType = XT::Grid::IntersectionFilter<AssemblyGridViewType>;
  using ApplyOnAllElements = XT::Grid::ApplyOn::AllElements<AssemblyGridViewType>;
  using ApplyOnAllIntersection = XT::Grid::ApplyOn::AllIntersections<AssemblyGridViewType>;

  /**
   * \name Ctors which accept an existing vector into which to assemble.
   * \{
   */

  VectorBasedFunctional(AssemblyGridViewType assembly_grid_view,
                        const SourceSpaceType& source_spc,
                        SourceVectorType& vec)
    : VectorStorage(vec)
    , FunctionalBaseType(source_spc, VectorStorage::access())
    , WalkerBaseType(assembly_grid_view)
    , assembled_(false)
  {
    // to detect assembly
    this->append(
        [](/*prepare nothing*/) {}, [](const auto&) { /*apply nothing*/ }, [&](/*finalize*/) { assembled_ = true; });
  }

  /**
   * \}
   * \name Ctors which create an appropriate vector into which to assemble (which is accessible via vector()).
   * \{
   */

  VectorBasedFunctional(AssemblyGridViewType assembly_grid_view, const SourceSpaceType& source_spc)
    : VectorStorage(new SourceVectorType(source_spc.mapper().size(), 0))
    , FunctionalBaseType(source_spc, VectorStorage::access())
    , WalkerBaseType(assembly_grid_view)
    , assembled_(false)
  {
    // to detect assembly
    this->append(
        [](/*prepare nothing*/) {}, [](const auto&) { /*apply nothing*/ }, [&](/*finalize*/) { assembled_ = true; });
  }

  /// \}

  VectorBasedFunctional(const ThisType&) = default;
  VectorBasedFunctional(ThisType&&) = default;

  using FunctionalBaseType::vector;

  SourceVectorType& vector()
  {
    return VectorStorage::access();
  }

  using WalkerBaseType::append;

  ThisType& append(const LocalElementFunctionalInterface<ElementType, r, rC, F, DofFieldType>& local_functional,
                   const XT::Common::Parameter& param = {},
                   const ElementFilterType& filter = ApplyOnAllElements())
  {
    LocalElementFunctionalAssembler<V, AssemblyGridViewType, r, rC, F, GV> tmp(
        this->source_space(), local_functional, VectorStorage::access(), param);
    this->append(tmp, filter);
    return *this;
  }

  ThisType& append(const LocalElementFunctionalInterface<ElementType, r, rC, F, DofFieldType>& local_functional,
                   const XT::Common::Parameter& param,
                   std::function<bool(const AssemblyGridViewType&, const ElementType&)> filter_lambda)
  {
    return append(
        local_functional, param, XT::Grid::ApplyOn::GenericFilteredElements<AssemblyGridViewType>(filter_lambda));
  }

  // similar append for LocalIntersectionFunctionalInterface ...

  void assemble(const bool use_tbb = false) override final
  {
    if (assembled_)
      return;
    // This clears all appended functionals, which is ok, since we are done after assembling once!
    this->walk(use_tbb);
    assembled_ = true;
  }

  template <class EntityRange>
  void assemble_range(const EntityRange& entity_range)
  {
    if (assembled_)
      return;
    // This clears all appended functionals, which is ok, since we are done after assembling once!
    this->walk_range(entity_range);
    assembled_ = true;
  }

  FieldType apply(const SourceType& source) const override final
  {
    DUNE_THROW_IF(!assembled_,
                  Exceptions::functional_error,
                  "You have to call assemble() first, or append this "
                  "functional to an existing XT::Grid::Walker or "
                  "GlobalAssembler!");
    return FunctionalBaseType::apply(source);
  } // ... apply(...)

private:
  bool assembled_;
}; // class VectorBasedFunctional


template <class GV, size_t r, size_t rC, class F, class V>
VectorBasedFunctional<typename XT::LA::VectorInterface<V>::derived_type, GV, r, rC, F>
make_vector_functional(const SpaceInterface<GV, r, rC, F>& space, XT::LA::VectorInterface<V>& vector)
{
  return VectorBasedFunctional<typename XT::LA::VectorInterface<V>::derived_type, GV, r, rC, F>(
      space.grid_view(), space, vector.as_imp());
}

template <class VectorType, class GV, size_t r, size_t rC, class F>
typename std::enable_if<XT::LA::is_vector<VectorType>::value, VectorBasedFunctional<VectorType, GV, r, rC, F>>::type
make_vector_functional(const SpaceInterface<GV, r, rC, F>& space)
{
  return VectorBasedFunctional<VectorType, GV, r, rC, F>(space.grid_view(), space);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_FUNCTIONALS_VECTOR_BASED_HH
