// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_LOCAL_ASSEMBLER_FUNCTIONAL_ASSEMBLERS_HH
#define DUNE_GDT_LOCAL_ASSEMBLER_FUNCTIONAL_ASSEMBLERS_HH

#include <dune/xt/grid/functors/interfaces.hh>
#include <dune/xt/la/container/vector-interface.hh>

#include <dune/gdt/local/functionals/interfaces.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {


template <class Vector, class GridView, size_t r, size_t rC, class R = double, class SpaceGridView = GridView>
class LocalElementFunctionalAssembler : public XT::Grid::ElementFunctor<GridView>
{
  static_assert(XT::LA::is_vector<Vector>::value, "");
  static_assert(XT::Grid::is_view<GridView>::value, "");
  static_assert(XT::Grid::is_view<SpaceGridView>::value, "");

  using ThisType = LocalElementFunctionalAssembler<Vector, GridView, r, rC, R, SpaceGridView>;
  using BaseType = XT::Grid::ElementFunctor<GridView>;

public:
  using typename BaseType::ElementType;
  using VectorType = Vector;
  using FieldType = typename VectorType::ScalarType;
  using SpaceType = SpaceInterface<SpaceGridView, r, rC, R>;
  using LocalFunctionalType = LocalElementFunctionalInterface<ElementType, r, rC, R, FieldType>;

  LocalElementFunctionalAssembler(const SpaceType& space,
                                  const LocalFunctionalType& local_functional,
                                  VectorType& global_vector,
                                  const XT::Common::Parameter& param = {})
    : space_(space)
    , local_functional_(local_functional.copy())
    , global_vector_(global_vector)
    , param_(param)
    , local_vector_(space_.mapper().max_local_size())
    , global_indices_(space_.mapper().max_local_size())
    , basis_(space_.basis().localize())
  {
    DUNE_THROW_IF(global_vector_.size() != space_.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "global_vector.size() = " << global_vector_.size() << "\n  "
                                            << "space.mapper().size()"
                                            << space_.mapper().size());
  }

  LocalElementFunctionalAssembler(const ThisType& other)
    : BaseType()
    , space_(other.space_)
    , local_functional_(other.local_functional_->copy())
    , global_vector_(other.global_vector_)
    , param_(other.param_)
    , local_vector_(other.local_vector_)
    , global_indices_(other.global_indices_)
    , basis_(space_.basis().localize())
  {
  }

  LocalElementFunctionalAssembler(ThisType&& source) = default;

  BaseType* copy() override final
  {
    return new ThisType(*this);
  }

  void apply_local(const ElementType& element) override final
  {
    // apply functional
    basis_->bind(element);
    local_functional_->apply(*basis_, local_vector_, param_);
    // copy local vector to global
    space_.mapper().global_indices(element, global_indices_);
    for (size_t jj = 0; jj < basis_->size(param_); ++jj)
      global_vector_.add_to_entry(global_indices_[jj], local_vector_[jj]);
  }

private:
  const SpaceType& space_;
  std::unique_ptr<LocalFunctionalType> local_functional_;
  VectorType& global_vector_;
  XT::Common::Parameter param_;
  DynamicVector<FieldType> local_vector_;
  DynamicVector<size_t> global_indices_;
  mutable std::unique_ptr<typename SpaceType::GlobalBasisType::LocalizedBasisType> basis_;
}; // class LocalElementFunctionalAssembler


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_ASSEMBLER_FUNCTIONAL_ASSEMBLERS_HH
