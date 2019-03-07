// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_OPERATORS_OSWALD_INTERPOLATION_HH
#define DUNE_GDT_OPERATORS_OSWALD_INTERPOLATION_HH

#include <map>
#include <set>
#include <vector>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/boundaryinfo/allneumann.hh>
#include <dune/xt/grid/boundaryinfo/interfaces.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/exceptions.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
#include <dune/gdt/tools/dirichlet-constraints.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


template <class M,
          class AssemblyGridView,
          size_t dim = 1,
          size_t dim_cols = 1,
          class SGV = AssemblyGridView,
          class RGV = AssemblyGridView>
class OswaldInterpolationOperator : public OperatorInterface<M, SGV, dim, dim_cols, dim, dim_cols, RGV>
{
  static_assert(XT::Grid::is_view<AssemblyGridView>::value, "");
  static_assert(dim == 1, "I did not think about this yet, feel free to implement!");
  static_assert(dim_cols == 1, "I did not think about this yet, feel free to implement!");
  using BaseType = OperatorInterface<M, SGV, dim, dim_cols, dim, dim_cols, RGV>;
  using ThisType = OswaldInterpolationOperator<M, AssemblyGridView, dim, dim_cols, SGV, RGV>;

public:
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  using AGV = AssemblyGridView;
  using AssemblyGridViewType = AGV;
  using I = XT::Grid::extract_intersection_t<AssemblyGridViewType>;

  /**
   * \param boundary_info To determine the Dirichlet boundary DoFs on which to set the range to zero.
   */
  OswaldInterpolationOperator(AssemblyGridViewType assembly_grid_view,
                              const SourceSpaceType& src_spc,
                              const RangeSpaceType& rng_spc,
                              const XT::Grid::BoundaryInfo<I>& boundary_info)
    : assembly_grid_view_(assembly_grid_view)
    , source_space_(src_spc)
    , range_space_(rng_spc)
    , boundary_info_(boundary_info)
    , assembled_(false)
    , global_DoF_id_to_global_LP_id_map_()
  {
    DUNE_THROW_IF(!range_space_.is_lagrangian(), Exceptions::operator_error, "This does not make any sense!");
    DUNE_THROW_IF(range_space_.continuous(0), Exceptions::operator_error, "This does not make any sense!");
    DUNE_THROW_IF(
        range_space_.min_polorder() != range_space_.max_polorder(), Exceptions::operator_error, "Not implemented yet!");
  }

  /**
   * Does not set any boundary DoFs to zero.
   */
  OswaldInterpolationOperator(const AssemblyGridViewType& assembly_grid_view,
                              const SourceSpaceType& src_spc,
                              const RangeSpaceType& rng_spc)
    : assembly_grid_view_(assembly_grid_view)
    , source_space_(src_spc)
    , range_space_(rng_spc)
    , boundary_info_(new XT::Grid::AllNeumannBoundaryInfo<I>()) // Anything without Dirichlet
    , assembled_(false)
  {
    DUNE_THROW_IF(!range_space_.is_lagrangian(), Exceptions::operator_error, "This does not make any sense!");
    DUNE_THROW_IF(range_space_.continuous(0), Exceptions::operator_error, "This does not make any sense!");
    DUNE_THROW_IF(
        range_space_.min_polorder() != range_space_.max_polorder(), Exceptions::operator_error, "Not implemented yet!");
  }

  bool linear() const override final
  {
    return true;
  }

  const SourceSpaceType& source_space() const override final
  {
    return source_space_;
  }

  const RangeSpaceType& range_space() const override final
  {
    return range_space_;
  }

  ThisType& assemble(const bool use_tbb = false) override final
  {
    if (assembled_)
      return *this;
    // create conforming space of same order to be used as mapper for the global Lagrange points
    const auto order = range_space_.min_polorder();
    DUNE_THROW_IF(
        range_space_.max_polorder() != order, Exceptions::operator_error, "Not implemented yet for variable orders!");
    auto cg_space = make_continuous_lagrange_space(assembly_grid_view_, order);
    // determine Dirichlet DoFs
    DirichletConstraints<I, decltype(cg_space)> dirichlet_constraints(boundary_info_.access(), cg_space);
    auto walker = XT::Grid::make_walker(assembly_grid_view_);
    walker.append(dirichlet_constraints);
    walker.walk(use_tbb);
    boundary_LPs_ = std::move(dirichlet_constraints.dirichlet_DoFs());
    global_LP_id_to_global_DoF_id_map_.resize(cg_space.mapper().size());
    global_DoF_id_to_global_LP_id_map_.resize(range_space_.mapper().size(), std::numeric_limits<size_t>::max());
    DynamicVector<size_t> global_lagrange_point_indices(cg_space.mapper().max_local_size());
    DynamicVector<size_t> global_DoF_indices(range_space_.mapper().max_local_size());
    // walk the grid
    for (auto&& element : elements(assembly_grid_view_)) {
      const auto& lagrange_points = cg_space.finite_elements().get(element.geometry().type(), order).lagrange_points();
      DUNE_THROW_IF(range_space_.finite_elements().get(element.geometry().type(), order).lagrange_points().size()
                        != lagrange_points.size(),
                    Exceptions::operator_error,
                    "This should not happen, the Lagrange points should coincide for Lagrange spaces of same order!\n"
                        << "range_space_.finite_element(element.geometry().type(), order).lagrange_points().size() = "
                        << range_space_.finite_elements().get(element.geometry().type(), order).lagrange_points().size()
                        << "\nlagrange_points.size() = " << lagrange_points.size());
      cg_space.mapper().global_indices(element, global_lagrange_point_indices);
      range_space_.mapper().global_indices(element, global_DoF_indices);
      for (size_t ii = 0; ii < lagrange_points.size(); ++ii) {
        const auto global_LP_id = global_lagrange_point_indices[ii];
        const auto global_DoF_id = global_DoF_indices[ii];
        global_LP_id_to_global_DoF_id_map_[global_LP_id].insert(global_DoF_id);
        global_DoF_id_to_global_LP_id_map_[global_DoF_id] = global_LP_id;
      }
    }
    assembled_ = true;
    return *this;
  } // ... assemble(...)

  using BaseType::apply;

  void
  apply(const VectorType& source, VectorType& range, const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    // some checks
    DUNE_THROW_IF(!assembled_, Exceptions::operator_error, "You need to call assemble() first!");
    DUNE_THROW_IF(!source.valid(), Exceptions::operator_error, "source contains inf or nan!");
    DUNE_THROW_IF(!source_space_.contains(source), Exceptions::operator_error, "");
    DUNE_THROW_IF(!range_space_.contains(range), Exceptions::operator_error, "");
    const auto source_function = make_discrete_function(source_space_, source);
    auto local_source = source_function.local_function();
    DynamicVector<size_t> global_DoF_indices(range_space_.mapper().max_local_size());
    // clear range on those DoFs associated with assembly_grid_view_ individually
    // (might only be a subset, range *= 0 would clear too much)
    for (const auto& DoF_id : global_DoF_id_to_global_LP_id_map_)
      if (DoF_id != std::numeric_limits<size_t>::max())
        range[DoF_id] = 0;
    // walk the grid to average on all inner Lagrange points
    auto range_basis = range_space_.basis().localize();
    for (auto&& element : elements(assembly_grid_view_)) {
      local_source->bind(element);
      range_basis->bind(element);
      const auto& lagrange_points = range_basis->finite_element().lagrange_points();
      range_space_.mapper().global_indices(element, global_DoF_indices);
      DUNE_THROW_IF(global_DoF_indices.size() < lagrange_points.size(),
                    Exceptions::operator_error,
                    "This should not happen, the range_space is broken:\n"
                        << "global_DoF_indices.size() = " << global_DoF_indices.size() << "\n"
                        << "lagrange_points.size() = " << lagrange_points.size());
      for (size_t ii = 0; ii < lagrange_points.size(); ++ii) {
        const auto& lagrange_point = lagrange_points[ii];
        const auto global_DoF_id = global_DoF_indices[ii];
        const auto global_LP_id = global_DoF_id_to_global_LP_id_map_.at(global_DoF_id);
        const auto& DoFs_per_global_LP = global_LP_id_to_global_DoF_id_map_.at(global_LP_id);
        const auto source_value = local_source->evaluate(lagrange_point)[0] / DoFs_per_global_LP.size();
        for (const auto& DoF_id : DoFs_per_global_LP) {
          range[DoF_id] += source_value;
        }
      }
    }
    // set Dirichlet DoFs to zero
    for (const auto& global_LP_id : boundary_LPs_)
      for (const auto& global_DoF_id : global_LP_id_to_global_DoF_id_map_.at(global_LP_id))
        range[global_DoF_id] = 0;
  } // ... apply(...)

private:
  const AssemblyGridViewType assembly_grid_view_;
  const SourceSpaceType& source_space_;
  const RangeSpaceType& range_space_;
  const XT::Common::ConstStorageProvider<XT::Grid::BoundaryInfo<I>> boundary_info_;
  bool assembled_;
  std::vector<size_t> global_DoF_id_to_global_LP_id_map_;
  std::vector<std::set<size_t>> global_LP_id_to_global_DoF_id_map_;
  std::set<size_t> boundary_LPs_;
}; // class OswaldInterpolationOperator


template <class MatrixType, class AssemblyGridView, class SGV, size_t dim, size_t dim_cols, class RGV>
OswaldInterpolationOperator<MatrixType, AssemblyGridView, dim, dim_cols, SGV, RGV> make_oswald_interpolation_operator(
    const AssemblyGridView& assembly_grid_view,
    const SpaceInterface<SGV, dim, dim_cols>& source_space,
    const SpaceInterface<RGV, dim, dim_cols>& range_space,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<AssemblyGridView>>& boundary_info)
{
  return OswaldInterpolationOperator<MatrixType, AssemblyGridView, dim, dim_cols, SGV, RGV>(
      assembly_grid_view, source_space, range_space, boundary_info);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_OSWALD_INTERPOLATION_HH
