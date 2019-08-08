// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#include "config.h"

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/exceptions.bindings.hh>

#include <dune/xt/la/container/istl.hh>
#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/dd/glued.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/indicator.hh>
#include <dune/xt/functions/interfaces/function.hh>
#include <dune/xt/functions/spe10/model1.hh>

#include <dune/gdt/functionals/vector-based.hh>
#include <dune/gdt/interpolations.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/integrands/conversion.hh>
#include <dune/gdt/local/integrands/elliptic-ipdg.hh>
#include <dune/gdt/local/integrands/elliptic.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/tools/dirichlet-constraints.hh>
#include <dune/gdt/tools/grid-quality-estimates.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/h1/continuous-flattop.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
#include <dune/gdt/spaces/l2/discontinuous-lagrange.hh>
#include <dune/gdt/spaces/l2/finite-volume.hh>

using namespace Dune;
using namespace Dune::GDT;

using G = YASP_2D_EQUIDISTANT_OFFSET;

using GP = XT::Grid::GridProvider<G>;
using GV = typename G::LeafGridView;
using E = XT::Grid::extract_entity_t<GV>;

static const constexpr size_t d = G::dimension;
using M = XT::LA::IstlRowMajorSparseMatrix<double>;
using V = XT::LA::IstlDenseVector<double>;

static const LocalEllipticIpdgIntegrands::Method ipdg_variant =
    LocalEllipticIpdgIntegrands::Method::swipdg_affine_factor;

static std::string default_space_type()
{
  return "cg_p1";
}


template <class GV>
std::unique_ptr<GDT::SpaceInterface<GV>> make_subdomain_space(GV subdomain_grid_view, const std::string& space_type)
{
  if (space_type.size() >= 4 && space_type.substr(0, 4) == "cg_p") {
    const auto order = XT::Common::from_string<int>(space_type.substr(4));
    return std::make_unique<ContinuousLagrangeSpace<GV>>(subdomain_grid_view, order);
  } else if (space_type.size() >= 4 && space_type.substr(0, 4) == "dg_p") {
    const auto order = XT::Common::from_string<int>(space_type.substr(4));
    return std::make_unique<DiscontinuousLagrangeSpace<GV>>(subdomain_grid_view, order);
  } else
    DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
               "space_type = " << space_type << "\n   has to be 'cg_pX' or 'dg_pX' for some order X!");
  return nullptr;
}


/**
 * Inherits all types and methods from the coupling intersection, but uses the macro intersection to provide a correctly
 * oriented normal.
 *
 * \attention Presumes that the coupling intersection lies exactly within the macro intersection!
 */
template <class CouplingIntersectionType, class MacroIntersectionType>
class CouplingIntersectionWithCorrectNormal : public CouplingIntersectionType
{
  using BaseType = CouplingIntersectionType;

public:
  using typename BaseType::GlobalCoordinate;
  using typename BaseType::LocalCoordinate;

  CouplingIntersectionWithCorrectNormal(const CouplingIntersectionType& coupling_intersection,
                                        const MacroIntersectionType& macro_intersection)
    : BaseType(coupling_intersection)
    , macro_intersection_(macro_intersection)
  {}

  GlobalCoordinate outerNormal(const LocalCoordinate& local) const
  {
    global_ = this->geometry().global(local);
    local_ = macro_intersection_.geometry().local(global_);
    return macro_intersection_.outerNormal(local_);
  }

  GlobalCoordinate integrationOuterNormal(const LocalCoordinate& local) const
  {
    auto normal = this->unitOuterNormal(local);
    const auto integration_element = BaseType::integrationOuterNormal(local).two_norm();
    normal *= integration_element;
    return normal;
  }

  GlobalCoordinate unitOuterNormal(const LocalCoordinate& local) const
  {
    global_ = this->geometry().global(local);
    local_ = macro_intersection_.geometry().local(global_);
    return macro_intersection_.unitOuterNormal(local_);
  }

  GlobalCoordinate centerUnitOuterNormal() const
  {
    global_ = this->geometry().center();
    local_ = macro_intersection_.geometry().local(global_);
    return macro_intersection_.unitOuterNormal(local_);
  }

private:
  const MacroIntersectionType& macro_intersection_;
  mutable GlobalCoordinate global_;
  mutable LocalCoordinate local_;
}; // class CouplingIntersectionWithCorrectNormal


namespace Dune {
namespace XT {
namespace Grid {


template <class C, class I>
struct is_intersection<CouplingIntersectionWithCorrectNormal<C, I>> : public Dune::XT::Grid::is_intersection<C>
{
  using GridType = typename Dune::XT::Grid::is_intersection<C>::GridType;
  using InsideElementType = typename Dune::XT::Grid::is_intersection<C>::InsideElementType;
  using OutsideElementType = typename Dune::XT::Grid::is_intersection<C>::OutsideElementType;
};


} // namespace Grid
} // namespace XT
} // namespace Dune


class DomainDecomposition
{
public:
  using DdGridType = XT::Grid::DD::Glued<G, G, XT::Grid::Layers::leaf>;

private:
  using GV = typename DdGridType::LocalViewType;

public:
  DomainDecomposition(const std::array<unsigned int, d> num_macro_elements_per_dim,
                      const size_t num_refinements_per_subdomain,
                      const XT::Common::FieldVector<double, d> ll,
                      const XT::Common::FieldVector<double, d> ur)
    : macro_grid_(XT::Grid::make_cube_grid<G>(ll, ur, num_macro_elements_per_dim))
    , dd_grid(macro_grid_, num_refinements_per_subdomain, false, true)
    , vtk_writer_(dd_grid)
    , local_spaces_(dd_grid.num_subdomains(), nullptr)
    , local_discrete_functions_(dd_grid.num_subdomains(), nullptr)
  {}

  void
  visualize_indicators(const std::vector<double>& indicators, const std::string& filename, const std::string& plot_name)
  {
    auto fv_space = make_finite_volume_space<GV>(macro_grid_.leaf_view());
    DUNE_THROW_IF(indicators.size() != fv_space.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "indicators.size() = " << indicators.size()
                                         << "\n   fv_space.mapper().size() = " << fv_space.mapper().size());
    make_discrete_function(fv_space, V(indicators), plot_name).visualize(filename);
  } // ... visualize_indicators(...)

  void visualize_local(const std::string& filename_prefix,
                       const size_t ss,
                       const V& vec,
                       const std::string& space_type,
                       const std::string& name)
  {
    DUNE_THROW_IF(ss >= dd_grid.num_subdomains(),
                  XT::Common::Exceptions::index_out_of_range,
                  "ss = " << ss << "\n   dd_grid.num_subdomains() = " << dd_grid.num_subdomains());
    for (auto&& macro_element : elements(dd_grid.macro_grid_view())) {
      if (dd_grid.subdomain(macro_element) == ss) { // this is the subdomain we are interested in
        auto local_space = make_subdomain_space(dd_grid.local_grid(macro_element).leaf_view(), space_type);
        auto local_discrete_function = make_discrete_function(*local_space, V(vec), name);
        local_discrete_function.visualize(filename_prefix);
        break;
      }
    }
  } // ... add_local_visualization(...) */

  void add_local_visualization(const size_t ss, const V& vec, const std::string& space_type, const std::string& name)
  {
    DUNE_THROW_IF(ss >= dd_grid.num_subdomains(),
                  XT::Common::Exceptions::index_out_of_range,
                  "ss = " << ss << "\n   dd_grid.num_subdomains() = " << dd_grid.num_subdomains());
    for (auto&& macro_element : elements(dd_grid.macro_grid_view())) {
      if (dd_grid.subdomain(macro_element) == ss) { // this is the subdomain we are interested in
        if (!local_spaces_[ss]) { // we are the first to do something here
          local_spaces_[ss] = make_subdomain_space(dd_grid.local_grid(macro_element).leaf_view(), space_type);
          local_discrete_functions_[ss] = std::make_shared<std::vector<DiscreteFunction<V, GV>>>();
        }
        // create and stash the discrete function for later
        local_discrete_functions_[ss]->emplace_back(make_discrete_function(*local_spaces_[ss], V(vec), name));
        auto visualization_adapter = std::make_shared<XT::Functions::VisualizationAdapter<GV, 1, 1, double>>(
            local_discrete_functions_[ss]->back(), name);
        vtk_writer_.addVertexData(ss, visualization_adapter);
        break;
      }
    }
  } // ... add_local_visualization(...)

  void add_local_visualization(const size_t subdomain, const XT::Functions::FunctionInterface<d>& func)
  {
    auto visualization_adapter = std::make_shared<XT::Functions::VisualizationAdapter<GV, 1, 1, double>>(
        func.template as_grid_function<typename GV::template Codim<0>::Entity>(), func.name());
    vtk_writer_.addVertexData(subdomain, visualization_adapter);
  }

  void write_visualization(const std::string& filename_prefix)
  {
    vtk_writer_.write(filename_prefix, VTK::appendedraw);
    vtk_writer_.clear();
    local_discrete_functions_ =
        std::vector<std::shared_ptr<std::vector<DiscreteFunction<V, GV>>>>(dd_grid.num_subdomains(), nullptr);
    local_spaces_ = std::vector<std::shared_ptr<SpaceInterface<GV>>>(dd_grid.num_subdomains(), nullptr);
  } // ... write_visualization(...)

private:
  XT::Grid::GridProvider<G> macro_grid_;

public:
  DdGridType dd_grid;

private:
  XT::Grid::DD::GluedVTKWriter<G, G, XT::Grid::Layers::leaf> vtk_writer_;
  std::vector<std::shared_ptr<SpaceInterface<GV>>> local_spaces_;
  std::vector<std::shared_ptr<std::vector<DiscreteFunction<V, GV>>>> local_discrete_functions_;
}; // class DomainDecomposition


class PartitionOfUnityBase
{
public:
  PartitionOfUnityBase(const DomainDecomposition& dd, SpaceInterface<GV>*&& space)
    : dd_(dd)
    , space_(std::move(space))
    , patches_(space_->mapper().size())
    , local_indices_(dd_.dd_grid.num_subdomains())
  {
    DynamicVector<size_t> global_indices(space_->mapper().max_local_size());
    for (auto&& subdomain : elements(space_->grid_view())) {
      space_->mapper().global_indices(subdomain, global_indices);
      const auto ss = dd_.dd_grid.subdomain(subdomain);
      for (size_t ii = 0; ii < space_->mapper().local_size(subdomain); ++ii) {
        const auto II = global_indices[ii];
        patches_[II].insert(ss);
        local_indices_[ss][II] = ii;
      }
      // build inverse index mapping
    }
  } // ... ContinuousLagrangePartitionOfUnity(...)

  size_t size() const
  {
    return space_->mapper().size();
  }

  std::set<size_t> patch(const size_t ii) const
  {
    DUNE_THROW_IF(
        ii >= size(), XT::Common::Exceptions::index_out_of_range, "ii = " << ii << "\n   size() = " << size());
    return patches_.at(ii);
  }

  std::vector<V> on_subdomain(const size_t ss, const std::string space_type)
  {
    auto coarse_basis = space_->basis().localize();
    std::vector<V> interpolated_basis;
    for (auto&& macro_element : elements(dd_.dd_grid.macro_grid_view())) {
      if (dd_.dd_grid.subdomain(macro_element) == ss) {
        // this is the subdomain we are interested in, create space
        auto subdomain_grid_view = dd_.dd_grid.local_grid(macro_element).leaf_view();
        auto subdomain_space = make_subdomain_space(subdomain_grid_view, space_type);
        coarse_basis->bind(macro_element);
        for (size_t ii = 0; ii < coarse_basis->size(); ++ii)
          interpolated_basis.push_back(
              default_interpolation<V>(coarse_basis->order(),
                                       [&](const auto& point_in_physical_coordinates, const auto&) {
                                         const auto point_macro_reference_element =
                                             macro_element.geometry().local(point_in_physical_coordinates);
                                         return coarse_basis->evaluate_set(point_macro_reference_element)[ii];
                                       },
                                       *subdomain_space)
                  .dofs()
                  .vector());
        break;
      }
    }
    DUNE_THROW_IF(interpolated_basis.size() == 0, InvalidStateException, "This should not happen, ss = " << ss);
    return interpolated_basis;
  } // ... on_subdomain(...)

  size_t local_index(const size_t ii, const size_t ss)
  {
    DUNE_THROW_IF(
        ii >= size(), XT::Common::Exceptions::index_out_of_range, "ii = " << ii << "\n   size() = " << size());
    DUNE_THROW_IF(ss >= local_indices_.size(),
                  XT::Common::Exceptions::index_out_of_range,
                  "ss = " << ss << "\n   local_indices_.size() = " << local_indices_.size());
    auto search_result = local_indices_[ss].find(ii);
    DUNE_THROW_IF(search_result == local_indices_[ss].end(),
                  XT::Common::Exceptions::index_out_of_range,
                  "Global index " << ii << " is not associated with subdomain " << ss << "!");
    return search_result->second;
  } // ... local_index(...)

protected:
  const DomainDecomposition& dd_;
  const std::unique_ptr<SpaceInterface<GV>> space_;
  std::vector<std::set<size_t>> patches_;
  std::vector<std::map<size_t, size_t>> local_indices_;
}; // class PartitionOfUnityBase


class ContinuousLagrangePartitionOfUnity : public PartitionOfUnityBase
{
public:
  ContinuousLagrangePartitionOfUnity(const DomainDecomposition& dd)
    : PartitionOfUnityBase(dd, new ContinuousLagrangeSpace<GV>(dd.dd_grid.macro_grid_view(), /*polorder=*/1))
  {}
};


class ContinuousFlatTopPartitionOfUnity : public PartitionOfUnityBase
{
public:
  ContinuousFlatTopPartitionOfUnity(const DomainDecomposition& dd, const double overlap = 0.5)
    : PartitionOfUnityBase(dd, new ContinuousFlatTopSpace<GV>(dd.dd_grid.macro_grid_view(), /*polorder=*/1, overlap))
  {}
};


template <class MacroGV, class MicroGV>
class MacroGridBasedBoundaryInfo : public XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<MicroGV>>
{
  using BaseType = XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<MicroGV>>;

public:
  using typename BaseType::IntersectionType;
  using MacroBoundaryInfoType = XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<MacroGV>>;
  using MacroElementType = XT::Grid::extract_entity_t<MacroGV>;

  MacroGridBasedBoundaryInfo(const MacroGV& macro_grid_view,
                             const MacroElementType& macro_element,
                             const MacroBoundaryInfoType& macro_boundary_info)
    : macro_grid_view_(macro_grid_view)
    , macro_element_(macro_element)
    , macro_boundary_info_(macro_boundary_info)
  {}

  const XT::Grid::BoundaryType& type(const IntersectionType& intersection) const override final
  {
    // find out if this micro intersection lies within the macro element or on a macro intersection
    for (auto&& macro_intersection : intersections(macro_grid_view_, macro_element_)) {
      const size_t num_corners = intersection.geometry().corners();
      size_t num_corners_inside = 0;
      size_t num_corners_outside = 0;
      for (size_t cc = 0; cc < num_corners; ++cc) {
        const auto micro_corner = intersection.geometry().corner(cc);
        if (XT::Grid::contains(macro_intersection, micro_corner))
          ++num_corners_inside;
        else
          ++num_corners_outside;
      }
      if (num_corners_inside == num_corners && num_corners_outside == 0) {
        // we found the macro intersection this micro intersection belongs to
        return macro_boundary_info_.type(macro_intersection);
      }
    }
    // we could not find a macro intersection this micro intersection belongs to
    return no_boundary_;
  } // ... type(...)

  const MacroGV& macro_grid_view_;
  const MacroElementType& macro_element_;
  const MacroBoundaryInfoType& macro_boundary_info_;
  const XT::Grid::NoBoundary no_boundary_;
}; // class MacroGridBasedBoundaryInfo


std::unique_ptr<M> assemble_local_system_matrix(const XT::Functions::GridFunctionInterface<E>& diffusion_factor,
                                                const XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
                                                const DomainDecomposition& domain_decomposition,
                                                const size_t ss,
                                                const std::string space_type)
{
  DUNE_THROW_IF(ss >= domain_decomposition.dd_grid.num_subdomains(),
                XT::Common::Exceptions::index_out_of_range,
                "ss = " << ss << "\n   domain_decomposition.dd_grid.num_subdomains() = "
                        << domain_decomposition.dd_grid.num_subdomains());
  std::unique_ptr<M> subdomain_matrix;
  bool found_subdomain = false;
  for (auto&& macro_element : elements(domain_decomposition.dd_grid.macro_grid_view())) {
    if (domain_decomposition.dd_grid.subdomain(macro_element) == ss) {
      // this is the subdomain we are interested in, create space
      found_subdomain = true;
      auto subdomain_grid_view = domain_decomposition.dd_grid.local_grid(macro_element).leaf_view();
      using GV = decltype(subdomain_grid_view);
      using E = typename GV::template Codim<0>::Entity;
      using I = typename GV::Intersection;
      auto subdomain_space = make_subdomain_space(subdomain_grid_view, space_type);
      // create operator
      auto subdomain_operator = make_matrix_operator<M>(*subdomain_space, Stencil::element_and_intersection);
      const LocalElementIntegralBilinearForm<E> element_bilinear_form(
          LocalEllipticIntegrand<E>(diffusion_factor, diffusion_tensor));
      subdomain_operator.append(element_bilinear_form);
      if (!subdomain_space->continuous(0)) {
        const LocalIntersectionIntegralBilinearForm<I> coupling_bilinear_form(
            LocalEllipticIpdgIntegrands::Inner<I, double, ipdg_variant>(diffusion_factor, diffusion_tensor));
        subdomain_operator.append(coupling_bilinear_form, {}, XT::Grid::ApplyOn::InnerIntersectionsOnce<GV>());
      }
      subdomain_operator.assemble();
      subdomain_matrix = std::make_unique<M>(subdomain_operator.matrix());
      break;
    }
  }
  DUNE_THROW_IF(!found_subdomain, InvalidStateException, "ss = " << ss);
  return std::move(subdomain_matrix);
} // ... assemble_local_system_matrix(...)


std::unique_ptr<V> assemble_local_rhs(const XT::Functions::GridFunctionInterface<E>& force,
                                      const DomainDecomposition& domain_decomposition,
                                      const size_t ss,
                                      const std::string space_type)
{
  DUNE_THROW_IF(ss >= domain_decomposition.dd_grid.num_subdomains(),
                XT::Common::Exceptions::index_out_of_range,
                "ss = " << ss << "\n   domain_decomposition.dd_grid.num_subdomains() = "
                        << domain_decomposition.dd_grid.num_subdomains());
  std::unique_ptr<V> subdomain_vector;
  bool found_subdomain = false;
  for (auto&& macro_element : elements(domain_decomposition.dd_grid.macro_grid_view())) {
    if (domain_decomposition.dd_grid.subdomain(macro_element) == ss) {
      // this is the subdomain we are interested in, create space
      found_subdomain = true;
      auto subdomain_grid_view = domain_decomposition.dd_grid.local_grid(macro_element).leaf_view();
      using GV = decltype(subdomain_grid_view);
      using E = typename GV::template Codim<0>::Entity;
      auto subdomain_space = make_subdomain_space(subdomain_grid_view, space_type);
      // create functional
      auto subdomain_functional = make_vector_functional<V>(*subdomain_space);
      const LocalElementIntegralFunctional<E> element_functional(
          local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E>(), force));
      subdomain_functional.append(element_functional);
      subdomain_functional.assemble();
      subdomain_vector = std::make_unique<V>(subdomain_functional.vector());
      break;
    }
  }
  DUNE_THROW_IF(!found_subdomain, InvalidStateException, "ss = " << ss);
  return std::move(subdomain_vector);
} // ... assemble_local_rhs(...)


std::tuple<std::unique_ptr<M>, std::unique_ptr<M>, std::unique_ptr<M>, std::unique_ptr<M>>
assemble_coupling_matrices(const XT::Functions::GridFunctionInterface<E>& diffusion_factor,
                           const XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
                           DomainDecomposition& domain_decomposition,
                           const size_t ss,
                           const size_t nn,
                           const std::string space_type)
{
  DUNE_THROW_IF(ss >= domain_decomposition.dd_grid.num_subdomains(),
                XT::Common::Exceptions::index_out_of_range,
                "ss = " << ss << "\n   domain_decomposition.dd_grid.num_subdomains() = "
                        << domain_decomposition.dd_grid.num_subdomains());
  std::unique_ptr<M> coupling_matrix_in_in;
  std::unique_ptr<M> coupling_matrix_in_out;
  std::unique_ptr<M> coupling_matrix_out_in;
  std::unique_ptr<M> coupling_matrix_out_out;
  for (auto&& inside_macro_element : elements(domain_decomposition.dd_grid.macro_grid_view())) {
    if (domain_decomposition.dd_grid.subdomain(inside_macro_element) == ss) {
      // this is the subdomain we are interested in, create space
      auto inner_subdomain_grid_view = domain_decomposition.dd_grid.local_grid(inside_macro_element).leaf_view();
      auto inner_subdomain_space = make_subdomain_space(inner_subdomain_grid_view, space_type);
      bool found_correct_macro_intersection = false;
      for (auto&& macro_intersection :
           intersections(domain_decomposition.dd_grid.macro_grid_view(), inside_macro_element)) {
        if (macro_intersection.neighbor()) {
          const auto outside_macro_element = macro_intersection.outside();
          if (domain_decomposition.dd_grid.subdomain(outside_macro_element) == nn) {
            found_correct_macro_intersection = true;
            // these are the subdomains we are interested in
            auto outer_subdomain_grid_view = domain_decomposition.dd_grid.local_grid(outside_macro_element).leaf_view();
            auto outer_subdomain_space = make_subdomain_space(outer_subdomain_grid_view, space_type);
            // walk the coupling once to create the sparsity patterns
            XT::LA::SparsityPatternDefault pattern_in_in(inner_subdomain_space->mapper().size());
            XT::LA::SparsityPatternDefault pattern_in_out(inner_subdomain_space->mapper().size());
            XT::LA::SparsityPatternDefault pattern_out_in(outer_subdomain_space->mapper().size());
            XT::LA::SparsityPatternDefault pattern_out_out(outer_subdomain_space->mapper().size());
            const auto& coupling =
                domain_decomposition.dd_grid.coupling(inside_macro_element, -1, outside_macro_element, -1, true);
            DynamicVector<size_t> global_indices_in(inner_subdomain_space->mapper().max_local_size());
            DynamicVector<size_t> global_indices_out(outer_subdomain_space->mapper().max_local_size());
            auto inside_basis = inner_subdomain_space->basis().localize();
            auto outside_basis = outer_subdomain_space->basis().localize();
            const auto coupling_intersection_it_end = coupling.template iend<0>();
            for (auto coupling_intersection_it = coupling.template ibegin<0>();
                 coupling_intersection_it != coupling_intersection_it_end;
                 ++coupling_intersection_it) {
              const auto& coupling_intersection = *coupling_intersection_it;
              const auto inside_element = coupling_intersection.inside();
              const auto outside_element = coupling_intersection.outside();
              inside_basis->bind(inside_element);
              outside_basis->bind(outside_element);
              inner_subdomain_space->mapper().global_indices(inside_element, global_indices_in);
              outer_subdomain_space->mapper().global_indices(outside_element, global_indices_out);
              for (size_t ii = 0; ii < inside_basis->size(); ++ii) {
                for (size_t jj = 0; jj < inside_basis->size(); ++jj)
                  pattern_in_in.insert(global_indices_in[ii], global_indices_in[jj]);
                for (size_t jj = 0; jj < outside_basis->size(); ++jj)
                  pattern_in_out.insert(global_indices_in[ii], global_indices_out[jj]);
              }
              for (size_t ii = 0; ii < outside_basis->size(); ++ii) {
                for (size_t jj = 0; jj < inside_basis->size(); ++jj)
                  pattern_out_in.insert(global_indices_out[ii], global_indices_in[jj]);
                for (size_t jj = 0; jj < outside_basis->size(); ++jj)
                  pattern_out_out.insert(global_indices_out[ii], global_indices_out[jj]);
              }
            }
            // we need to ensure at least one entry per row
            for (size_t ii = 0; ii < inner_subdomain_space->mapper().size(); ++ii) {
              pattern_in_in.insert(ii, 0);
              pattern_in_out.insert(ii, 0);
            }
            for (size_t ii = 0; ii < outer_subdomain_space->mapper().size(); ++ii) {
              pattern_out_in.insert(ii, 0);
              pattern_out_out.insert(ii, 0);
            }
            pattern_in_in.sort();
            pattern_in_out.sort();
            pattern_out_in.sort();
            pattern_out_out.sort();
            coupling_matrix_in_in = std::make_unique<M>(
                inner_subdomain_space->mapper().size(), inner_subdomain_space->mapper().size(), pattern_in_in);
            coupling_matrix_in_out = std::make_unique<M>(
                inner_subdomain_space->mapper().size(), outer_subdomain_space->mapper().size(), pattern_in_out);
            coupling_matrix_out_in = std::make_unique<M>(
                outer_subdomain_space->mapper().size(), inner_subdomain_space->mapper().size(), pattern_out_in);
            coupling_matrix_out_out = std::make_unique<M>(
                outer_subdomain_space->mapper().size(), outer_subdomain_space->mapper().size(), pattern_out_out);
            // walk the coupling for the second time to assemble
            DynamicMatrix<double> local_matrix_in_in(inner_subdomain_space->mapper().max_local_size(),
                                                     inner_subdomain_space->mapper().max_local_size());
            DynamicMatrix<double> local_matrix_in_out(inner_subdomain_space->mapper().max_local_size(),
                                                      outer_subdomain_space->mapper().max_local_size());
            DynamicMatrix<double> local_matrix_out_in(outer_subdomain_space->mapper().max_local_size(),
                                                      inner_subdomain_space->mapper().max_local_size());
            DynamicMatrix<double> local_matrix_out_out(outer_subdomain_space->mapper().max_local_size(),
                                                       outer_subdomain_space->mapper().max_local_size());
            using MacroI = decltype(macro_intersection);
            using CouplingI = typename DomainDecomposition::DdGridType::GlueType::Intersection;
            using I = CouplingIntersectionWithCorrectNormal<CouplingI, MacroI>;
            using E = typename I::InsideEntity;
            const LocalIntersectionIntegralBilinearForm<I> intersection_bilinear_form(
                LocalEllipticIpdgIntegrands::Inner<I, double, ipdg_variant>(diffusion_factor, diffusion_tensor));
            for (auto coupling_intersection_it = coupling.template ibegin<0>();
                 coupling_intersection_it != coupling_intersection_it_end;
                 ++coupling_intersection_it) {
              const I coupling_intersection(*coupling_intersection_it, macro_intersection);
              const auto inside_element = coupling_intersection.inside();
              const auto outside_element = coupling_intersection.outside();
              inside_basis->bind(inside_element);
              outside_basis->bind(outside_element);
              inner_subdomain_space->mapper().global_indices(inside_element, global_indices_in);
              outer_subdomain_space->mapper().global_indices(outside_element, global_indices_out);
              intersection_bilinear_form.apply2(coupling_intersection,
                                                *inside_basis,
                                                *inside_basis,
                                                *outside_basis,
                                                *outside_basis,
                                                local_matrix_in_in,
                                                local_matrix_in_out,
                                                local_matrix_out_in,
                                                local_matrix_out_out);
              for (size_t ii = 0; ii < inside_basis->size(); ++ii) {
                for (size_t jj = 0; jj < inside_basis->size(); ++jj)
                  coupling_matrix_in_in->add_to_entry(
                      global_indices_in[ii], global_indices_in[jj], local_matrix_in_in[ii][jj]);
                for (size_t jj = 0; jj < outside_basis->size(); ++jj)
                  coupling_matrix_in_out->add_to_entry(
                      global_indices_in[ii], global_indices_out[jj], local_matrix_in_out[ii][jj]);
              }
              for (size_t ii = 0; ii < outside_basis->size(); ++ii) {
                for (size_t jj = 0; jj < inside_basis->size(); ++jj)
                  coupling_matrix_out_in->add_to_entry(
                      global_indices_out[ii], global_indices_in[jj], local_matrix_out_in[ii][jj]);
                for (size_t jj = 0; jj < outside_basis->size(); ++jj)
                  coupling_matrix_out_out->add_to_entry(
                      global_indices_out[ii], global_indices_out[jj], local_matrix_out_out[ii][jj]);
              }
            }
            break;
          }
        }
      }
      DUNE_THROW_IF(!found_correct_macro_intersection,
                    XT::Common::Exceptions::index_out_of_range,
                    "ss = " << ss << "\n   nn = " << nn);
    }
  }
  return std::make_tuple(std::move(coupling_matrix_in_in),
                         std::move(coupling_matrix_in_out),
                         std::move(coupling_matrix_out_in),
                         std::move(coupling_matrix_out_out));
} // ... assemble_coupling_matrices(...)


std::unique_ptr<M> assemble_boundary_matrix(const XT::Functions::GridFunctionInterface<E>& diffusion_factor,
                                            const XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
                                            DomainDecomposition& domain_decomposition,
                                            const size_t ss,
                                            const std::string space_type)
{
  DUNE_THROW_IF(ss >= domain_decomposition.dd_grid.num_subdomains(),
                XT::Common::Exceptions::index_out_of_range,
                "ss = " << ss << "\n   domain_decomposition.dd_grid.num_subdomains() = "
                        << domain_decomposition.dd_grid.num_subdomains());
  using MGV = typename DomainDecomposition::DdGridType::MacroGridViewType;
  const XT::Grid::AllDirichletBoundaryInfo<XT::Grid::extract_intersection_t<MGV>> macro_boundary_info;
  std::unique_ptr<M> subdomain_matrix;
  for (auto&& macro_element : elements(domain_decomposition.dd_grid.macro_grid_view())) {
    if (domain_decomposition.dd_grid.subdomain(macro_element) == ss) {
      // this is the subdomain we are interested in, create space
      auto subdomain_grid_view = domain_decomposition.dd_grid.local_grid(macro_element).leaf_view();
      using GV = decltype(subdomain_grid_view);
      using E = typename GV::template Codim<0>::Entity;
      using I = typename GV::Intersection;
      auto subdomain_space = make_subdomain_space(subdomain_grid_view, space_type);
      const MacroGridBasedBoundaryInfo<MGV, GV> subdomain_boundary_info(
          domain_decomposition.dd_grid.macro_grid_view(), macro_element, macro_boundary_info);
      // create operator
      auto subdomain_operator = make_matrix_operator<M>(*subdomain_space, Stencil::element);
      const LocalIntersectionIntegralBilinearForm<I> dirichlet_bilinear_form(
          LocalEllipticIpdgIntegrands::DirichletBoundaryLhs<I, double, ipdg_variant>(diffusion_factor,
                                                                                     diffusion_tensor));
      subdomain_operator.append(dirichlet_bilinear_form,
                                {},
                                XT::Grid::ApplyOn::CustomBoundaryIntersections<GV>(subdomain_boundary_info,
                                                                                   new XT::Grid::DirichletBoundary()));
      subdomain_operator.assemble();
      subdomain_matrix = std::make_unique<M>(subdomain_operator.matrix());
      break;
    }
  }
  return std::move(subdomain_matrix);
} // ... assemble_boundary_matrix(...)


std::unique_ptr<M>
assemble_local_product_contributions(const XT::Functions::GridFunctionInterface<E>& diffusion_factor,
                                     const XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
                                     DomainDecomposition& domain_decomposition,
                                     const size_t ss,
                                     const std::string space_type)
{
  DUNE_THROW_IF(ss >= domain_decomposition.dd_grid.num_subdomains(),
                XT::Common::Exceptions::index_out_of_range,
                "ss = " << ss << "\n   domain_decomposition.dd_grid.num_subdomains() = "
                        << domain_decomposition.dd_grid.num_subdomains());
  std::unique_ptr<M> subdomain_matrix;
  for (auto&& macro_element : elements(domain_decomposition.dd_grid.macro_grid_view())) {
    if (domain_decomposition.dd_grid.subdomain(macro_element) == ss) {
      // this is the subdomain we are interested in, create space
      auto subdomain_grid_view = domain_decomposition.dd_grid.local_grid(macro_element).leaf_view();
      using GV = decltype(subdomain_grid_view);
      using E = typename GV::template Codim<0>::Entity;
      using I = typename GV::Intersection;
      auto subdomain_space = make_subdomain_space(subdomain_grid_view, space_type);
      // create operator
      auto subdomain_operator = make_matrix_operator<M>(*subdomain_space, Stencil::element_and_intersection);
      const LocalElementIntegralBilinearForm<E> element_bilinear_form(
          LocalEllipticIntegrand<E>(diffusion_factor, diffusion_tensor));
      subdomain_operator.append(element_bilinear_form);
      if (!subdomain_space->continuous(0)) {
        const LocalIntersectionIntegralBilinearForm<I> coupling_bilinear_form(
            LocalEllipticIpdgIntegrands::InnerOnlyPenalty<I, double, ipdg_variant>(diffusion_factor, diffusion_tensor));
        subdomain_operator.append(coupling_bilinear_form, {}, XT::Grid::ApplyOn::InnerIntersectionsOnce<GV>());
      }
      subdomain_operator.assemble();
      subdomain_matrix = std::make_unique<M>(subdomain_operator.matrix());
      break;
    }
  }
  return std::move(subdomain_matrix);
} // ... assemble_local_product_contributions(...)

std::tuple<std::unique_ptr<M>, std::unique_ptr<M>, std::unique_ptr<M>, std::unique_ptr<M>>
assemble_coupling_product_contributions(const XT::Functions::GridFunctionInterface<E>& diffusion_factor,
                                        const XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
                                        DomainDecomposition& domain_decomposition,
                                        const size_t ss,
                                        const size_t nn,
                                        const std::string space_type)
{
  DUNE_THROW_IF(ss >= domain_decomposition.dd_grid.num_subdomains(),
                XT::Common::Exceptions::index_out_of_range,
                "ss = " << ss << "\n   domain_decomposition.dd_grid.num_subdomains() = "
                        << domain_decomposition.dd_grid.num_subdomains());
  std::unique_ptr<M> coupling_matrix_in_in;
  std::unique_ptr<M> coupling_matrix_in_out;
  std::unique_ptr<M> coupling_matrix_out_in;
  std::unique_ptr<M> coupling_matrix_out_out;
  for (auto&& inside_macro_element : elements(domain_decomposition.dd_grid.macro_grid_view())) {
    if (domain_decomposition.dd_grid.subdomain(inside_macro_element) == ss) {
      // this is the subdomain we are interested in, create space
      auto inner_subdomain_grid_view = domain_decomposition.dd_grid.local_grid(inside_macro_element).leaf_view();
      auto inner_subdomain_space = make_subdomain_space(inner_subdomain_grid_view, space_type);
      bool found_correct_macro_intersection = false;
      for (auto&& macro_intersection :
           intersections(domain_decomposition.dd_grid.macro_grid_view(), inside_macro_element)) {
        if (macro_intersection.neighbor()) {
          const auto outside_macro_element = macro_intersection.outside();
          if (domain_decomposition.dd_grid.subdomain(outside_macro_element) == nn) {
            found_correct_macro_intersection = true;
            // these are the subdomains we are interested in
            auto outer_subdomain_grid_view = domain_decomposition.dd_grid.local_grid(outside_macro_element).leaf_view();
            auto outer_subdomain_space = make_subdomain_space(outer_subdomain_grid_view, space_type);
            // walk the coupling once to create the sparsity patterns
            XT::LA::SparsityPatternDefault pattern_in_in(inner_subdomain_space->mapper().size());
            XT::LA::SparsityPatternDefault pattern_in_out(inner_subdomain_space->mapper().size());
            XT::LA::SparsityPatternDefault pattern_out_in(outer_subdomain_space->mapper().size());
            XT::LA::SparsityPatternDefault pattern_out_out(outer_subdomain_space->mapper().size());
            const auto& coupling =
                domain_decomposition.dd_grid.coupling(inside_macro_element, -1, outside_macro_element, -1, true);
            DynamicVector<size_t> global_indices_in(inner_subdomain_space->mapper().max_local_size());
            DynamicVector<size_t> global_indices_out(outer_subdomain_space->mapper().max_local_size());
            auto inside_basis = inner_subdomain_space->basis().localize();
            auto outside_basis = outer_subdomain_space->basis().localize();
            const auto coupling_intersection_it_end = coupling.template iend<0>();
            for (auto coupling_intersection_it = coupling.template ibegin<0>();
                 coupling_intersection_it != coupling_intersection_it_end;
                 ++coupling_intersection_it) {
              const auto& coupling_intersection = *coupling_intersection_it;
              const auto inside_element = coupling_intersection.inside();
              const auto outside_element = coupling_intersection.outside();
              inside_basis->bind(inside_element);
              outside_basis->bind(outside_element);
              inner_subdomain_space->mapper().global_indices(inside_element, global_indices_in);
              outer_subdomain_space->mapper().global_indices(outside_element, global_indices_out);
              for (size_t ii = 0; ii < inside_basis->size(); ++ii) {
                for (size_t jj = 0; jj < inside_basis->size(); ++jj)
                  pattern_in_in.insert(global_indices_in[ii], global_indices_in[jj]);
                for (size_t jj = 0; jj < outside_basis->size(); ++jj)
                  pattern_in_out.insert(global_indices_in[ii], global_indices_out[jj]);
              }
              for (size_t ii = 0; ii < outside_basis->size(); ++ii) {
                for (size_t jj = 0; jj < inside_basis->size(); ++jj)
                  pattern_out_in.insert(global_indices_out[ii], global_indices_in[jj]);
                for (size_t jj = 0; jj < outside_basis->size(); ++jj)
                  pattern_out_out.insert(global_indices_out[ii], global_indices_out[jj]);
              }
            }
            // we need to ensure at least one entry per row
            for (size_t ii = 0; ii < inner_subdomain_space->mapper().size(); ++ii) {
              pattern_in_in.insert(ii, 0);
              pattern_in_out.insert(ii, 0);
            }
            for (size_t ii = 0; ii < outer_subdomain_space->mapper().size(); ++ii) {
              pattern_out_in.insert(ii, 0);
              pattern_out_out.insert(ii, 0);
            }
            pattern_in_in.sort();
            pattern_in_out.sort();
            pattern_out_in.sort();
            pattern_out_out.sort();
            coupling_matrix_in_in = std::make_unique<M>(
                inner_subdomain_space->mapper().size(), inner_subdomain_space->mapper().size(), pattern_in_in);
            coupling_matrix_in_out = std::make_unique<M>(
                inner_subdomain_space->mapper().size(), outer_subdomain_space->mapper().size(), pattern_in_out);
            coupling_matrix_out_in = std::make_unique<M>(
                outer_subdomain_space->mapper().size(), inner_subdomain_space->mapper().size(), pattern_out_in);
            coupling_matrix_out_out = std::make_unique<M>(
                outer_subdomain_space->mapper().size(), outer_subdomain_space->mapper().size(), pattern_out_out);
            // walk the coupling for the second time to assemble
            DynamicMatrix<double> local_matrix_in_in(inner_subdomain_space->mapper().max_local_size(),
                                                     inner_subdomain_space->mapper().max_local_size());
            DynamicMatrix<double> local_matrix_in_out(inner_subdomain_space->mapper().max_local_size(),
                                                      outer_subdomain_space->mapper().max_local_size());
            DynamicMatrix<double> local_matrix_out_in(outer_subdomain_space->mapper().max_local_size(),
                                                      inner_subdomain_space->mapper().max_local_size());
            DynamicMatrix<double> local_matrix_out_out(outer_subdomain_space->mapper().max_local_size(),
                                                       outer_subdomain_space->mapper().max_local_size());
            using I = typename DomainDecomposition::DdGridType::GlueType::Intersection;
            using E = typename I::InsideEntity;
            const LocalIntersectionIntegralBilinearForm<I> intersection_bilinear_form(
                LocalEllipticIpdgIntegrands::InnerOnlyPenalty<I, double, ipdg_variant>(diffusion_factor,
                                                                                       diffusion_tensor));
            for (auto coupling_intersection_it = coupling.template ibegin<0>();
                 coupling_intersection_it != coupling_intersection_it_end;
                 ++coupling_intersection_it) {
              const auto& coupling_intersection = *coupling_intersection_it;
              const auto inside_element = coupling_intersection.inside();
              const auto outside_element = coupling_intersection.outside();
              inside_basis->bind(inside_element);
              outside_basis->bind(outside_element);
              inner_subdomain_space->mapper().global_indices(inside_element, global_indices_in);
              outer_subdomain_space->mapper().global_indices(outside_element, global_indices_out);
              intersection_bilinear_form.apply2(coupling_intersection,
                                                *inside_basis,
                                                *inside_basis,
                                                *outside_basis,
                                                *outside_basis,
                                                local_matrix_in_in,
                                                local_matrix_in_out,
                                                local_matrix_out_in,
                                                local_matrix_out_out);
              for (size_t ii = 0; ii < inside_basis->size(); ++ii) {
                for (size_t jj = 0; jj < inside_basis->size(); ++jj)
                  coupling_matrix_in_in->add_to_entry(
                      global_indices_in[ii], global_indices_in[jj], local_matrix_in_in[ii][jj]);
                for (size_t jj = 0; jj < outside_basis->size(); ++jj)
                  coupling_matrix_in_out->add_to_entry(
                      global_indices_in[ii], global_indices_out[jj], local_matrix_in_out[ii][jj]);
              }
              for (size_t ii = 0; ii < outside_basis->size(); ++ii) {
                for (size_t jj = 0; jj < inside_basis->size(); ++jj)
                  coupling_matrix_out_in->add_to_entry(
                      global_indices_out[ii], global_indices_in[jj], local_matrix_out_in[ii][jj]);
                for (size_t jj = 0; jj < outside_basis->size(); ++jj)
                  coupling_matrix_out_out->add_to_entry(
                      global_indices_out[ii], global_indices_out[jj], local_matrix_out_out[ii][jj]);
              }
            }
            break;
          }
        }
      }
      DUNE_THROW_IF(!found_correct_macro_intersection,
                    XT::Common::Exceptions::index_out_of_range,
                    "ss = " << ss << "\n   nn = " << nn);
    }
  }
  return std::make_tuple(std::move(coupling_matrix_in_in),
                         std::move(coupling_matrix_in_out),
                         std::move(coupling_matrix_out_in),
                         std::move(coupling_matrix_out_out));
} // ... assemble_coupling_product_contributions(...)

std::unique_ptr<M>
assemble_boundary_product_contributions(const XT::Functions::GridFunctionInterface<E>& diffusion_factor,
                                        const XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
                                        DomainDecomposition& domain_decomposition,
                                        const size_t ss,
                                        const std::string space_type)
{
  DUNE_THROW_IF(ss >= domain_decomposition.dd_grid.num_subdomains(),
                XT::Common::Exceptions::index_out_of_range,
                "ss = " << ss << "\n   domain_decomposition.dd_grid.num_subdomains() = "
                        << domain_decomposition.dd_grid.num_subdomains());
  using MGV = typename DomainDecomposition::DdGridType::MacroGridViewType;
  const XT::Grid::AllDirichletBoundaryInfo<XT::Grid::extract_intersection_t<MGV>> macro_boundary_info;
  std::unique_ptr<M> subdomain_matrix;
  for (auto&& macro_element : elements(domain_decomposition.dd_grid.macro_grid_view())) {
    if (domain_decomposition.dd_grid.subdomain(macro_element) == ss) {
      // this is the subdomain we are interested in, create space
      auto subdomain_grid_view = domain_decomposition.dd_grid.local_grid(macro_element).leaf_view();
      using GV = decltype(subdomain_grid_view);
      using E = typename GV::template Codim<0>::Entity;
      using I = typename GV::Intersection;
      auto subdomain_space = make_subdomain_space(subdomain_grid_view, space_type);
      const MacroGridBasedBoundaryInfo<MGV, GV> subdomain_boundary_info(
          domain_decomposition.dd_grid.macro_grid_view(), macro_element, macro_boundary_info);
      // create operator
      auto subdomain_operator = make_matrix_operator<M>(*subdomain_space, Stencil::element);
      const LocalIntersectionIntegralBilinearForm<I> dirichlet_bilinear_form(
          LocalEllipticIpdgIntegrands::DirichletBoundaryLhsOnlyPenalty<I, double, ipdg_variant>(diffusion_factor,
                                                                                                diffusion_tensor));
      subdomain_operator.append(dirichlet_bilinear_form,
                                {},
                                XT::Grid::ApplyOn::CustomBoundaryIntersections<GV>(subdomain_boundary_info,
                                                                                   new XT::Grid::DirichletBoundary()));
      subdomain_operator.assemble();
      subdomain_matrix = std::make_unique<M>(subdomain_operator.matrix());
      break;
    }
  }
  return std::move(subdomain_matrix);
} // ... assemble_boundary_product_contributions(...)


PYBIND11_MODULE(usercode, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  Dune::XT::Common::bindings::add_initialization(m, "dune.gdt", "usercode");

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");
  py::module::import("dune.gdt.discretefunction");

  // these might already be defined
  XT::Common::bindings::try_register(m, [](auto& m_) {
    py::class_<XT::Functions::FunctionAsGridFunctionWrapper<E, 1, 1, double>,
               XT::Functions::GridFunctionInterface<E, 1, 1, double>>
        scalar_wrapper(m_, "ScalarFunctionAsGridFunctionWrapper", "ScalarFunctionAsGridFunctionWrapper");
    scalar_wrapper.def(py::init([](const XT::Functions::FunctionInterface<d>& func) {
                         return new XT::Functions::FunctionAsGridFunctionWrapper<E, 1, 1, double>(func);
                       }),
                       py::keep_alive<1, 2>(),
                       "scalar_function"_a);
  });
  XT::Common::bindings::try_register(m, [](auto& m_) {
    py::class_<XT::Functions::FunctionAsGridFunctionWrapper<E, d, d, double>,
               XT::Functions::GridFunctionInterface<E, d, d, double>>
        matrix_wrapper(m_, "MatrixFunctionAsGridFunctionWrapper", "MatrixFunctionAsGridFunctionWrapper");
    matrix_wrapper.def(py::init([](const XT::Functions::FunctionInterface<d, d, d>& func) {
                         return new XT::Functions::FunctionAsGridFunctionWrapper<E, d, d, double>(func);
                       }),
                       py::keep_alive<1, 2>(),
                       "matrix_function"_a);
  });

  m.def("function_to_grid_function",
        [](XT::Functions::FunctionInterface<d>& func) {
          return std::make_unique<XT::Functions::FunctionAsGridFunctionWrapper<E, 1, 1, double>>(func);
        },
        py::keep_alive<0, 1>(),
        "scalar_function"_a);
  m.def("function_to_grid_function",
        [](XT::Functions::FunctionInterface<d, d, d>& func) {
          return std::make_unique<XT::Functions::FunctionAsGridFunctionWrapper<E, d, d, double>>(func);
        },
        py::keep_alive<0, 1>(),
        "matrix_function"_a);

  py::class_<DomainDecomposition> domain_decomposition(m, "DomainDecomposition", "DomainDecomposition");
  domain_decomposition.def(py::init([](const std::array<unsigned int, d> num_macro_elements_per_dim,
                                       const size_t num_refinements_per_subdomain,
                                       const FieldVector<double, d> ll,
                                       const FieldVector<double, d> ur) {
                             return new DomainDecomposition(
                                 num_macro_elements_per_dim, num_refinements_per_subdomain, ll, ur);
                           }),
                           "num_macro_elements_per_dim"_a,
                           "num_refinements_per_subdomain"_a,
                           "lower_left"_a,
                           "upper_right"_a);
  domain_decomposition.def_property_readonly("num_subdomains",
                                             [](DomainDecomposition& self) { return self.dd_grid.num_subdomains(); });
  domain_decomposition.def_property_readonly("boundary_subdomains", [](DomainDecomposition& self) {
    std::vector<size_t> boundary_subdomains;
    for (auto&& macro_element : elements(self.dd_grid.macro_grid_view()))
      if (self.dd_grid.boundary(macro_element))
        boundary_subdomains.push_back(self.dd_grid.subdomain(macro_element));
    return boundary_subdomains;
  });
  domain_decomposition.def("neighbors", [](DomainDecomposition& self, const size_t ss) {
    DUNE_THROW_IF(ss >= self.dd_grid.num_subdomains(),
                  XT::Common::Exceptions::index_out_of_range,
                  "ss = " << ss << "\n   self.dd_grid.num_subdomains() = " << self.dd_grid.num_subdomains());
    std::vector<size_t> neighboring_subdomains;
    for (auto&& macro_element : elements(self.dd_grid.macro_grid_view())) {
      if (self.dd_grid.subdomain(macro_element) == ss) {
        for (auto&& macro_intersection : intersections(self.dd_grid.macro_grid_view(), macro_element))
          if (macro_intersection.neighbor())
            neighboring_subdomains.push_back(self.dd_grid.subdomain(macro_intersection.outside()));
        break;
      }
    }
    return neighboring_subdomains;
  });
  domain_decomposition.def(
      "visualize_indicators",
      [](DomainDecomposition& self,
         const std::vector<double>& indicators,
         const std::string& filename,
         const std::string& plot_name) { self.visualize_indicators(indicators, filename, plot_name); },
      py::call_guard<py::gil_scoped_release>(),
      "indicators"_a,
      "filename"_a,
      "plot_name"_a = std::string("indicators"));
  domain_decomposition.def(
      "visualize_local",
      [](DomainDecomposition& self,
         const std::string& filename_prefix,
         const size_t ss,
         const V& vec,
         const std::string& space_type,
         const std::string& name) { self.visualize_local(filename_prefix, ss, vec, space_type, name); },
      py::call_guard<py::gil_scoped_release>(),
      "filename_prefix"_a,
      "ss"_a,
      "subdomain_vector"_a,
      "space_type"_a = default_space_type(),
      "name"_a = "STATE");
  domain_decomposition.def("add_local_visualization",
                           [](DomainDecomposition& self,
                              const size_t ss,
                              const V& vec,
                              const std::string& space_type,
                              const std::string& name) { self.add_local_visualization(ss, vec, space_type, name); },
                           "ss"_a,
                           "subdomain_vector"_a,
                           "space_type"_a = default_space_type(),
                           "name"_a = "STATE");
  domain_decomposition.def(
      "add_local_visualization",
      [](DomainDecomposition& self, const size_t ss, const XT::Functions::FunctionInterface<d>& func) {
        self.add_local_visualization(ss, func);
      },
      "ss"_a,
      "function"_a);
  domain_decomposition.def(
      "write_visualization",
      [](DomainDecomposition& self, const std::string& filename_prefix) { self.write_visualization(filename_prefix); },
      "filename_prefix"_a);

  py::class_<ContinuousLagrangePartitionOfUnity> cg_pou(
      m, "ContinuousLagrangePartitionOfUnity", "ContinuousLagrangePartitionOfUnity");
  cg_pou.def(py::init([](const DomainDecomposition& dd) { return new ContinuousLagrangePartitionOfUnity(dd); }),
             "domain_decomposition"_a,
             py::keep_alive<1, 2>());
  cg_pou.def_property_readonly("size", [](ContinuousLagrangePartitionOfUnity& self) { return self.size(); });
  cg_pou.def("patch", [](ContinuousLagrangePartitionOfUnity& self, const size_t ii) { return self.patch(ii); }, "ii"_a);
  cg_pou.def("local_index",
             [](ContinuousLagrangePartitionOfUnity& self, const size_t ii, const size_t ss) {
               return self.local_index(ii, ss);
             },
             "patch_index"_a,
             "subdomain"_a);
  cg_pou.def("on_subdomain",
             [](ContinuousLagrangePartitionOfUnity& self, const size_t ss, const std::string space_type) {
               return self.on_subdomain(ss, space_type);
             },
             "ss"_a,
             "space_type"_a = default_space_type());

  py::class_<ContinuousFlatTopPartitionOfUnity> flattop_pou(
      m, "ContinuousFlatTopPartitionOfUnity", "ContinuousFlatTopPartitionOfUnity");
  flattop_pou.def(py::init([](const DomainDecomposition& dd, const double& overlap) {
                    return new ContinuousFlatTopPartitionOfUnity(dd, overlap);
                  }),
                  "domain_decomposition"_a,
                  "overlap"_a = 0.5,
                  py::keep_alive<1, 2>());
  flattop_pou.def_property_readonly("size", [](ContinuousFlatTopPartitionOfUnity& self) { return self.size(); });
  flattop_pou.def(
      "patch", [](ContinuousFlatTopPartitionOfUnity& self, const size_t ii) { return self.patch(ii); }, "ii"_a);
  flattop_pou.def("local_index",
                  [](ContinuousFlatTopPartitionOfUnity& self, const size_t ii, const size_t ss) {
                    return self.local_index(ii, ss);
                  },
                  "patch_index"_a,
                  "subdomain"_a);
  flattop_pou.def("on_subdomain",
                  [](ContinuousFlatTopPartitionOfUnity& self, const size_t ss, const std::string space_type) {
                    return self.on_subdomain(ss, space_type);
                  },
                  "ss"_a,
                  "space_type"_a = default_space_type());

  m.def("assemble_local_system_matrix",
        [](XT::Functions::FunctionInterface<d>& diffusion_factor,
           DomainDecomposition& domain_decomposition,
           const size_t ss,
           const std::string space_type) {
          const XT::Functions::ConstantFunction<d, d, d> diffusion_tensor(
              XT::Common::from_string<typename XT::Functions::ConstantFunction<d, d, d>::RangeReturnType>(
                  "[1 0 0; 0 1 0; 0 0 1]"));
          return std::move(assemble_local_system_matrix(diffusion_factor.as_grid_function<E>(),
                                                        diffusion_tensor.as_grid_function<E>(),
                                                        domain_decomposition,
                                                        ss,
                                                        space_type));
        },
        py::call_guard<py::gil_scoped_release>(),
        "diffusion_factor"_a,
        "domain_decomposition"_a,
        "ss"_a,
        "space_type"_a = default_space_type());
  m.def("assemble_local_system_matrix",
        [](XT::Functions::GridFunctionInterface<E>& diffusion_factor,
           XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
           DomainDecomposition& domain_decomposition,
           const size_t ss,
           const std::string space_type) {
          return std::move(
              assemble_local_system_matrix(diffusion_factor, diffusion_tensor, domain_decomposition, ss, space_type));
        },
        py::call_guard<py::gil_scoped_release>(),
        "diffusion_factor"_a,
        "diffusion_tensor"_a,
        "domain_decomposition"_a,
        "ss"_a,
        "space_type"_a = default_space_type());

  m.def("assemble_local_l2_matrix",
        [](DomainDecomposition& domain_decomposition, const size_t ss, const std::string space_type) {
          DUNE_THROW_IF(ss >= domain_decomposition.dd_grid.num_subdomains(),
                        XT::Common::Exceptions::index_out_of_range,
                        "ss = " << ss << "\n   domain_decomposition.dd_grid.num_subdomains() = "
                                << domain_decomposition.dd_grid.num_subdomains());
          std::unique_ptr<M> subdomain_matrix;
          for (auto&& macro_element : elements(domain_decomposition.dd_grid.macro_grid_view())) {
            if (domain_decomposition.dd_grid.subdomain(macro_element) == ss) {
              // this is the subdomain we are interested in, create space
              auto subdomain_grid_view = domain_decomposition.dd_grid.local_grid(macro_element).leaf_view();
              using GV = decltype(subdomain_grid_view);
              using E = typename GV::template Codim<0>::Entity;
              using I = typename GV::Intersection;
              auto subdomain_space = make_subdomain_space(subdomain_grid_view, space_type);
              // create operator
              auto subdomain_operator = make_matrix_operator<M>(*subdomain_space, Stencil::element_and_intersection);
              const LocalElementIntegralBilinearForm<E> element_bilinear_form(LocalElementProductIntegrand<E>{});
              subdomain_operator.append(element_bilinear_form);
              subdomain_operator.assemble();
              subdomain_matrix = std::make_unique<M>(subdomain_operator.matrix());
              break;
            }
          }
          return std::move(subdomain_matrix);
        },
        py::call_guard<py::gil_scoped_release>(),
        "domain_decomposition"_a,
        "ss"_a,
        "space_type"_a = default_space_type());

  m.def("assemble_local_rhs",
        [](XT::Functions::FunctionInterface<d>& force,
           DomainDecomposition& domain_decomposition,
           const size_t ss,
           const std::string space_type) {
          return std::move(assemble_local_rhs(force.as_grid_function<E>(), domain_decomposition, ss, space_type));
        },
        py::call_guard<py::gil_scoped_release>(),
        "force"_a,
        "domain_decomposition"_a,
        "ss"_a,
        "space_type"_a = default_space_type());
  m.def("assemble_local_rhs",
        [](XT::Functions::GridFunctionInterface<E>& force,
           DomainDecomposition& domain_decomposition,
           const size_t ss,
           const std::string space_type) {
          return std::move(assemble_local_rhs(force, domain_decomposition, ss, space_type));
        },
        py::call_guard<py::gil_scoped_release>(),
        "force"_a,
        "domain_decomposition"_a,
        "ss"_a,
        "space_type"_a = default_space_type());

  m.def("assemble_coupling_matrices",
        [](XT::Functions::FunctionInterface<d>& diffusion_factor,
           DomainDecomposition& domain_decomposition,
           const size_t ss,
           const size_t nn,
           const std::string space_type) {
          const XT::Functions::ConstantFunction<d, d, d> diffusion_tensor(
              XT::Common::from_string<typename XT::Functions::ConstantFunction<d, d, d>::RangeReturnType>(
                  "[1 0 0; 0 1 0; 0 0 1]"));
          return assemble_coupling_matrices(diffusion_factor.as_grid_function<E>(),
                                            diffusion_tensor.as_grid_function<E>(),
                                            domain_decomposition,
                                            ss,
                                            nn,
                                            space_type);
        },
        py::call_guard<py::gil_scoped_release>(),
        "diffusion_factor"_a,
        "domain_decomposition"_a,
        "ss"_a,
        "nn"_a,
        "space_type"_a = default_space_type());
  m.def("assemble_coupling_matrices",
        [](XT::Functions::GridFunctionInterface<E>& diffusion_factor,
           XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
           DomainDecomposition& domain_decomposition,
           const size_t ss,
           const size_t nn,
           const std::string space_type) {
          return assemble_coupling_matrices(
              diffusion_factor, diffusion_tensor, domain_decomposition, ss, nn, space_type);
        },
        py::call_guard<py::gil_scoped_release>(),
        "diffusion_factor"_a,
        "diffusion_tensor"_a,
        "domain_decomposition"_a,
        "ss"_a,
        "nn"_a,
        "space_type"_a = default_space_type());

  m.def("assemble_boundary_matrix",
        [](XT::Functions::FunctionInterface<d>& diffusion_factor,
           DomainDecomposition& domain_decomposition,
           const size_t ss,
           const std::string space_type) {
          const XT::Functions::ConstantFunction<d, d, d> diffusion_tensor(
              XT::Common::from_string<typename XT::Functions::ConstantFunction<d, d, d>::RangeReturnType>(
                  "[1 0 0; 0 1 0; 0 0 1]"));
          return std::move(assemble_boundary_matrix(diffusion_factor.as_grid_function<E>(),
                                                    diffusion_tensor.as_grid_function<E>(),
                                                    domain_decomposition,
                                                    ss,
                                                    space_type));
        },
        py::call_guard<py::gil_scoped_release>(),
        "diffusion_factor"_a,
        "domain_decomposition"_a,
        "ss"_a,
        "space_type"_a = default_space_type());
  m.def("assemble_boundary_matrix",
        [](XT::Functions::GridFunctionInterface<E>& diffusion_factor,
           XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
           DomainDecomposition& domain_decomposition,
           const size_t ss,
           const std::string space_type) {
          return std::move(
              assemble_boundary_matrix(diffusion_factor, diffusion_tensor, domain_decomposition, ss, space_type));
        },
        py::call_guard<py::gil_scoped_release>(),
        "diffusion_factor"_a,
        "diffusion_tensor"_a,
        "domain_decomposition"_a,
        "ss"_a,
        "space_type"_a = default_space_type());

  m.def("assemble_local_product_contributions",
        [](DomainDecomposition& domain_decomposition, const size_t ss, const std::string space_type) {
          const XT::Functions::ConstantFunction<d> diffusion_factor(1.);
          const XT::Functions::ConstantFunction<d, d, d> diffusion_tensor(
              XT::Common::from_string<typename XT::Functions::ConstantFunction<d, d, d>::RangeReturnType>(
                  "[1 0 0; 0 1 0; 0 0 1]"));
          return std::move(assemble_local_product_contributions(diffusion_factor.as_grid_function<E>(),
                                                                diffusion_tensor.as_grid_function<E>(),
                                                                domain_decomposition,
                                                                ss,
                                                                space_type));
        },
        py::call_guard<py::gil_scoped_release>(),
        "domain_decomposition"_a,
        "ss"_a,
        "space_type"_a = default_space_type());
  m.def("assemble_local_product_contributions",
        [](XT::Functions::GridFunctionInterface<E>& diffusion_factor,
           XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
           DomainDecomposition& domain_decomposition,
           const size_t ss,
           const std::string space_type) {
          return std::move(assemble_local_product_contributions(
              diffusion_factor, diffusion_tensor, domain_decomposition, ss, space_type));
        },
        py::call_guard<py::gil_scoped_release>(),
        "diffusion_factor"_a,
        "diffusion_tensor"_a,
        "domain_decomposition"_a,
        "ss"_a,
        "space_type"_a = default_space_type());

  m.def("assemble_coupling_product_contributions",
        [](DomainDecomposition& domain_decomposition, const size_t ss, const size_t nn, const std::string space_type) {
          const XT::Functions::ConstantFunction<d> diffusion_factor(1.);
          const XT::Functions::ConstantFunction<d, d, d> diffusion_tensor(
              XT::Common::from_string<typename XT::Functions::ConstantFunction<d, d, d>::RangeReturnType>(
                  "[1 0 0; 0 1 0; 0 0 1]"));
          return std::move(assemble_coupling_product_contributions(diffusion_factor.as_grid_function<E>(),
                                                                   diffusion_tensor.as_grid_function<E>(),
                                                                   domain_decomposition,
                                                                   ss,
                                                                   nn,
                                                                   space_type));
        },
        py::call_guard<py::gil_scoped_release>(),
        "domain_decomposition"_a,
        "ss"_a,
        "nn"_a,
        "space_type"_a = default_space_type());
  m.def("assemble_coupling_product_contributions",
        [](XT::Functions::GridFunctionInterface<E>& diffusion_factor,
           XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
           DomainDecomposition& domain_decomposition,
           const size_t ss,
           const size_t nn,
           const std::string space_type) {
          return std::move(assemble_coupling_product_contributions(
              diffusion_factor, diffusion_tensor, domain_decomposition, ss, nn, space_type));
        },
        py::call_guard<py::gil_scoped_release>(),
        "diffusion_factor"_a,
        "diffusion_tensor"_a,
        "domain_decomposition"_a,
        "ss"_a,
        "nn"_a,
        "space_type"_a = default_space_type());

  m.def("assemble_boundary_product_contributions",
        [](DomainDecomposition& domain_decomposition, const size_t ss, const std::string space_type) {
          const XT::Functions::ConstantFunction<d> diffusion_factor(1.);
          const XT::Functions::ConstantFunction<d, d, d> diffusion_tensor(
              XT::Common::from_string<typename XT::Functions::ConstantFunction<d, d, d>::RangeReturnType>(
                  "[1 0 0; 0 1 0; 0 0 1]"));
          return std::move(assemble_boundary_product_contributions(diffusion_factor.as_grid_function<E>(),
                                                                   diffusion_tensor.as_grid_function<E>(),
                                                                   domain_decomposition,
                                                                   ss,
                                                                   space_type));
        },
        py::call_guard<py::gil_scoped_release>(),
        "domain_decomposition"_a,
        "ss"_a,
        "space_type"_a = default_space_type());
  m.def("assemble_boundary_product_contributions",
        [](XT::Functions::GridFunctionInterface<E>& diffusion_factor,
           XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor,
           DomainDecomposition& domain_decomposition,
           const size_t ss,
           const std::string space_type) {
          return std::move(assemble_boundary_product_contributions(
              diffusion_factor, diffusion_tensor, domain_decomposition, ss, space_type));
        },
        py::call_guard<py::gil_scoped_release>(),
        "diffusion_factor"_a,
        "diffusion_tensor"_a,
        "domain_decomposition"_a,
        "ss"_a,
        "space_type"_a = default_space_type());

  m.def("coupling_DoFs",
        [](DomainDecomposition& domain_decomposition, const size_t ss, const size_t nn, const std::string space_type) {
          DUNE_THROW_IF(ss >= domain_decomposition.dd_grid.num_subdomains(),
                        XT::Common::Exceptions::index_out_of_range,
                        "ss = " << ss << "\n   domain_decomposition.dd_grid.num_subdomains() = "
                                << domain_decomposition.dd_grid.num_subdomains());
          std::set<size_t> edge_DoFs;
          for (auto&& inside_macro_element : elements(domain_decomposition.dd_grid.macro_grid_view())) {
            if (domain_decomposition.dd_grid.subdomain(inside_macro_element) == ss) {
              // this is the subdomain we are interested in, create space
              auto inner_subdomain_grid_view =
                  domain_decomposition.dd_grid.local_grid(inside_macro_element).leaf_view();
              auto inner_subdomain_space = make_subdomain_space(inner_subdomain_grid_view, space_type);
              bool found_correct_macro_intersection = false;
              for (auto&& macro_intersection :
                   intersections(domain_decomposition.dd_grid.macro_grid_view(), inside_macro_element)) {
                if (macro_intersection.neighbor()) {
                  const auto outside_macro_element = macro_intersection.outside();
                  if (domain_decomposition.dd_grid.subdomain(outside_macro_element) == nn) {
                    found_correct_macro_intersection = true;
                    // these are the subdomains we are interested in
                    const auto outer_normal = macro_intersection.centerUnitOuterNormal();
                    using GV = decltype(inner_subdomain_grid_view);
                    using I = GV::Intersection;
                    using SpaceType = GDT::SpaceInterface<GV>;
                    XT::Grid::NormalBasedBoundaryInfo<XT::Grid::extract_intersection_t<GV>> macro_boundary_info;
                    macro_boundary_info.register_new_normal(outer_normal, new XT::Grid::DirichletBoundary());
                    DirichletConstraints<I, SpaceType> constraints(macro_boundary_info, *inner_subdomain_space);
                    const auto& coupling = domain_decomposition.dd_grid.coupling(
                        inside_macro_element, -1, outside_macro_element, -1, true);
                    DynamicVector<size_t> global_indices_in(inner_subdomain_space->mapper().max_local_size());
                    auto inside_basis = inner_subdomain_space->basis().localize();
                    const auto coupling_intersection_it_end = coupling.template iend<0>();
                    for (auto coupling_intersection_it = coupling.template ibegin<0>();
                         coupling_intersection_it != coupling_intersection_it_end;
                         ++coupling_intersection_it) {
                      const auto& coupling_intersection = *coupling_intersection_it;
                      const auto inside_element = coupling_intersection.inside();
                      constraints.apply_local(inside_element);
                    }
                    edge_DoFs = std::move(constraints.dirichlet_DoFs());
                    break;
                  }
                }
              }

              DUNE_THROW_IF(!found_correct_macro_intersection,
                            XT::Common::Exceptions::index_out_of_range,
                            "ss = " << ss << "\n   nn = " << nn);
            }
          }
          return edge_DoFs;
        },
        py::call_guard<py::gil_scoped_release>(),
        "domain_decomposition"_a,
        "ss"_a,
        "nn"_a,
        "space_type"_a = default_space_type());

  m.def("boundary_DoFs",
        [](DomainDecomposition& domain_decomposition, const size_t ss, const std::string space_type) {
          DUNE_THROW_IF(ss >= domain_decomposition.dd_grid.num_subdomains(),
                        XT::Common::Exceptions::index_out_of_range,
                        "ss = " << ss << "\n   domain_decomposition.dd_grid.num_subdomains() = "
                                << domain_decomposition.dd_grid.num_subdomains());
          using MGV = typename DomainDecomposition::DdGridType::MacroGridViewType;
          const XT::Grid::AllDirichletBoundaryInfo<XT::Grid::extract_intersection_t<MGV>> macro_boundary_info;
          std::set<size_t> boundary_DoFs;
          for (auto&& macro_element : elements(domain_decomposition.dd_grid.macro_grid_view())) {
            if (domain_decomposition.dd_grid.subdomain(macro_element) == ss) {
              // this is the subdomain we are interested in, create space
              auto subdomain_grid_view = domain_decomposition.dd_grid.local_grid(macro_element).leaf_view();
              using GV = decltype(subdomain_grid_view);
              using E = typename GV::template Codim<0>::Entity;
              using I = typename GV::Intersection;
              auto subdomain_space = make_subdomain_space(subdomain_grid_view, space_type);
              const MacroGridBasedBoundaryInfo<MGV, GV> subdomain_boundary_info(
                  domain_decomposition.dd_grid.macro_grid_view(), macro_element, macro_boundary_info);
              DirichletConstraints<I, std::decay_t<decltype(*subdomain_space)>> constraints(subdomain_boundary_info,
                                                                                            *subdomain_space);
              auto walker = XT::Grid::make_walker(subdomain_grid_view);
              walker.append(constraints, XT::Grid::ApplyOn::BoundaryElements<GV>());
              walker.walk(/*parallel=*/true);
              boundary_DoFs = std::move(constraints.dirichlet_DoFs());
              break;
            }
          }
          return boundary_DoFs;
        },
        py::call_guard<py::gil_scoped_release>(),
        "domain_decomposition"_a,
        "ss"_a,
        "space_type"_a = default_space_type());

  m.def("to_csr",
        [](M& matrix, const double& prune) {
          const auto& mat = matrix.backend();
          // The istl matrix reports the wrong numbers of nonzerose, so we start with std::vectors ...
          std::vector<double> data;
          std::vector<size_t> rows;
          std::vector<size_t> cols;
          if (prune > 0) {
            for (size_t ii = 0; ii < mat.N(); ++ii) {
              if (mat.getrowsize(ii) > 0) {
                const auto& row = mat[ii];
                const auto it_end = row.end();
                for (auto it = row.begin(); it != it_end; ++it) {
                  const auto val = it->operator[](0)[0];
                  if (XT::Common::FloatCmp::ne<XT::Common::FloatCmp::Style::absolute>(val, decltype(val)(0), prune)) {
                    data.push_back(val);
                    rows.push_back(ii);
                    cols.push_back(it.index());
                  }
                }
              }
            }
          } else {
            for (size_t ii = 0; ii < mat.N(); ++ii) {
              if (mat.getrowsize(ii) > 0) {
                const auto& row = mat[ii];
                const auto it_end = row.end();
                for (auto it = row.begin(); it != it_end; ++it) {
                  data.push_back(it->operator[](0)[0]);
                  rows.push_back(ii);
                  cols.push_back(it.index());
                }
              }
            }
          }
          // ... and convert afterwards.
          return std::make_tuple(std::make_unique<V>(data),
                                 std::make_unique<XT::LA::CommonDenseVector<size_t>>(rows),
                                 std::make_unique<XT::LA::CommonDenseVector<size_t>>(cols));
        },
        py::call_guard<py::gil_scoped_release>(),
        "matrix"_a,
        "prune"_a = 1e-15);

  m.def("estimate_inverse_inequality_constant",
        [](DomainDecomposition& domain_decomposition, const size_t ss, const std::string space_type) {
          DUNE_THROW_IF(ss >= domain_decomposition.dd_grid.num_subdomains(),
                        XT::Common::Exceptions::index_out_of_range,
                        "ss = " << ss << "\n   domain_decomposition.dd_grid.num_subdomains() = "
                                << domain_decomposition.dd_grid.num_subdomains());
          for (auto&& macro_element : elements(domain_decomposition.dd_grid.macro_grid_view())) {
            if (domain_decomposition.dd_grid.subdomain(macro_element) == ss) {
              // this is the subdomain we are interested in, create space
              auto subdomain_grid_view = domain_decomposition.dd_grid.local_grid(macro_element).leaf_view();
              auto subdomain_space = make_subdomain_space(subdomain_grid_view, space_type);
              return estimate_inverse_inequality_constant(*subdomain_space);
            }
          }
          DUNE_THROW(InvalidStateException, "This should not happen, ss = " << ss);
          return 0.;
        },
        py::call_guard<py::gil_scoped_release>(),
        "domain_decomposition"_a,
        "ss"_a,
        "space_type"_a = default_space_type());

  m.def("estimate_combined_inverse_trace_inequality_constant",
        [](DomainDecomposition& domain_decomposition, const size_t ss, const std::string space_type) {
          DUNE_THROW_IF(ss >= domain_decomposition.dd_grid.num_subdomains(),
                        XT::Common::Exceptions::index_out_of_range,
                        "ss = " << ss << "\n   domain_decomposition.dd_grid.num_subdomains() = "
                                << domain_decomposition.dd_grid.num_subdomains());
          for (auto&& macro_element : elements(domain_decomposition.dd_grid.macro_grid_view())) {
            if (domain_decomposition.dd_grid.subdomain(macro_element) == ss) {
              // this is the subdomain we are interested in, create space
              auto subdomain_grid_view = domain_decomposition.dd_grid.local_grid(macro_element).leaf_view();
              auto subdomain_space = make_subdomain_space(subdomain_grid_view, space_type);
              return estimate_combined_inverse_trace_inequality_constant(*subdomain_space);
            }
          }
          DUNE_THROW(InvalidStateException, "This should not happen, ss = " << ss);
          return 0.;
        },
        py::call_guard<py::gil_scoped_release>(),
        "domain_decomposition"_a,
        "ss"_a,
        "space_type"_a = default_space_type());

  m.def("estimate_element_to_intersection_equivalence_constant",
        [](DomainDecomposition& domain_decomposition, const size_t ss) {
          DUNE_THROW_IF(ss >= domain_decomposition.dd_grid.num_subdomains(),
                        XT::Common::Exceptions::index_out_of_range,
                        "ss = " << ss << "\n   domain_decomposition.dd_grid.num_subdomains() = "
                                << domain_decomposition.dd_grid.num_subdomains());
          for (auto&& macro_element : elements(domain_decomposition.dd_grid.macro_grid_view())) {
            if (domain_decomposition.dd_grid.subdomain(macro_element) == ss) {
              // this is the subdomain we are interested in, create space
              auto subdomain_grid_view = domain_decomposition.dd_grid.local_grid(macro_element).leaf_view();
              return estimate_element_to_intersection_equivalence_constant(subdomain_grid_view);
            }
          }
          DUNE_THROW(InvalidStateException, "This should not happen, ss = " << ss);
          return 0.;
        },
        py::call_guard<py::gil_scoped_release>(),
        "domain_decomposition"_a,
        "ss"_a);
} // PYBIND11_MODULE(usercode, ...)
