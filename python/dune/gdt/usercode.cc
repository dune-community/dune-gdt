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

#include <dune/gdt/functionals/vector-based.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/integrands/conversion.hh>
#include <dune/gdt/local/integrands/elliptic-ipdg.hh>
#include <dune/gdt/local/integrands/elliptic.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
#include <dune/gdt/spaces/l2/discontinuous-lagrange.hh>

using namespace Dune;
using namespace Dune::GDT;

using G = YASP_2D_EQUIDISTANT_OFFSET;
static const constexpr size_t d = G::dimension;
using M = XT::LA::IstlRowMajorSparseMatrix<double>;
using V = XT::LA::IstlDenseVector<double>;


class DomainDecomposition
{
public:
  using DdGridType = XT::Grid::DD::Glued<G, G, XT::Grid::Layers::leaf>;

  DomainDecomposition(const unsigned int num_macro_elements_per_dim, const size_t num_refinements_per_subdomain)
    : macro_grid(XT::Grid::make_cube_grid<G>(0., 1., num_macro_elements_per_dim))
    , dd_grid(macro_grid, num_refinements_per_subdomain, false, true)
  {
  }

private:
  XT::Grid::GridProvider<G> macro_grid;

public:
  DdGridType dd_grid;
}; // class DomainDecomposition


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
  {
  }

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


template <class GV>
std::unique_ptr<GDT::SpaceInterface<GV>> make_subdomain_space(GV subdomain_grid_view, const std::string& space_type)
{
  if (space_type == "continuous_lagrange")
    return std::make_unique<GDT::ContinuousLagrangeSpace<GV, /*order=*/1>>(subdomain_grid_view);
  else if (space_type == "discontinuous_lagrange")
    return std::make_unique<GDT::DiscontinuousLagrangeSpace<GV, /*order=*/1>>(subdomain_grid_view);
  else
    DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
               "space_type = " << space_type << "\n   has to be 'continuous_lagrange' or 'discontinuous_lagrange'!");
}


PYBIND11_PLUGIN(usercode)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  py::module m("usercode", "dune-gdt");

  Dune::XT::Common::bindings::addbind_exceptions(m);

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");

  py::class_<DomainDecomposition> domain_decomposition(m, "DomainDecomposition", "DomainDecomposition");
  domain_decomposition.def(
      py::init([](const unsigned int num_macro_elements_per_dim, const size_t num_refinements_per_subdomain) {
        return new DomainDecomposition(num_macro_elements_per_dim, num_refinements_per_subdomain);
      }),
      "num_macro_elements_per_dim"_a,
      "num_refinements_per_subdomain"_a);
  domain_decomposition.def_property_readonly("num_sudomains",
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

  m.def("assemble_local_system_matrix",
        [](XT::Functions::FunctionInterface<d>& diffusion_factor,
           DomainDecomposition& domain_decomposition,
           const size_t ss,
           const std::string space_type) {
          DUNE_THROW_IF(ss >= domain_decomposition.dd_grid.num_subdomains(),
                        XT::Common::Exceptions::index_out_of_range,
                        "ss = " << ss << "\n   domain_decomposition.dd_grid.num_subdomains() = "
                                << domain_decomposition.dd_grid.num_subdomains());
          const XT::Functions::ConstantFunction<d, d, d> diffusion_tensor(
              XT::Common::from_string<typename XT::Functions::ConstantFunction<d, d, d>::RangeReturnType>(
                  "[1 0 0; 0 1 0; 0 0 1]"));
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
              const LocalElementIntegralBilinearForm<E> element_bilinear_form(LocalEllipticIntegrand<E>(
                  diffusion_factor.as_grid_function<E>(), diffusion_tensor.as_grid_function<E>()));
              subdomain_operator.append(element_bilinear_form);
              if (!subdomain_space->continuous(0)) {
                const LocalIntersectionIntegralBilinearForm<I> coupling_bilinear_form(
                    LocalEllipticIpdgIntegrands::Inner<I>(diffusion_factor.as_grid_function<E>(),
                                                          diffusion_tensor.as_grid_function<E>()));
                subdomain_operator.append(coupling_bilinear_form, {}, XT::Grid::ApplyOn::InnerIntersectionsOnce<GV>());
              }
              subdomain_operator.assemble();
              subdomain_matrix = std::make_unique<M>(subdomain_operator.matrix());
              break;
            }
          }
          return std::move(subdomain_matrix);
        },
        "diffusion_factor"_a,
        "domain_decomposition"_a,
        "ss"_a,
        "space_type"_a = "discontinuous_lagrange");

  m.def("assemble_local_rhs",
        [](XT::Functions::FunctionInterface<d>& force,
           DomainDecomposition& domain_decomposition,
           const size_t ss,
           const std::string space_type) {
          DUNE_THROW_IF(ss >= domain_decomposition.dd_grid.num_subdomains(),
                        XT::Common::Exceptions::index_out_of_range,
                        "ss = " << ss << "\n   domain_decomposition.dd_grid.num_subdomains() = "
                                << domain_decomposition.dd_grid.num_subdomains());
          std::unique_ptr<V> subdomain_vector;
          for (auto&& macro_element : elements(domain_decomposition.dd_grid.macro_grid_view())) {
            if (domain_decomposition.dd_grid.subdomain(macro_element) == ss) {
              // this is the subdomain we are interested in, create space
              auto subdomain_grid_view = domain_decomposition.dd_grid.local_grid(macro_element).leaf_view();
              using GV = decltype(subdomain_grid_view);
              using E = typename GV::template Codim<0>::Entity;
              auto subdomain_space = make_subdomain_space(subdomain_grid_view, space_type);
              // create functional
              auto subdomain_functional = make_vector_functional<V>(*subdomain_space);
              const LocalElementIntegralFunctional<E> element_functional(local_binary_to_unary_element_integrand(
                  force.as_grid_function<E>(), LocalElementProductIntegrand<E>()));
              subdomain_functional.append(element_functional);
              subdomain_functional.assemble();
              subdomain_vector = std::make_unique<V>(subdomain_functional.vector());
              break;
            }
          }
          return std::move(subdomain_vector);
        },
        "force"_a,
        "domain_decomposition"_a,
        "ss"_a,
        "space_type"_a = "discontinuous_lagrange");

  m.def("assemble_coupling_matrices",
        [](XT::Functions::FunctionInterface<d>& diffusion_factor,
           DomainDecomposition& domain_decomposition,
           const size_t ss,
           const size_t nn,
           const std::string space_type) {
          DUNE_THROW_IF(ss >= domain_decomposition.dd_grid.num_subdomains(),
                        XT::Common::Exceptions::index_out_of_range,
                        "ss = " << ss << "\n   domain_decomposition.dd_grid.num_subdomains() = "
                                << domain_decomposition.dd_grid.num_subdomains());
          const XT::Functions::ConstantFunction<d, d, d> diffusion_tensor(
              XT::Common::from_string<typename XT::Functions::ConstantFunction<d, d, d>::RangeReturnType>(
                  "[1 0 0; 0 1 0; 0 0 1]"));
          std::unique_ptr<M> coupling_matrix_in_in;
          std::unique_ptr<M> coupling_matrix_in_out;
          std::unique_ptr<M> coupling_matrix_out_in;
          std::unique_ptr<M> coupling_matrix_out_out;
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
                    auto outer_subdomain_grid_view =
                        domain_decomposition.dd_grid.local_grid(outside_macro_element).leaf_view();
                    auto outer_subdomain_space = make_subdomain_space(outer_subdomain_grid_view, space_type);
                    // walk the coupling once to create the sparsity patterns
                    XT::LA::SparsityPatternDefault pattern_in_in(inner_subdomain_space->mapper().size());
                    XT::LA::SparsityPatternDefault pattern_in_out(inner_subdomain_space->mapper().size());
                    XT::LA::SparsityPatternDefault pattern_out_in(outer_subdomain_space->mapper().size());
                    XT::LA::SparsityPatternDefault pattern_out_out(outer_subdomain_space->mapper().size());
                    const auto& coupling = domain_decomposition.dd_grid.coupling(
                        inside_macro_element, -1, outside_macro_element, -1, true);
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
                        outer_subdomain_space->mapper().size(), outer_subdomain_space->mapper().size(), pattern_in_in);
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
                        LocalEllipticIpdgIntegrands::Inner<I>(diffusion_factor.as_grid_function<E>(),
                                                              diffusion_tensor.as_grid_function<E>()));
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
        },
        "diffusion_factor"_a,
        "domain_decomposition"_a,
        "ss"_a,
        "nn"_a,
        "space_type"_a = "discontinuous_lagrange");

  m.def("assemble_boundary_matrices",
        [](XT::Functions::FunctionInterface<d>& diffusion_factor,
           DomainDecomposition& domain_decomposition,
           const size_t ss,
           const std::string space_type) {
          DUNE_THROW_IF(ss >= domain_decomposition.dd_grid.num_subdomains(),
                        XT::Common::Exceptions::index_out_of_range,
                        "ss = " << ss << "\n   domain_decomposition.dd_grid.num_subdomains() = "
                                << domain_decomposition.dd_grid.num_subdomains());
          using MGV = typename DomainDecomposition::DdGridType::MacroGridViewType;
          const XT::Grid::AllDirichletBoundaryInfo<XT::Grid::extract_intersection_t<MGV>> macro_boundary_info;
          const XT::Functions::ConstantFunction<d, d, d> diffusion_tensor(
              XT::Common::from_string<typename XT::Functions::ConstantFunction<d, d, d>::RangeReturnType>(
                  "[1 0 0; 0 1 0; 0 0 1]"));
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
                  LocalEllipticIpdgIntegrands::DirichletBoundaryLhs<I>(diffusion_factor.as_grid_function<E>(),
                                                                       diffusion_tensor.as_grid_function<E>()));
              subdomain_operator.append(dirichlet_bilinear_form,
                                        {},
                                        XT::Grid::ApplyOn::CustomBoundaryIntersections<GV>(
                                            subdomain_boundary_info, new XT::Grid::DirichletBoundary()));
              subdomain_operator.assemble();
              subdomain_matrix = std::make_unique<M>(subdomain_operator.matrix());
              break;
            }
          }
          return std::move(subdomain_matrix);
        },
        "diffusion_factor"_a,
        "domain_decomposition"_a,
        "ss"_a,
        "space_type"_a = "discontinuous_lagrange");

  m.def("assemble_local_product_contributions",
        [](DomainDecomposition& domain_decomposition, const size_t ss, const std::string space_type) {
          DUNE_THROW_IF(ss >= domain_decomposition.dd_grid.num_subdomains(),
                        XT::Common::Exceptions::index_out_of_range,
                        "ss = " << ss << "\n   domain_decomposition.dd_grid.num_subdomains() = "
                                << domain_decomposition.dd_grid.num_subdomains());
          const XT::Functions::ConstantFunction<d> diffusion_factor(1.);
          const XT::Functions::ConstantFunction<d, d, d> diffusion_tensor(
              XT::Common::from_string<typename XT::Functions::ConstantFunction<d, d, d>::RangeReturnType>(
                  "[1 0 0; 0 1 0; 0 0 1]"));
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
              const LocalElementIntegralBilinearForm<E> element_bilinear_form(LocalEllipticIntegrand<E>(
                  diffusion_factor.as_grid_function<E>(), diffusion_tensor.as_grid_function<E>()));
              subdomain_operator.append(element_bilinear_form);
              if (!subdomain_space->continuous(0)) {
                const LocalIntersectionIntegralBilinearForm<I> coupling_bilinear_form(
                    LocalEllipticIpdgIntegrands::InnerOnlyPenalty<I>(diffusion_factor.as_grid_function<E>(),
                                                                     diffusion_tensor.as_grid_function<E>()));
                subdomain_operator.append(coupling_bilinear_form, {}, XT::Grid::ApplyOn::InnerIntersectionsOnce<GV>());
              }
              subdomain_operator.assemble();
              subdomain_matrix = std::make_unique<M>(subdomain_operator.matrix());
              break;
            }
          }
          return std::move(subdomain_matrix);
        },
        "domain_decomposition"_a,
        "ss"_a,
        "space_type"_a = "discontinuous_lagrange");

  m.def("assemble_coupling_product_contributions",
        [](DomainDecomposition& domain_decomposition, const size_t ss, const size_t nn, const std::string space_type) {
          DUNE_THROW_IF(ss >= domain_decomposition.dd_grid.num_subdomains(),
                        XT::Common::Exceptions::index_out_of_range,
                        "ss = " << ss << "\n   domain_decomposition.dd_grid.num_subdomains() = "
                                << domain_decomposition.dd_grid.num_subdomains());
          const XT::Functions::ConstantFunction<d> diffusion_factor(1.);
          const XT::Functions::ConstantFunction<d, d, d> diffusion_tensor(
              XT::Common::from_string<typename XT::Functions::ConstantFunction<d, d, d>::RangeReturnType>(
                  "[1 0 0; 0 1 0; 0 0 1]"));
          std::unique_ptr<M> coupling_matrix_in_in;
          std::unique_ptr<M> coupling_matrix_in_out;
          std::unique_ptr<M> coupling_matrix_out_in;
          std::unique_ptr<M> coupling_matrix_out_out;
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
                    auto outer_subdomain_grid_view =
                        domain_decomposition.dd_grid.local_grid(outside_macro_element).leaf_view();
                    auto outer_subdomain_space = make_subdomain_space(outer_subdomain_grid_view, space_type);
                    // walk the coupling once to create the sparsity patterns
                    XT::LA::SparsityPatternDefault pattern_in_in(inner_subdomain_space->mapper().size());
                    XT::LA::SparsityPatternDefault pattern_in_out(inner_subdomain_space->mapper().size());
                    XT::LA::SparsityPatternDefault pattern_out_in(outer_subdomain_space->mapper().size());
                    XT::LA::SparsityPatternDefault pattern_out_out(outer_subdomain_space->mapper().size());
                    const auto& coupling = domain_decomposition.dd_grid.coupling(
                        inside_macro_element, -1, outside_macro_element, -1, true);
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
                        outer_subdomain_space->mapper().size(), outer_subdomain_space->mapper().size(), pattern_in_in);
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
                        LocalEllipticIpdgIntegrands::InnerOnlyPenalty<I>(diffusion_factor.as_grid_function<E>(),
                                                                         diffusion_tensor.as_grid_function<E>()));
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
        },
        "domain_decomposition"_a,
        "ss"_a,
        "nn"_a,
        "space_type"_a = "discontinuous_lagrange");

  m.def("assemble_boundary_product_contributions",
        [](DomainDecomposition& domain_decomposition, const size_t ss, const std::string space_type) {
          DUNE_THROW_IF(ss >= domain_decomposition.dd_grid.num_subdomains(),
                        XT::Common::Exceptions::index_out_of_range,
                        "ss = " << ss << "\n   domain_decomposition.dd_grid.num_subdomains() = "
                                << domain_decomposition.dd_grid.num_subdomains());
          using MGV = typename DomainDecomposition::DdGridType::MacroGridViewType;
          const XT::Grid::AllDirichletBoundaryInfo<XT::Grid::extract_intersection_t<MGV>> macro_boundary_info;
          const XT::Functions::ConstantFunction<d> diffusion_factor(1.);
          const XT::Functions::ConstantFunction<d, d, d> diffusion_tensor(
              XT::Common::from_string<typename XT::Functions::ConstantFunction<d, d, d>::RangeReturnType>(
                  "[1 0 0; 0 1 0; 0 0 1]"));
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
                  LocalEllipticIpdgIntegrands::DirichletBoundaryLhsOnlyPenalty<I>(
                      diffusion_factor.as_grid_function<E>(), diffusion_tensor.as_grid_function<E>()));
              subdomain_operator.append(dirichlet_bilinear_form,
                                        {},
                                        XT::Grid::ApplyOn::CustomBoundaryIntersections<GV>(
                                            subdomain_boundary_info, new XT::Grid::DirichletBoundary()));
              subdomain_operator.assemble();
              subdomain_matrix = std::make_unique<M>(subdomain_operator.matrix());
              break;
            }
          }
          return std::move(subdomain_matrix);
        },
        "domain_decomposition"_a,
        "ss"_a,
        "space_type"_a = "discontinuous_lagrange");

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
          return std::make_tuple(std::make_unique<XT::LA::CommonDenseVector<double>>(data),
                                 std::make_unique<XT::LA::CommonDenseVector<size_t>>(rows),
                                 std::make_unique<XT::LA::CommonDenseVector<size_t>>(cols));
        },
        "matrix"_a,
        "prune"_a = 1e-15);

  Dune::XT::Common::bindings::add_initialization(m, "dune.gdt");
  return m.ptr();
} // PYBIND11_PLUGIN(...)
