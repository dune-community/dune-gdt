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

#include "usercode.hh"


PYBIND11_MODULE(usercode, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");
  py::module::import("dune.gdt.discretefunction");

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
  domain_decomposition.def(
      "num_elements",
      [](DomainDecomposition& self, const size_t ss) {
        DUNE_THROW_IF(ss >= self.dd_grid.num_subdomains(),
                      XT::Common::Exceptions::index_out_of_range,
                      "ss = " << ss << "\n   self.dd_grid.num_subdomains() = " << self.dd_grid.num_subdomains());
        size_t result = 0;
        for (auto&& macro_element : elements(self.dd_grid.macro_grid_view()))
          if (self.dd_grid.subdomain(macro_element) == ss) {
            result = make_finite_volume_space(self.dd_grid.local_grid(macro_element).leaf_view()).mapper().size();
            break;
          }
        return result;
      },
      py::call_guard<py::gil_scoped_release>(),
      "subdomain"_a);
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
  domain_decomposition.def("max_subdomain_diameter",
                           [](DomainDecomposition& self) {
                             double H = 0.;
                             for (auto&& macro_element : elements(self.dd_grid.macro_grid_view()))
                               H = std::max(H, XT::Grid::diameter(macro_element));
                             return H;
                           },
                           py::call_guard<py::gil_scoped_release>());
  domain_decomposition.def(
      "max_element_diameter",
      [](DomainDecomposition& self, const size_t ss) {
        DUNE_THROW_IF(ss >= self.dd_grid.num_subdomains(),
                      XT::Common::Exceptions::index_out_of_range,
                      "ss = " << ss << "\n   self.dd_grid.num_subdomains() = " << self.dd_grid.num_subdomains());
        double h = 0.;
        for (auto&& macro_element : elements(self.dd_grid.macro_grid_view()))
          if (self.dd_grid.subdomain(macro_element) == ss) {
            for (auto&& micro_element : elements(self.dd_grid.local_grid(macro_element).leaf_view()))
              h = std::max(h, XT::Grid::diameter(micro_element));
            break;
          }
        return h;
      },
      py::call_guard<py::gil_scoped_release>(),
      "subdomain"_a);
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
  cg_pou.def(py::init([](DomainDecomposition& dd) { return new ContinuousLagrangePartitionOfUnity(dd); }),
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
  cg_pou.def("visualize",
             [](ContinuousLagrangePartitionOfUnity& self, const std::string& filename, const std::string space_type) {
               self.visualize(filename, "CgPoU", space_type);
             },
             py::call_guard<py::gil_scoped_release>(),
             "filename"_a,
             "space_type"_a = default_space_type());

  py::class_<ContinuousFlatTopPartitionOfUnity> flattop_pou(
      m, "ContinuousFlatTopPartitionOfUnity", "ContinuousFlatTopPartitionOfUnity");
  flattop_pou.def(py::init([](DomainDecomposition& dd, const double& overlap) {
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
  flattop_pou.def("visualize",
                  [](ContinuousFlatTopPartitionOfUnity& self,
                     const std::string& filename,
                     const std::string space_type) { self.visualize(filename, "FlatTopPoU", space_type); },
                  py::call_guard<py::gil_scoped_release>(),
                  "filename"_a,
                  "space_type"_a = default_space_type());

  m.def("assemble_local_system_matrix",
        [](XT::Functions::GridFunctionInterface<E>& diffusion,
           const double& penalty,
           XT::Functions::GridFunctionInterface<E>& weight,
           DomainDecomposition& domain_decomposition,
           const size_t ss,
           const std::string space_type) {
          return std::move(
              assemble_local_system_matrix(diffusion, penalty, weight, domain_decomposition, ss, space_type));
        },
        py::call_guard<py::gil_scoped_release>(),
        "diffusion"_a,
        "penalty"_a,
        "weight"_a,
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
        [](XT::Functions::GridFunctionInterface<E>& diffusion,
           const double& penalty,
           XT::Functions::GridFunctionInterface<E>& weight,
           DomainDecomposition& domain_decomposition,
           const size_t ss,
           const size_t nn,
           const std::string space_type) {
          return assemble_coupling_matrices(diffusion, penalty, weight, domain_decomposition, ss, nn, space_type);
        },
        py::call_guard<py::gil_scoped_release>(),
        "diffusion"_a,
        "penalty"_a,
        "weight"_a,
        "domain_decomposition"_a,
        "ss"_a,
        "nn"_a,
        "space_type"_a = default_space_type());
  m.def("assemble_boundary_matrix",
        [](XT::Functions::GridFunctionInterface<E>& diffusion,
           const double& penalty,
           XT::Functions::GridFunctionInterface<E>& weight,
           DomainDecomposition& domain_decomposition,
           const size_t ss,
           const std::string space_type) {
          return std::move(assemble_boundary_matrix(diffusion, penalty, weight, domain_decomposition, ss, space_type));
        },
        py::call_guard<py::gil_scoped_release>(),
        "diffusion"_a,
        "penalty"_a,
        "weight"_a,
        "domain_decomposition"_a,
        "ss"_a,
        "space_type"_a = default_space_type());

  m.def("assemble_local_product_contributions",
        [](const double& penalty,
           XT::Functions::GridFunctionInterface<E>& weight,
           DomainDecomposition& domain_decomposition,
           const size_t ss,
           const std::string space_type) {
          return std::move(assemble_local_product_contributions(penalty, weight, domain_decomposition, ss, space_type));
        },
        py::call_guard<py::gil_scoped_release>(),
        "penalty"_a,
        "weight"_a,
        "domain_decomposition"_a,
        "ss"_a,
        "space_type"_a = default_space_type());

  m.def("assemble_coupling_product_contributions",
        [](const double& penalty,
           XT::Functions::GridFunctionInterface<E>& weight,
           DomainDecomposition& domain_decomposition,
           const size_t ss,
           const size_t nn,
           const std::string space_type) {
          return std::move(
              assemble_coupling_product_contributions(penalty, weight, domain_decomposition, ss, nn, space_type));
        },
        py::call_guard<py::gil_scoped_release>(),
        "penalty"_a,
        "weight"_a,
        "domain_decomposition"_a,
        "ss"_a,
        "nn"_a,
        "space_type"_a = default_space_type());

  m.def("assemble_boundary_product_contributions",
        [](const double& penalty,
           XT::Functions::GridFunctionInterface<E>& weight,
           DomainDecomposition& domain_decomposition,
           const size_t ss,
           const std::string space_type) {
          return std::move(
              assemble_boundary_product_contributions(penalty, weight, domain_decomposition, ss, space_type));
        },
        py::call_guard<py::gil_scoped_release>(),
        "penalty"_a,
        "weight"_a,
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
