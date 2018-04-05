// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef PYTHON_DUNE_GDT_PLAYGROUND_SPACES_BLOCK_BINDINGS_HH
#define PYTHON_DUNE_GDT_PLAYGROUND_SPACES_BLOCK_BINDINGS_HH

#include <dune/pybindxi/pybind11.h>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/la/container.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/projections.hh>
#include <python/dune/gdt/spaces/interface.hh>
#include <dune/gdt/type_traits.hh>

#include <dune/gdt/playground/spaces/block.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class SP>
class BlockMapper
{
  typedef typename SP::type S;
  static_assert(is_space<S>::value, "");
  using G = XT::Grid::extract_grid_t<typename S::GridLayerType>;

public:
  typedef GDT::BlockMapper<S> type;
  typedef pybind11::class_<type> bound_type;

  template <XT::LA::Backends la>
  static void addbind_matrix(bound_type& c)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    typedef typename XT::LA::Container<double, la>::MatrixType M;

    c.def("copy_local_to_global",
          [](const type& self,
             const M& local_matrix,
             const XT::LA::SparsityPatternDefault& local_pattern,
             const ssize_t block,
             M& global_matrix) {
            auto bb = XT::Common::numeric_cast<size_t>(block);
            for (size_t local_ii = 0; local_ii < local_pattern.size(); ++local_ii) {
              const size_t global_ii = self.mapToGlobal(bb, local_ii);
              for (const size_t& local_jj : local_pattern.inner(local_ii)) {
                const size_t global_jj = self.mapToGlobal(bb, local_jj);
                global_matrix.add_to_entry(global_ii, global_jj, local_matrix.get_entry(local_ii, local_jj));
              }
            }
          },
          "local_matrix"_a,
          "local_sparsity_pattern"_a,
          "block"_a,
          "global_matrix"_a);
    c.def("copy_local_to_global",
          [](const type& self,
             const M& local_matrix,
             const XT::LA::SparsityPatternDefault& local_pattern,
             const ssize_t test_block,
             const ssize_t ansatz_block,
             M& global_matrix) {
            auto tt = XT::Common::numeric_cast<size_t>(test_block);
            auto aa = XT::Common::numeric_cast<size_t>(ansatz_block);
            for (size_t local_ii = 0; local_ii < local_pattern.size(); ++local_ii) {
              const size_t global_ii = self.mapToGlobal(tt, local_ii);
              for (const size_t& local_jj : local_pattern.inner(local_ii)) {
                const size_t global_jj = self.mapToGlobal(aa, local_jj);
                global_matrix.add_to_entry(global_ii, global_jj, local_matrix.get_entry(local_ii, local_jj));
              }
            }
          },
          "local_matrix"_a,
          "local_sparsity_pattern"_a,
          "test_block"_a,
          "ansatz_block"_a,
          "global_matrix"_a);
  } // ... addbind_matrix(...)

  template <XT::LA::Backends la>
  static void addbind_vector(bound_type& c)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    typedef typename XT::LA::Container<double, la>::VectorType V;

    c.def("copy_local_to_global",
          [](const type& self, const V& local_vector, const ssize_t block, V& global_vector) {
            auto bb = XT::Common::numeric_cast<size_t>(block);
            for (size_t local_ii = 0; local_ii < local_vector.size(); ++local_ii) {
              const size_t global_ii = self.mapToGlobal(bb, local_ii);
              global_vector.add_to_entry(global_ii, local_vector.get_entry(local_ii));
            }
          },
          "local_vector"_a,
          "block"_a,
          "global_vector"_a);
  } // ... addbind_vector(...)

  static bound_type bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = XT::Common::to_camel_case("block_" + space_name<SP>::value() + "_mapper");

    bound_type c(m, ClassName.c_str());

    c.def_property_readonly("size", [](const type& self) { return self.size(); });

    // sparsity patterns
    c.def("copy_local_to_global",
          [](const type& self,
             const XT::LA::SparsityPatternDefault& local_pattern,
             const ssize_t block,
             XT::LA::SparsityPatternDefault& global_pattern) {
            auto bb = XT::Common::numeric_cast<size_t>(block);
            for (size_t local_ii = 0; local_ii < local_pattern.size(); ++local_ii) {
              const size_t global_ii = self.mapToGlobal(bb, local_ii);
              const auto& local_rows = local_pattern.inner(local_ii);
              for (const auto& local_jj : local_rows) {
                const size_t global_jj = self.mapToGlobal(bb, local_jj);
                global_pattern.insert(global_ii, global_jj);
              }
            }
          },
          "local_sparsity_pattern"_a,
          "block"_a,
          "global_sparsity_pattern"_a);
    c.def("copy_local_to_global",
          [](const type& self,
             const XT::LA::SparsityPatternDefault& local_pattern,
             const ssize_t test_block,
             const ssize_t ansatz_block,
             XT::LA::SparsityPatternDefault& global_pattern) {
            auto tt = XT::Common::numeric_cast<size_t>(test_block);
            auto aa = XT::Common::numeric_cast<size_t>(ansatz_block);
            for (size_t local_ii = 0; local_ii < local_pattern.size(); ++local_ii) {
              const size_t global_ii = self.mapToGlobal(tt, local_ii);
              const auto& local_rows = local_pattern.inner(local_ii);
              for (const auto& local_jj : local_rows) {
                const size_t global_jj = self.mapToGlobal(aa, local_jj);
                global_pattern.insert(global_ii, global_jj);
              }
            }
          },
          "local_sparsity_pattern"_a,
          "test_block"_a,
          "ansatz_block"_a,
          "global_sparsity_pattern"_a);

// matrices
//    addbind_matrix<XT::LA::Backends::common_dense>(c);
//    addbind_matrix<XT::LA::Backends::common_sparse>(c);
//#if HAVE_EIGEN
//    addbind_matrix<XT::LA::Backends::eigen_dense>(c);
//    addbind_matrix<XT::LA::Backends::eigen_sparse>(c);
//#endif
#if HAVE_DUNE_ISTL
    addbind_matrix<XT::LA::Backends::istl_sparse>(c);
#endif

// vectors
//    addbind_vector<XT::LA::Backends::common_dense>(c);
//#if HAVE_EIGEN
//    addbind_vector<XT::LA::Backends::eigen_dense>(c);
//#endif
#if HAVE_DUNE_ISTL
    addbind_vector<XT::LA::Backends::istl_dense>(c);
#endif

    return c;
  } // ... bind(...)
}; // class BlockMapper


template <class SP>
class BlockSpace
{
  typedef typename SP::type S;
  static_assert(is_space<S>::value, "");
  using G = XT::Grid::extract_grid_t<typename S::GridLayerType>;

public:
  typedef GDT::BlockSpace<S> type;
  typedef pybind11::class_<type> bound_type;

private:
  template <bool is_dg = (SP::space_type == SpaceType::dg), bool anything = false>
  struct projector
  {
    template <class V>
    static V project(const type& self, const std::vector<V>& local_vectors, const std::vector<size_t>& subdomains)
    {
      V neighborhood_vector(self.mapper().size(), 0.);
      assert(local_vectors.size() == subdomains.size() && "This should not happen, blame the caller!");
      for (size_t ii = 0; ii < local_vectors.size(); ++ii) {
        const auto subdomain = subdomains[ii];
        const auto& local_vector = local_vectors[ii];
        for (size_t jj = 0; jj < local_vector.size(); ++jj)
          neighborhood_vector[self.mapper().mapToGlobal(subdomain, jj)] = local_vector[jj];
      }
      return neighborhood_vector;
    }
  }; // struct projector<true, ...>

  template <bool anything>
  struct projector<false, anything>
  {
    template <class V>
    static V project(const type&, const std::vector<V>&, const std::vector<size_t>&)
    {
      static_assert(AlwaysFalse<V>::value, "Not implemented for non DG spaces");
    }
  };

  template <XT::LA::Backends la>
  static void addbind_vector(bound_type& c)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    typedef typename XT::LA::Container<double, la>::VectorType V;

    c.def("project_onto_neighborhood",
          [](const type& self, const std::vector<V>& local_vectors, const std::vector<ssize_t>& neighborhood) {
            if (local_vectors.size() != neighborhood.size())
              DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                         "local_vectors.size(): " << local_vectors.size() << "\n   neighborhood.size(): "
                                                  << neighborhood.size());
            std::vector<size_t> subdomains(neighborhood.size());
            for (size_t ii = 0; ii < neighborhood.size(); ++ii)
              subdomains[ii] = XT::Common::numeric_cast<size_t>(neighborhood[ii]);
            std::vector<std::shared_ptr<const S>> neighborhood_spaces(self.dd_grid().size(), nullptr);
            for (const auto& subdomain : subdomains) {
              if (self.backend()[subdomain] == nullptr)
                DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
                           "This BlockSpace (restricted to a neighborhood) does not have a local space for subdomain "
                               << subdomain
                               << "!");
              neighborhood_spaces[subdomain] = self.backend()[subdomain];
            }
            const type neighborhood_space(self.dd_grid(), neighborhood_spaces);
            return projector<>::project(neighborhood_space, local_vectors, subdomains);
          },
          "local_vectors"_a,
          "neighborhood"_a);
  } // ... addbind_vector(...)

public:
  static bound_type bind(pybind11::module& m)
  {
    try {
      BlockMapper<SP>::bind(m);
    } catch (const std::runtime_error&) {
      // already exists
    }

    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = XT::Common::to_camel_case("block_" + space_name<SP>::value());

    bound_type c(m, ClassName.c_str(), ClassName.c_str());

    c.def("__init__",
          [](type& self, XT::Grid::GridProvider<G, XT::Grid::DD::SubdomainGrid<G>>& dd_grid_provider) {
            const auto& dd_grid = dd_grid_provider.dd_grid();
            std::vector<std::shared_ptr<const S>> local_spaces(dd_grid.size());
            for (size_t ss = 0; ss < dd_grid.size(); ++ss)
              local_spaces[ss] = std::make_shared<S>(SP::create(dd_grid_provider, boost::numeric_cast<int>(ss)));
            try {
              new (&self) type(dd_grid, local_spaces);
            } catch (...) {
              self.~type();
              throw;
            }
          },
          "dd_grid"_a,
          py::keep_alive<1, 2>());
    c.def_property_readonly("dimDomain", [](const type& /*self*/) { return S::dimDomain; });
    c.def_property_readonly("dimRange", [](const type& /*self*/) { return S::dimRange; });
    c.def_property_readonly("dimRangeCols", [](const type& /*self*/) { return S::dimRangeCols; });
    c.def_property_readonly("polOrder", [](const type& /*self*/) { return S::polOrder; });
    c.def_property_readonly("num_blocks", [](const type& self) { return self.num_blocks(); });
    c.def_property_readonly("mapper",
                            [](const type& self) {
                              // we get a segfault without the explicit copy
                              return std::decay_t<decltype(self.mapper())>(self.mapper());
                            },
                            py::keep_alive<0, 1>());
    // these need to be defined *after* their non static counterparts
    c.def_property_readonly_static("dimDomain", [](const type& /*self*/) { return S::dimDomain; });
    c.def_property_readonly_static("dimRange", [](const type& /*self*/) { return S::dimRange; });
    c.def_property_readonly_static("dimRangeCols", [](const type& /*self*/) { return S::dimRangeCols; });
    c.def_property_readonly_static("polOrder", [](const type& /*self*/) { return S::polOrder; });
    c.def("local_space",
          [](const type& self, ssize_t block) { return self.local_space(XT::Common::numeric_cast<size_t>(block)); },
          "block"_a);
    c.def("compute_boundary_pattern",
          [](const type& self, const ssize_t block, const std::string tp) {
            auto bb = XT::Common::numeric_cast<size_t>(block);
            if (tp == "volume")
              return self.local_space(bb).compute_volume_pattern(self.dd_grid().boundary_grid_view(bb));
            else if (tp == "face")
              return self.local_space(bb).compute_face_pattern(self.dd_grid().boundary_grid_view(bb));
            else if (tp == "face_and_volume")
              return self.local_space(bb).compute_face_and_volume_pattern(self.dd_grid().boundary_grid_view(bb));
            else
              DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
                         "  type has to be one of ('volume', 'face', 'face_and_volume'), is '" << tp << "'!");
            // we will never get here
            return XT::LA::SparsityPatternDefault();
          },
          "block"_a,
          "type"_a);
    c.def("compute_coupling_pattern",
          [](const type& self, const ssize_t subdomain, const ssize_t neighbor, const std::string tp) {
            auto ss = XT::Common::numeric_cast<size_t>(subdomain);
            auto nn = XT::Common::numeric_cast<size_t>(neighbor);
            if (tp == "volume")
              return self.local_space(ss).compute_volume_pattern(self.dd_grid().coupling_grid_view(ss, nn),
                                                                 self.local_space(nn));
            else if (tp == "face")
              return self.local_space(ss).compute_face_pattern(self.dd_grid().coupling_grid_view(ss, nn),
                                                               self.local_space(nn));
            else if (tp == "face_and_volume")
              return self.local_space(ss).compute_face_and_volume_pattern(self.dd_grid().coupling_grid_view(ss, nn),
                                                                          self.local_space(nn));
            else
              DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
                         "  type has to be one of ('volume', 'face', 'face_and_volume'), is '" << tp << "'!");
            // we will never get here
            return XT::LA::SparsityPatternDefault();
          },
          "block_subdomain"_a,
          "neighboring_subdomain"_a,
          "type"_a);
    c.def("boundary_assembler",
          [](const type& self, const ssize_t subdomain) {
            auto ss = XT::Common::numeric_cast<size_t>(subdomain);
            auto boundary_grid_part = self.dd_grid().boundary_grid_view(ss);
            return new GDT::SystemAssembler<S, decltype(boundary_grid_part), S>(self.local_space(ss), // see below for
                                                                                boundary_grid_part); //  the 'new'
          },
          "subdomain"_a);
    c.def("coupling_assembler",
          [](const type& self, const ssize_t subdomain, const ssize_t neighbor) {
            auto ss = XT::Common::numeric_cast<size_t>(subdomain);
            auto nn = XT::Common::numeric_cast<size_t>(neighbor);
            auto coupling_grid_part = self.dd_grid().coupling_grid_view(ss, nn);
            return new GDT::SystemAssembler<S, decltype(coupling_grid_part), S>(coupling_grid_part, //   SystemAssembler
                                                                                self.local_space(ss), // is not copyable
                                                                                self.local_space(ss), // or movable,
                                                                                self.local_space(nn), // thus the raw
                                                                                self.local_space(nn)); // pointer
          },
          "subdomain"_a,
          "neighbor"_a);
    c.def("restricted_to_neighborhood",
          [](const type& self, const std::vector<ssize_t>& neighborhood) {
            std::vector<std::shared_ptr<const S>> neighborhood_spaces(self.dd_grid().size(), nullptr);
            for (const auto& subdomain : neighborhood)
              neighborhood_spaces[subdomain] = self.backend()[XT::Common::numeric_cast<size_t>(subdomain)];
            return type(self.dd_grid(), neighborhood_spaces);
          },
          "neighborhood"_a);

//    addbind_vector<XT::LA::Backends::common_dense>(c);
//#if HAVE_EIGEN
//    addbind_vector<XT::LA::Backends::eigen_dense>(c);
//#endif
#if HAVE_DUNE_ISTL
    addbind_vector<XT::LA::Backends::istl_dense>(c);
#endif

    const std::string factory_method_name = "make_block_" + space_name<SP>::value_wo_grid();
    m.def(factory_method_name.c_str(),
          [](XT::Grid::GridProvider<G, XT::Grid::DD::SubdomainGrid<G>>& dd_grid_provider) {
            const auto& dd_grid = dd_grid_provider.dd_grid();
            std::vector<std::shared_ptr<const S>> local_spaces(dd_grid.size());
            for (size_t ss = 0; ss < dd_grid.size(); ++ss)
              local_spaces[ss] = std::make_shared<S>(SP::create(dd_grid_provider, boost::numeric_cast<int>(ss)));
            return type(dd_grid, local_spaces);
          },
          "dd_grid_provider"_a);

    return c;
  } // ... bind(...)
}; // class BlockSpace


} // namespace bindings
} // namespace GDT
} // namespace Dune


#endif // PYTHON_DUNE_GDT_PLAYGROUND_SPACES_BLOCK_BINDINGS_HH
