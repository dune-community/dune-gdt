// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   Rene Milk       (2018)

#ifndef PYTHON_DUNE_GDT_SPACES_INTERFACE_BINDINGS_HH
#define PYTHON_DUNE_GDT_SPACES_INTERFACE_BINDINGS_HH

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/string.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/layers.bindings.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/playground/spaces/restricted.hh>
#include <dune/gdt/spaces.hh>
#include <dune/gdt/type_traits.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/cg.hh>
#include <dune/gdt/spaces/dg.hh>
#include <dune/gdt/spaces/fv.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <Backends backend>
struct backend_name
{
  static_assert(AlwaysFalse<typename internal::backend_dependent_typename<backend>::type>::value,
                "Please add a specialization for this backend!");

  static std::string value()
  {
    return "";
  }
};

template <>
struct backend_name<Backends::gdt>
{
  static std::string value()
  {
    return "gdt";
  }
};


template <SpaceType tp>
struct space_type_name
{
  static_assert(AlwaysFalse<typename internal::space_type_dependent_typename<tp>::type>::value,
                "Please add a specialization for this space type!");

  static std::string value()
  {
    return "";
  }
};

template <>
struct space_type_name<SpaceType::cg>
{
  static std::string value()
  {
    return "cg";
  }
};

template <>
struct space_type_name<SpaceType::block_cg>
{
  static std::string value()
  {
    return "block_cg";
  }
};

template <>
struct space_type_name<SpaceType::dg>
{
  static std::string value()
  {
    return "dg";
  }
};

template <>
struct space_type_name<SpaceType::block_dg>
{
  static std::string value()
  {
    return "block_dg";
  }
};

template <>
struct space_type_name<SpaceType::fv>
{
  static std::string value()
  {
    return "fv";
  }
};

template <>
struct space_type_name<SpaceType::block_fv>
{
  static std::string value()
  {
    return "block_fv";
  }
};

template <>
struct space_type_name<SpaceType::rt>
{
  static std::string value()
  {
    return "rt";
  }
};

template <>
struct space_type_name<SpaceType::block_rt>
{
  static std::string value()
  {
    return "block_rt";
  }
};


namespace internal {


template <class G, XT::Grid::Layers layer, Backends backend, size_t r, size_t rC, XT::Grid::Backends g>
struct space_name_base
{
  static std::string value_wo_grid()
  {
    using XT::Common::to_string;
    return XT::Grid::bindings::layer_name<layer>::value() + "_" + XT::Grid::bindings::backend_name<g>::value() + "_to_"
           + to_string(r) + "x" + to_string(rC) + "_" + backend_name<backend>::value();
  }

  static std::string value()
  {
    return XT::Grid::bindings::grid_name<G>::value() + "_" + value_wo_grid();
  }
};


} // namespace internal


template <class P>
struct space_name
{
  static_assert(AlwaysFalse<P>::value, "Please add a specialization for this space provider!");

  static std::string value()
  {
    return "";
  }
};

template <class G, XT::Grid::Layers layer, Backends backend, int p, size_t r, size_t rC, XT::Grid::Backends g>
struct space_name<CgSpaceProvider<G, layer, backend, p, double, r, rC, g>>
{
  static std::string value()
  {
    return space_type_name<SpaceType::cg>::value() + "_"
           + internal::space_name_base<G, layer, backend, r, rC, g>::value() + "_p" + XT::Common::to_string(p)
           + "_space";
  }

  static std::string value_wo_grid()
  {
    return space_type_name<SpaceType::cg>::value() + "_"
           + internal::space_name_base<G, layer, backend, r, rC, g>::value_wo_grid() + "_p" + XT::Common::to_string(p)
           + "_space";
  }
};

template <class G, XT::Grid::Layers layer, Backends backend, int p, size_t r, size_t rC, XT::Grid::Backends g>
struct space_name<BlockCgSpaceProvider<G, layer, backend, p, double, r, rC, g>>
{
  static std::string value()
  {
    return space_type_name<SpaceType::block_cg>::value() + "_"
           + internal::space_name_base<G, layer, backend, r, rC, g>::value() + "_p" + XT::Common::to_string(p)
           + "_space";
  }

  static std::string value_wo_grid()
  {
    return space_type_name<SpaceType::block_cg>::value() + "_"
           + internal::space_name_base<G, layer, backend, r, rC, g>::value_wo_grid() + "_p" + XT::Common::to_string(p)
           + "_space";
  }
};

template <class G, XT::Grid::Layers layer, Backends backend, int p, size_t r, size_t rC, XT::Grid::Backends g>
struct space_name<DgSpaceProvider<G, layer, backend, p, double, r, rC, g>>
{
  static std::string value()
  {
    return space_type_name<SpaceType::dg>::value() + "_"
           + internal::space_name_base<G, layer, backend, r, rC, g>::value() + "_p" + XT::Common::to_string(p)
           + "_space";
  }

  static std::string value_wo_grid()
  {
    return space_type_name<SpaceType::dg>::value() + "_"
           + internal::space_name_base<G, layer, backend, r, rC, g>::value_wo_grid() + "_p" + XT::Common::to_string(p)
           + "_space";
  }
};

template <class G, XT::Grid::Layers layer, Backends backend, int p, size_t r, size_t rC, XT::Grid::Backends g>
struct space_name<BlockDgSpaceProvider<G, layer, backend, p, double, r, rC, g>>
{
  static std::string value()
  {
    return space_type_name<SpaceType::block_dg>::value() + "_"
           + internal::space_name_base<G, layer, backend, r, rC, g>::value() + "_p" + XT::Common::to_string(p)
           + "_space";
  }

  static std::string value_wo_grid()
  {
    return space_type_name<SpaceType::block_dg>::value() + "_"
           + internal::space_name_base<G, layer, backend, r, rC, g>::value_wo_grid() + "_p" + XT::Common::to_string(p)
           + "_space";
  }
};

template <class G, XT::Grid::Layers layer, Backends backend, size_t r, size_t rC, XT::Grid::Backends g>
struct space_name<FvSpaceProvider<G, layer, backend, double, r, rC, g>>
{
  static std::string value()
  {
    return space_type_name<SpaceType::fv>::value() + "_"
           + internal::space_name_base<G, layer, backend, r, rC, g>::value() + "_space";
  }

  static std::string value_wo_grid()
  {
    return space_type_name<SpaceType::fv>::value() + "_"
           + internal::space_name_base<G, layer, backend, r, rC, g>::value_wo_grid() + "_space";
  }
};

template <class G, XT::Grid::Layers layer, Backends backend, size_t r, size_t rC, XT::Grid::Backends g>
struct space_name<BlockFvSpaceProvider<G, layer, backend, double, r, rC, g>>
{
  static std::string value()
  {
    return space_type_name<SpaceType::block_fv>::value() + "_"
           + internal::space_name_base<G, layer, backend, r, rC, g>::value() + "_space";
  }

  static std::string value_wo_grid()
  {
    return space_type_name<SpaceType::block_fv>::value() + "_"
           + internal::space_name_base<G, layer, backend, r, rC, g>::value_wo_grid() + "_space";
  }
};

template <class G, XT::Grid::Layers layer, Backends backend, int p, size_t r, size_t rC, XT::Grid::Backends g>
struct space_name<RtSpaceProvider<G, layer, backend, p, double, r, rC, g>>
{
  static std::string value()
  {
    return space_type_name<SpaceType::rt>::value() + "_"
           + internal::space_name_base<G, layer, backend, r, rC, g>::value() + "_p" + XT::Common::to_string(p)
           + "_space";
  }

  static std::string value_wo_grid()
  {
    return space_type_name<SpaceType::rt>::value() + "_"
           + internal::space_name_base<G, layer, backend, r, rC, g>::value_wo_grid() + "_p" + XT::Common::to_string(p)
           + "_space";
  }
};

template <class G, XT::Grid::Layers layer, Backends backend, int p, size_t r, size_t rC, XT::Grid::Backends g>
struct space_name<BlockRtSpaceProvider<G, layer, backend, p, double, r, rC, g>>
{
  static std::string value()
  {
    return space_type_name<SpaceType::block_rt>::value() + "_"
           + internal::space_name_base<G, layer, backend, r, rC, g>::value() + "_p" + XT::Common::to_string(p)
           + "_space";
  }

  static std::string value_wo_grid()
  {
    return space_type_name<SpaceType::block_rt>::value() + "_"
           + internal::space_name_base<G, layer, backend, r, rC, g>::value_wo_grid() + "_p" + XT::Common::to_string(p)
           + "_space";
  }
};

template <class G,
          XT::Grid::Layers l,
          SpaceType tp,
          Backends b,
          int p,
          class R,
          size_t r,
          size_t rC,
          XT::Grid::Backends g>
struct space_name<SpaceProvider<G, l, tp, b, p, R, r, rC, g>>
{
  static std::string value()
  {
    return space_type_name<tp>::value() + "_" + internal::space_name_base<G, l, b, r, rC, g>::value() + "_p"
           + XT::Common::to_string(p) + "_space";
  }

  static std::string value_wo_grid()
  {
    return space_type_name<tp>::value() + "_" + internal::space_name_base<G, l, b, r, rC, g>::value_wo_grid() + "_p"
           + XT::Common::to_string(p) + "_space";
  }
};


template <class S>
class SpaceInterfaceWoFactory
{
  static_assert(is_space<S>::value, "");

public:
  typedef S type;
  typedef pybind11::class_<type> bound_type;

  static bound_type bind(pybind11::module& m, const std::string& space_name)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = XT::Common::to_camel_case(space_name /*space_name<SP>::value()*/);

    bound_type c(m, ClassName.c_str(), ClassName.c_str());

    c.def_property_readonly("dimDomain", [](const type& /*self*/) { return S::dimDomain; });
    c.def_property_readonly("dimRange", [](const type& /*self*/) { return S::dimRange; });
    c.def_property_readonly("dimRangeCols", [](const type& /*self*/) { return S::dimRangeCols; });
    c.def_property_readonly("polOrder", [](const type& /*self*/) { return S::polOrder; });
    c.def_property_readonly_static("dimDomain", [](const type& /*self*/) { return S::dimDomain; });
    c.def_property_readonly_static("dimRange", [](const type& /*self*/) { return S::dimRange; });
    c.def_property_readonly_static("dimRangeCols", [](const type& /*self*/) { return S::dimRangeCols; });
    c.def_property_readonly_static("polOrder", [](const type& /*self*/) { return S::polOrder; });

    c.def("size", [](const type& self) { return self.mapper().size(); });
    c.def("visualize",
          [](const type& self, const std::string& filename) { self.visualize(filename); },
          "filename"_a = "");
    c.def("compute_pattern",
          [](const type& self, const std::string tp) {
            if (tp == "default")
              return self.compute_pattern();
            else if (tp == "volume")
              return self.compute_volume_pattern();
            else if (tp == "face")
              return self.compute_face_pattern();
            else if (tp == "face_and_volume")
              return self.compute_face_and_volume_pattern();
            else
              DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
                         "  type has to be one of ('default', volume', 'face', 'face_and_volume'), is '" << tp << "'!");
            // we will never get here
            return XT::LA::SparsityPatternDefault();
          },
          "type"_a = "default");

    return c;
  } // ... bind(...)
}; // class SpaceInterfaceWoFactory


template <class SP>
class SpaceInterface
{
  typedef typename SP::type S;
  static_assert(is_space<S>::value, "");
  using G = XT::Grid::extract_grid_t<typename S::GridLayerType>;

  template <bool is_dd_subdomain_layer = (SP::grid_layer == XT::Grid::Layers::dd_subdomain), bool anything = true>
  struct factory_methods
  {
    static void addbind(pybind11::module& m)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;
      const std::string factory_method_name = "make_" + space_name<SP>::value_wo_grid();

      m.def(factory_method_name.c_str(),
            [](XT::Grid::GridProvider<G, XT::Grid::DD::SubdomainGrid<G>>& grid_provider, int level) {
              return SP::create(grid_provider, level);
            },
            "grid_provider"_a,
            "level"_a = 0,
            py::keep_alive<0, 1>());
    }
  };

  template <bool anything>
  struct factory_methods<false, anything>
  {
    static void addbind(pybind11::module& m)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;
      const std::string factory_method_name = "make_" + space_name<SP>::value_wo_grid();

      m.def(factory_method_name.c_str(),
            [](XT::Grid::GridProvider<G>& grid_provider, int level) { return SP::create(grid_provider, level); },
            "grid_provider"_a,
            "level"_a = 0,
            py::keep_alive<0, 1>());
      m.def(factory_method_name.c_str(),
            [](XT::Grid::GridProvider<G, XT::Grid::DD::SubdomainGrid<G>>& grid_provider, int level) {
              return SP::create(grid_provider, level);
            },
            "grid_provider"_a,
            "level"_a = 0,
            py::keep_alive<0, 1>());
    }
  };

public:
  typedef S type;
  typedef pybind11::class_<type> bound_type;

private:
  template <XT::Grid::Backends backend,
            XT::Grid::Layers layer,
            bool is_dd = (layer == XT::Grid::Layers::dd_subdomain || layer == XT::Grid::Layers::dd_subdomain_boundary
                          || layer == XT::Grid::Layers::dd_subdomain_coupling
                          || layer == XT::Grid::Layers::dd_subdomain_oversampled)>
  struct restriction_methods // true
  {
    static void addbind(bound_type& c)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      typedef GDT::RestrictedSpace<S, typename XT::Grid::Layer<G, layer, backend>::type> RestrictedSpaceType;

      c.def(std::string("restrict_to_" + XT::Grid::bindings::layer_name<layer>::value() + "_"
                        + XT::Grid::bindings::backend_name<backend>::value())
                .c_str(),
            [](type& self,
               XT::Grid::GridProvider<G, XT::Grid::DD::SubdomainGrid<G>>& grid_provider,
               const int level_or_subdomain = -1) {
              return RestrictedSpaceType(self, grid_provider.template layer<layer, backend>(level_or_subdomain));
            },
            "grid_provider"_a,
            "level_or_subdomain"_a,
            py::keep_alive<1, 0>());
    }
  }; // struct restriction_methods

  template <XT::Grid::Backends backend, XT::Grid::Layers layer>
  struct restriction_methods<backend, layer, false>
  {
    static void addbind(bound_type& c)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      typedef GDT::RestrictedSpace<S, typename XT::Grid::Layer<G, layer, backend>::type> RestrictedSpaceType;

      c.def(std::string("restrict_to_" + XT::Grid::bindings::layer_name<layer>::value() + "_"
                        + XT::Grid::bindings::backend_name<backend>::value())
                .c_str(),
            [](type& self, XT::Grid::GridProvider<G>& grid_provider, const int level = -1) {
              return RestrictedSpaceType(self, grid_provider.template layer<layer, backend>(level));
            },
            "grid_provider"_a,
            "level"_a,
            py::keep_alive<1, 0>());
      restriction_methods<backend, layer, true>::addbind(c);
    }
  }; // struct restriction_methods<..., false>

  template <XT::Grid::Backends backend, XT::Grid::Layers layer>
  static void addbind_restricted(pybind11::module& m, bound_type& c, const std::string sp_name)
  {
    using namespace pybind11::literals;

    typedef GDT::RestrictedSpace<S, typename XT::Grid::Layer<G, layer, backend>::type> RestrictedSpaceType;

    XT::Common::bindings::try_register(m, [&](pybind11::module& mod) {
      const auto restricted_space_name = sp_name + "_restricted_to_" + XT::Grid::bindings::layer_name<layer>::value()
                                         + "_" + XT::Grid::bindings::backend_name<backend>::value();
      auto restricted_space = SpaceInterfaceWoFactory<RestrictedSpaceType>::bind(mod, restricted_space_name);
      restricted_space.def("restrict",
                           [](RestrictedSpaceType& self, const XT::LA::IstlDenseVector<double>& unrestricted_vector) {
                             return self.mapper().restrict(unrestricted_vector);
                           },
                           "unrestricted_vector"_a);
      restricted_space.def("extend",
                           [](RestrictedSpaceType& self, const XT::LA::IstlDenseVector<double>& restricted_vector) {
                             return self.mapper().extend(restricted_vector);
                           },
                           "restricted_vector"_a);
    });

    restriction_methods<backend, layer>::addbind(c);
  } // ... addbind_restricted(...)

public:
  static bound_type bind(pybind11::module& m)
  {
    const auto sp_name = space_name<SP>::value();
    auto c = SpaceInterfaceWoFactory<S>::bind(m, sp_name);

    addbind_restricted<XT::Grid::Backends::view, XT::Grid::Layers::dd_subdomain>(m, c, sp_name);
    addbind_restricted<XT::Grid::Backends::view, XT::Grid::Layers::dd_subdomain_boundary>(m, c, sp_name);
    addbind_restricted<XT::Grid::Backends::view, XT::Grid::Layers::dd_subdomain_coupling>(m, c, sp_name);
    addbind_restricted<XT::Grid::Backends::view, XT::Grid::Layers::dd_subdomain_oversampled>(m, c, sp_name);
    addbind_restricted<XT::Grid::Backends::view, XT::Grid::Layers::leaf>(m, c, sp_name);
    addbind_restricted<XT::Grid::Backends::view, XT::Grid::Layers::level>(m, c, sp_name);

    factory_methods<>::addbind(m);
    return c;
  } // ... bind(...)
}; // class SpaceInterface


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // PYTHON_DUNE_GDT_SPACES_INTERFACE_BINDINGS_HH
