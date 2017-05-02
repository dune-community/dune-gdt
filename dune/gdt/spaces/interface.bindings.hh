// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_SPACES_INTERFACE_BINDINGS_HH
#define DUNE_GDT_SPACES_INTERFACE_BINDINGS_HH
//#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/layers.bindings.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/spaces.hh>

#include "interface.hh"
#include "cg.hh"
#include "dg.hh"
#include "fv.hh"

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
struct backend_name<Backends::fem>
{
  static std::string value()
  {
    return "fem";
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

template <>
struct backend_name<Backends::pdelab>
{
  static std::string value()
  {
    return "pdelab";
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
struct space_type_name<SpaceType::dg>
{
  static std::string value()
  {
    return "dg";
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
struct space_type_name<SpaceType::rt>
{
  static std::string value()
  {
    return "rt";
  }
};


namespace internal {


template <class G, XT::Grid::Layers layer, Backends backend, size_t r, size_t rC>
struct space_name_base
{
  static std::string value_wo_grid()
  {
    using XT::Common::to_string;
    return XT::Grid::bindings::layer_name<layer>::value() + "_to_" + to_string(r) + "x" + to_string(rC) + "_"
           + backend_name<backend>::value();
  }

  static std::string value()
  {
    using XT::Common::to_string;
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

template <class G, XT::Grid::Layers layer, Backends backend, int p, size_t r, size_t rC>
struct space_name<CgSpaceProvider<G, layer, backend, p, double, r, rC>>
{
  static std::string value()
  {
    return std::string("cg_") + internal::space_name_base<G, layer, backend, r, rC>::value() + "_p"
           + XT::Common::to_string(p) + "_space";
  }

  static std::string value_wo_grid()
  {
    return std::string("cg_") + internal::space_name_base<G, layer, backend, r, rC>::value_wo_grid() + "_p"
           + XT::Common::to_string(p) + "_space";
  }
};

template <class G, XT::Grid::Layers layer, Backends backend, int p, size_t r, size_t rC>
struct space_name<DgSpaceProvider<G, layer, backend, p, double, r, rC>>
{
  static std::string value()
  {
    return std::string("dg_") + internal::space_name_base<G, layer, backend, r, rC>::value() + "_p"
           + XT::Common::to_string(p) + "_space";
  }

  static std::string value_wo_grid()
  {
    return std::string("dg_") + internal::space_name_base<G, layer, backend, r, rC>::value_wo_grid() + "_p"
           + XT::Common::to_string(p) + "_space";
  }
};

template <class G, XT::Grid::Layers layer, Backends backend, size_t r, size_t rC>
struct space_name<FvSpaceProvider<G, layer, backend, double, r, rC>>
{
  static std::string value()
  {
    return std::string("fv_") + internal::space_name_base<G, layer, backend, r, rC>::value() + "_space";
  }

  static std::string value_wo_grid()
  {
    return std::string("fv_") + internal::space_name_base<G, layer, backend, r, rC>::value_wo_grid() + "_space";
  }
};

template <class G, XT::Grid::Layers l, SpaceType tp, Backends backend, int p, class R, size_t r, size_t rC>
struct space_name<SpaceProvider<G, l, tp, backend, p, R, r, rC>>
{
  static std::string value()
  {
    return space_type_name<tp>::value() + "_" + internal::space_name_base<G, l, backend, r, rC>::value() + "_p"
           + XT::Common::to_string(p) + "_space";
  }

  static std::string value_wo_grid()
  {
    return space_type_name<tp>::value() + "_" + internal::space_name_base<G, l, backend, r, rC>::value_wo_grid() + "_p"
           + XT::Common::to_string(p) + "_space";
  }
};


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
      using namespace pybind11::literals;
      const std::string factory_method_name = "make_" + space_name<SP>::value_wo_grid();

      m.def(factory_method_name.c_str(),
            [](XT::Grid::GridProvider<G, XT::Grid::DD::SubdomainGrid<G>>& grid_provider, int level) {
              return SP::create(grid_provider, level);
            },
            "grid_provider"_a,
            "level"_a = 0);
    }
  };

  template <bool anything>
  struct factory_methods<false, anything>
  {
    static void addbind(pybind11::module& m)
    {
      using namespace pybind11::literals;
      const std::string factory_method_name = "make_" + space_name<SP>::value_wo_grid();

      m.def(factory_method_name.c_str(),
            [](XT::Grid::GridProvider<G>& grid_provider, int level) { return SP::create(grid_provider, level); },
            "grid_provider"_a,
            "level"_a = 0);
      m.def(factory_method_name.c_str(),
            [](XT::Grid::GridProvider<G, XT::Grid::DD::SubdomainGrid<G>>& grid_provider, int level) {
              return SP::create(grid_provider, level);
            },
            "grid_provider"_a,
            "level"_a = 0);
    }
  };

public:
  typedef S type;
  typedef pybind11::class_<type> bound_type;

  static bound_type bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = XT::Common::to_camel_case(space_name<SP>::value());

    bound_type c(m, ClassName.c_str(), ClassName.c_str(), py::metaclass());

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
    c.def("compute_pattern", [](const type& self) { return self.compute_pattern(); });
    c.def("compute_volume_pattern", [](const type& self) { return self.compute_volume_pattern(); });
    c.def("compute_face_pattern", [](const type& self) { return self.compute_face_pattern(); });
    c.def("compute_face_and_volume_pattern", [](const type& self) { return self.compute_face_and_volume_pattern(); });

    factory_methods<>::addbind(m);

    return c;
  } // ... bind(...)
}; // class SpaceInterface


} // namespace bindings
} // namespace GDT
} // namespace Dune

//#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_SPACES_INTERFACE_BINDINGS_HH
