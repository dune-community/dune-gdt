// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_OPERATORS_BASE_BINDINGS_HH
#define DUNE_GDT_OPERATORS_BASE_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/la/container.hh>

#include <dune/gdt/type_traits.hh>
#include <dune/gdt/assembler/system.bindings.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace bindings {


template <class OperatorType>
class MatrixOperatorBase
{
  static_assert(is_matrix_operator<OperatorType>::value, "");

public:
  typedef OperatorType type;
  using R = typename OperatorType::RangeSpaceType;
  using S = typename OperatorType::RangeSpaceType;
  using GL = typename OperatorType::GridLayerType;
  typedef GDT::SystemAssembler<R, GL, S> BaseType;
  typedef pybind11::class_<type, BaseType> bound_type;

private:
  typedef typename XT::LA::Container<typename type::FieldType, type::MatrixType::vector_type>::VectorType V;

  template <bool same_spaces = std::is_same<R, S>::value,
            bool same_layer =
                std::is_same<GL, typename R::GridLayerType>::value&& std::is_same<GL, typename S::GridLayerType>::value,
            bool anything = true>
  struct addbind_switch
  {
    static void induced_norm(bound_type& /*c*/)
    {
    }

    static void pattern(bound_type& /*c*/)
    {
    }
  };

  template <bool anything>
  struct addbind_switch<true, true, anything>
  {
    static void induced_norm(bound_type& c)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      c.def("induced_norm",
            [](type& self, const V& range) {
              py::gil_scoped_release DUNE_UNUSED(release);
              return self.induced_norm(range);
            },
            "range"_a);
      c.def("induced_norm",
            [](type& self, const GDT::ConstDiscreteFunction<S, V>& range) {
              py::gil_scoped_release DUNE_UNUSED(release);
              return self.induced_norm(range);
            },
            "range"_a);
    } // ... induced_norm(...)

    static void pattern(bound_type& c)
    {
      c.def_static("pattern", [](const R& space) { return type::pattern(space); });
      addbind_switch<false, true>::induced_norm(c);
    }
  }; // struct addbind_switch<true, true, ...>

  template <bool anything>
  struct addbind_switch<false, true, anything>
  {
    static void induced_norm(bound_type& /*c*/)
    {
    }

    static void pattern(bound_type& c)
    {
      c.def_static("pattern", [](const R& range_space, const S& source_space) {
        return type::pattern(range_space, source_space);
      });
    }
  }; // struct addbind_switch<true, true, ...>

public:
  static bound_type bind(pybind11::module& m,
                         const std::string& class_id,
                         const std::string& test_space_name,
                         const std::string& ansatz_space_name,
                         const std::string& grid_layer_name)
  {
    try { //  we might not be the first to add this
      internal::SystemAssembler<R, GL, S>::bind(m, test_space_name, ansatz_space_name, grid_layer_name);
    } catch (const std::runtime_error&) {
    }

    namespace py = pybind11;
    using namespace pybind11::literals;

    bound_type c(m, std::string(class_id).c_str(), std::string(class_id).c_str());
    // from MatrixOperatorBase
    addbind_switch<>::pattern(c);
    c.def("pattern",
          [](type& self) { return self.pattern(self.range_space(), self.source_space(), self.grid_layer()); });
    c.def("matrix", [](type& self) { return self.matrix(); });
    c.def("source_space", [](type& self) { return self.source_space(); });
    c.def("range_space", [](type& self) { return self.range_space(); });
    c.def("apply",
          [](type& self, const V& source, V& range) {
            py::gil_scoped_release DUNE_UNUSED(release);
            self.apply(source, range);
          },
          "source"_a,
          "range"_a);
    c.def("apply",
          [](type& self, const GDT::ConstDiscreteFunction<S, V>& source, GDT::DiscreteFunction<R, V>& range) {
            py::gil_scoped_release DUNE_UNUSED(release);
            self.apply(source, range);
          },
          "source"_a,
          "range"_a);
    c.def("apply2",
          [](type& self, const V& range, const V& source) {
            py::gil_scoped_release DUNE_UNUSED(release);
            return self.apply2(range, source);
          },
          "range"_a,
          "source"_a);
    c.def(
        "apply2",
        [](type& self, const GDT::ConstDiscreteFunction<R, V>& range, const GDT::ConstDiscreteFunction<S, V>& source) {
          py::gil_scoped_release DUNE_UNUSED(release);
          return self.apply2(range, source);
        },
        "range"_a,
        "source"_a);
    c.def("apply_inverse",
          [](type& self, const V& range, V& source, const XT::Common::Configuration& opts) {
            py::gil_scoped_release DUNE_UNUSED(release);
            self.apply_inverse(range, source, opts);
          },
          "range"_a,
          "source"_a,
          "opts"_a);
    c.def("apply_inverse",
          [](type& self,
             const GDT::ConstDiscreteFunction<R, V>& range,
             GDT::ConstDiscreteFunction<S, V>& source,
             const XT::Common::Configuration& opts) {
            py::gil_scoped_release DUNE_UNUSED(release);
            self.apply_inverse(range, source, opts);
          },
          "range"_a,
          "source"_a,
          "opts"_a);
    c.def("invert_options", [](type& self) { return self.invert_options(); });
    c.def("invert_options", [](type& self, const std::string& type) { return self.invert_options(type); }, "type"_a);

    // from OperatorInterface
    c.def("apply_inverse",
          [](type& self, const V& range, V& source, const std::string& type) {
            py::gil_scoped_release DUNE_UNUSED(release);
            self.apply_inverse(range, source, type);
          },
          "range"_a,
          "source"_a,
          "type"_a);
    c.def("apply_inverse",
          [](type& self,
             const GDT::ConstDiscreteFunction<R, V>& range,
             GDT::ConstDiscreteFunction<S, V>& source,
             const std::string& type) {
            py::gil_scoped_release DUNE_UNUSED(release);
            self.apply_inverse(range, source, type);
          },
          "range"_a,
          "source"_a,
          "type"_a);
    c.def("apply_inverse",
          [](type& self, const V& range, V& source) {
            py::gil_scoped_release DUNE_UNUSED(release);
            self.apply_inverse(range, source);
          },
          "range"_a,
          "source"_a);
    c.def("apply_inverse",
          [](type& self, const GDT::ConstDiscreteFunction<R, V>& range, GDT::ConstDiscreteFunction<S, V>& source) {
            py::gil_scoped_release DUNE_UNUSED(release);
            self.apply_inverse(range, source);
          },
          "range"_a,
          "source"_a);
    addbind_switch<>::induced_norm(c);

    return c;
  } // ... bind(...)
}; // class MatrixOperatorBase


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_OPERATORS_BASE_BINDINGS_HH
