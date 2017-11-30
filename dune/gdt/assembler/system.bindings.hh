// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_ASSEMBLER_SYSTEM_BINDINGS_HH
#define DUNE_GDT_ASSEMBLER_SYSTEM_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/memory.hh>
#include <dune/xt/la/container.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/layers.bindings.hh>
#include <dune/xt/grid/walker.hh>

#include <dune/gdt/spaces.bindings.hh>
#include <dune/gdt/spaces/constraints.bindings.hh>
#include <dune/gdt/type_traits.hh>

#include "system.hh"

namespace Dune {
namespace GDT {
namespace bindings {


class ResultStorage
{
public:
  ResultStorage()
    : result_(0.)
  {
  }

  ResultStorage(const ResultStorage& other) = delete;
  ResultStorage(ResultStorage&& source) = delete;

  ResultStorage& operator=(const ResultStorage& other) = delete;
  ResultStorage& operator=(ResultStorage&& source) = delete;

  double& result()
  {
    return result_;
  }

  const double& result() const
  {
    return result_;
  }

private:
  double result_;
}; // class ResultStorage


template <class TP, XT::Grid::Layers grid_layer, XT::Grid::Backends grid_backend>
class SystemAssembler
{
  typedef typename TP::type T;
  static_assert(is_space<T>::value, "");
  typedef XT::Grid::extract_grid_t<typename T::GridLayerType> G;
  typedef typename XT::Grid::Layer<G, grid_layer, grid_backend, XT::Grid::DD::SubdomainGrid<G>>::type GL;
  typedef XT::Grid::extract_entity_t<GL> E;
  typedef typename G::ctype D;
  static const constexpr size_t d = G::dimension;

public:
  typedef GDT::SystemAssembler<T, GL> type;
  typedef pybind11::class_<type, XT::Grid::Walker<GL>> bound_type;

private:
  typedef typename type::TestSpaceType TestSpaceType;
  typedef typename type::GridLayerType GridLayerType;
  typedef typename type::AnsatzSpaceType AnsatzSpaceType;

  template <bool do_bind = (std::is_same<TestSpaceType, AnsatzSpaceType>::value
                            && std::is_same<GridLayerType, typename TestSpaceType::GridLayerType>::value),
            bool anything = true>
  struct addbind_ctor_single
  {
    void operator()(bound_type& c)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;
      c.def(py::init<T>(),
            "space"_a,
            "Uses given space as test and ansatz space, and the grid layer of the given space as grid layer.",
            py::keep_alive<1, 2>());
    }
  }; // struct addbind_ctor_single

  template <bool anything>
  struct addbind_ctor_single<false, anything>
  {
    void operator()(bound_type& /*c*/)
    {
    }
  }; // struct addbind_ctor_single

  template <bool same = std::is_same<GL, typename T::GridLayerType>::value, bool anything = true>
  struct space_and_layer_name
  {
    static std::string value()
    {
      return space_name<TP>::value();
    }
  };

  template <bool anything>
  struct space_and_layer_name<false, anything>
  {
    static std::string value()
    {
      return space_name<TP>::value() + "_" + XT::Grid::bindings::layer_name<grid_layer>::value() + "_"
             + XT::Grid::bindings::backend_name<grid_backend>::value();
    }
  };

  template <bool same = std::is_same<GL, typename T::GridLayerType>::value, bool anything = true>
  struct addbind_factory_methods
  {
    void operator()(pybind11::module& m)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      m.def("make_system_assembler",
            [](const TestSpaceType& space) { return new type(space); },
            "space"_a,
            py::keep_alive<0, 1>());
    }
  };

  template <bool anything>
  struct addbind_factory_methods<false, anything>
  {
    void operator()(pybind11::module& /*m*/)
    {
    }
  };

  template <XT::LA::Backends la>
  static void addaddbind_matrixatrix(bound_type& c)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    typedef typename XT::LA::Container<typename T::RangeFieldType, la>::MatrixType M;

    c.def("append",
          [](type& self,
             const GDT::LocalBoundaryTwoFormInterface<typename T::BaseFunctionSetType,
                                                      XT::Grid::extract_intersection_t<GL>>& local_boundary_two_form,
             M& matrix,
             const XT::Grid::ApplyOn::WhichIntersection<GL>& which_intersections) {
            self.append(local_boundary_two_form, matrix, which_intersections.copy());
          },
          "local_boundary_two_form"_a,
          "matrix"_a,
          "which_intersections"_a = XT::Grid::ApplyOn::AllIntersections<GL>(),
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());
    c.def("append",
          [](type& self,
             const GDT::LocalCouplingTwoFormInterface<typename T::BaseFunctionSetType,
                                                      XT::Grid::extract_intersection_t<GL>>& local_coupling_two_form,
             M& matrix,
             const XT::Grid::ApplyOn::WhichIntersection<GL>& which_intersections) {
            self.append(local_coupling_two_form, matrix, which_intersections.copy());
          },
          "local_coupling_two_form"_a,
          "matrix"_a,
          "which_intersections"_a = XT::Grid::ApplyOn::AllIntersections<GL>(),
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>());
    c.def("append",
          [](type& self,
             const GDT::LocalCouplingTwoFormInterface<typename T::BaseFunctionSetType,
                                                      XT::Grid::extract_intersection_t<GL>>& local_coupling_two_form,
             M& matrix_in_in,
             M& matrix_out_out,
             M& matrix_in_out,
             M& matrix_out_in,
             const XT::Grid::ApplyOn::WhichIntersection<GL>& which_intersections) {
            self.append(local_coupling_two_form,
                        matrix_in_in,
                        matrix_out_out,
                        matrix_in_out,
                        matrix_out_in,
                        which_intersections.copy());
          },
          "local_coupling_two_form"_a,
          "matrix_in_in"_a,
          "matrix_out_out"_a,
          "matrix_in_out"_a,
          "matrix_out_in"_a,
          "which_intersections"_a = XT::Grid::ApplyOn::AllIntersections<GL>(),
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>(),
          py::keep_alive<0, 4>(),
          py::keep_alive<0, 5>());
  } // ... addaddbind_matrixatrix(...)

public:
  static bound_type bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName =
        XT::Common::to_camel_case(std::string("system_assembler_") + space_and_layer_name<>::value());

    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    addbind_ctor_single<>()(c);

    c.def("append",
          [](type& self, type& other, const XT::Grid::ApplyOn::WhichIntersection<GL>& which_intersections) {
            self.append(other, which_intersections.copy());
          },
          "system_assembler"_a,
          "which_intersections"_a = XT::Grid::ApplyOn::AllIntersections<GL>(),
          py::keep_alive<1, 2>());
    c.def("append",
          [](type& self,
             XT::Grid::Walker<GridLayerType>& other,
             const XT::Grid::ApplyOn::WhichIntersection<GL>& which_intersections) {
            self.append(other, which_intersections.copy());
          },
          "grid_walker"_a,
          "which_intersections"_a = XT::Grid::ApplyOn::AllIntersections<GL>(),
          py::keep_alive<1, 2>());
    c.def("assemble",
          [](type& self, const bool use_tbb) {
            py::gil_scoped_release DUNE_UNUSED(release);
            self.assemble(use_tbb);
          },
          "use_tbb"_a = false);

    bindings::DirichletConstraints<XT::Grid::extract_intersection_t<typename type::GridLayerType>,
                                   XT::Grid::extract_grid_t<typename type::GridLayerType>>::addbind(c);

#if HAVE_DUNE_ISTL
    addaddbind_matrixatrix<XT::LA::Backends::istl_sparse>(c);
#endif

    c.def("append",
          [](type& self,
             const GDT::LocalVolumeTwoFormInterface<XT::Functions::LocalfunctionInterface<E, D, d, double, 1>,
                                                    XT::Functions::LocalfunctionInterface<E, D, d, double, 1>,
                                                    double>& local_volume_two_form,
             const XT::Functions::LocalizableFunctionInterface<E, D, d, double, 1>& test_function,
             const XT::Functions::LocalizableFunctionInterface<E, D, d, double, 1>& ansatz_function,
             ResultStorage& result /*,
             const XT::Grid::ApplyOn::WhichEntity<GL>& where*/) {
            self.append(local_volume_two_form, test_function, ansatz_function, result.result() /*, where.copy()*/);
          },
          "local_volume_two_form"_a,
          "test_function"_a,
          "ansatz_function"_a,
          "result"_a /*,
          "where"_a = XT::Grid::ApplyOn::AllEntities<GL>()*/,
          py::keep_alive<0, 1>(),
          py::keep_alive<0, 2>(),
          py::keep_alive<0, 3>());

    addbind_factory_methods<>()(m);

    return c;
  } // ... bind(...)
}; // class SystemAssembler


} // namespace bindings
} // namespace GDT
} // namespace Dune






  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_LIB_ALU_SPACE(_pre, leaf, view, _s_type, gdt, leaf, _p, 1, 1);                       \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_LIB_ALU(_pre, fv, 0, 1, 1);                                                          \
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_LIB_ALU(_pre, cg, 1, 1, 1)
  _DUNE_GDT_ASSEMBLER_SYSTEM_BIND_LIB_YASP(_pre, cg, 1, 1, 1)



#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_ASSEMBLER_SYSTEM_BINDINGS_HH
