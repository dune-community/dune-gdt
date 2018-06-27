// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   Rene Milk       (2018)

#ifndef PYTHON_DUNE_GDT_ASSEMBLER_SYSTEM_BINDINGS_HH
#define PYTHON_DUNE_GDT_ASSEMBLER_SYSTEM_BINDINGS_HH

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/common/memory.hh>
#include <dune/xt/la/container.hh>
#include <dune/xt/grid/grids.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/layers.bindings.hh>
#include <python/dune/xt/grid/walker.bindings.hh>

#include <python/dune/gdt/spaces/interface.hh>
#include <python/dune/gdt/spaces/constraints.hh>
#include <dune/gdt/type_traits.hh>

#include <dune/gdt/assembler/system.hh>

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


namespace internal {


template <XT::LA::Backends la, class type, class... Bases>
static void addaddbind_matrixatrix(pybind11::class_<type, Bases...>& c)
{
  namespace py = pybind11;
  using namespace pybind11::literals;
  using T = typename type::TestSpaceType;
  using GL = typename type::GridLayerType;

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

  // add constraints
  bindings::DirichletConstraints<XT::Grid::extract_intersection_t<typename type::GridLayerType>,
                                 XT::Grid::extract_grid_t<typename type::GridLayerType>>::addbind(c);


} // ... addaddbind_matrixatrix(...)


template <class AssemblerOrDerivedType, class... Bases>
static void bind_system_assembler_functions(pybind11::class_<AssemblerOrDerivedType, Bases...>& c)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  using type = AssemblerOrDerivedType;
  using GL = typename type::GridLayerType;
  XT::Grid::bindings::internal::bind_walker_functions(c);

  c.def(
      "append",
      [](type& self, XT::Grid::Walker<GL>& other, const XT::Grid::ApplyOn::WhichIntersection<GL>& which_intersections) {
        self.append(other, which_intersections.copy());
      },
      "grid_walker"_a,
      "which_intersections"_a = XT::Grid::ApplyOn::AllIntersections<GL>(),
      py::keep_alive<0, 1>(),
      py::keep_alive<0, 2>());
  // add rest
  c.def("append", [](type& self, type& other) { self.append(other); }, "system_assembler"_a, py::keep_alive<1, 2>());

  c.def("append",
        [](type& self, type& other, const XT::Grid::ApplyOn::WhichIntersection<GL>& which_intersections) {
          self.append(other, which_intersections.copy());
        },
        "system_assembler"_a,
        "which_intersections"_a,
        py::keep_alive<1, 2>());
  c.def("assemble",
        [](type& self, const bool use_tbb) {
          py::gil_scoped_release DUNE_UNUSED(release);
          self.assemble(use_tbb);
        },
        "use_tbb"_a = false);

  using G = XT::Grid::extract_grid_t<GL>;
  using I = XT::Grid::extract_intersection_t<GL>;
  typedef XT::Grid::extract_entity_t<GL> E;
  typedef typename G::ctype D;
  static const constexpr size_t d = G::dimension;
  bindings::DirichletConstraints<I, G>::addbind(c);
  addaddbind_matrixatrix<XT::LA::Backends::istl_sparse>(c);

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
}

template <class T, class GL = typename T::GridLayerType, class A = T>
class SystemAssembler
{
  static_assert(is_space<T>::value, "");
  static_assert(XT::Grid::is_layer<GL>::value, "");
  static_assert(is_space<A>::value, "");
  typedef XT::Grid::extract_grid_t<GL> G;
  typedef XT::Grid::extract_entity_t<GL> E;
  typedef typename G::ctype D;
  static const constexpr size_t d = G::dimension;

public:
  typedef GDT::SystemAssembler<T, GL, A> type;
  typedef pybind11::class_<type, XT::Grid::Walker<GL>> bound_type;

private:
  typedef typename type::TestSpaceType TestSpaceType;
  typedef typename type::GridLayerType GridLayerType;
  typedef typename type::AnsatzSpaceType AnsatzSpaceType;

  template <bool same_spaces = std::is_same<TestSpaceType, AnsatzSpaceType>::value,
            bool same_layer = std::is_same<GridLayerType, typename TestSpaceType::GridLayerType>::value&&
                std::is_same<GridLayerType, typename AnsatzSpaceType::GridLayerType>::value,
            bool anything = true>
  struct addbind_switch
  {
    static void ctors(bound_type& /*c*/)
    {
    }

    static void factory_methods(pybind11::module& /*m*/)
    {
    }
  };

  template <bool anything>
  struct addbind_switch<true, true, anything>
  {
    static void ctors(bound_type& c)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      c.def(py::init<T>(),
            "space"_a,
            "Uses space as test and ansatz space, and its grid layer as grid layer",
            py::keep_alive<1, 2>());
    }

    static void factory_methods(pybind11::module& m)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      m.def("make_system_assembler",
            [](TestSpaceType space) { return new type(space); },
            "space"_a,
            "Uses space as test and ansatz space, and its grid layer as grid layer",
            py::keep_alive<0, 1>());

      addbind_switch<false, true>::factory_methods(m);
    }
  }; // struct factory_methods<true, true, ...>

  template <bool anything>
  struct addbind_switch<false, true, anything>
  {
    static void ctors(bound_type& /*c*/)
    {
    }

    static void factory_methods(pybind11::module& m)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      m.def("make_system_assembler",
            [](const TestSpaceType& test_space, const TestSpaceType& ansatz_space) {
              return new type(test_space, ansatz_space);
            },
            "test_space"_a,
            "ansatz_space"_a,
            "Uses the grid layer of the test_space as grid layer.",
            py::keep_alive<0, 1>(),
            py::keep_alive<0, 2>());
    }
  }; // struct factory_methods<false, true, ...>

public:
  static bound_type bind(pybind11::module& m,
                         const std::string& test_space_name,
                         const std::string& ansatz_space_name,
                         const std::string& grid_layer_name)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    // add class
    std::string class_name = "system_assembler";
    class_name += "_wrt_" + test_space_name;
    if (ansatz_space_name != test_space_name)
      class_name += "_and_" + ansatz_space_name;
    class_name += "_on_" + grid_layer_name;

    bound_type c(m, XT::Common::to_camel_case(class_name).c_str(), XT::Common::to_camel_case(class_name).c_str());
    bind_system_assembler_functions(c);
    // add ctor
    addbind_switch<>::ctors(c);
    // add factory methods
    addbind_switch<>::factory_methods(m);
    // finished
    return c;
  } // ... bind(...)
}; // class SystemAssembler


} // namespace internal


template <class TP, XT::Grid::Layers grid_layer, XT::Grid::Backends grid_backend, class AP>
class SystemAssembler
{
  typedef typename TP::type T;
  static_assert(is_space<T>::value, "");
  typedef typename AP::type A;
  static_assert(is_space<A>::value, "");
  typedef XT::Grid::extract_grid_t<typename T::GridLayerType> G;
  typedef typename XT::Grid::Layer<G, grid_layer, grid_backend, XT::Grid::DD::SubdomainGrid<G>>::type GL;
  using binder = internal::SystemAssembler<T, GL, A>;

public:
  using type = typename binder::type;
  using bound_type = typename pybind11::class_<type>;

  static bound_type bind(pybind11::module& m)
  {
    return binder::bind(m,
                        space_name<TP>::value(),
                        space_name<AP>::value(),
                        XT::Grid::bindings::layer_name<grid_layer>::value() + "_"
                            + XT::Grid::bindings::backend_name<grid_backend>::value());
  }
}; // class SystemAssembler
} // namespace bindings
} // namespace GDT
} // namespace Dune


#endif // PYTHON_DUNE_GDT_ASSEMBLER_SYSTEM_BINDINGS_HH
