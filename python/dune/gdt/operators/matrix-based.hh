// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef DUNE_GDT_PYTHON_OPERATORS_MATRIX_BASED_HH
#define DUNE_GDT_PYTHON_OPERATORS_MATRIX_BASED_HH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/functors/interfaces.hh>
#include <dune/xt/grid/gridprovider/coupling.hh>
#include <dune/xt/grid/gridprovider/provider.hh>

#include <dune/xt/la/type_traits.hh>

#include <dune/gdt/operators/matrix.hh>
#include <dune/gdt/tools/sparsity-pattern.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/parameter.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/grid/walker.hh>
#include <python/dune/xt/la/traits.hh>
#include <python/dune/gdt/operators/interfaces.hh>


namespace Dune {
namespace GDT {
namespace bindings {


template <class M,
          class MatrixTag,
          class SparsityTag,
          class AGV,
          size_t s = 1,
          size_t r = s,
          class SGV = AGV,
          class RGV = AGV>
class MatrixOperator
{
  using G = std::decay_t<XT::Grid::extract_grid_t<AGV>>;
  static const size_t d = G::dimension;
  using GP = XT::Grid::GridProvider<G>;

public:
  using type = GDT::MatrixOperator<AGV, s, 1, r, 1, double, M, SGV, RGV>;
  using base_operator_type = GDT::ConstMatrixOperator<AGV, s, 1, r, 1, double, M, SGV, RGV>;
  using base_functor_type = Dune::XT::Grid::ElementAndIntersectionFunctor<AGV>;
  using bound_type = pybind11::class_<type, base_operator_type, base_functor_type>;

private:
  using SS = typename type::SourceSpaceType;
  using RS = typename type::RangeSpaceType;
  using E = typename type::E;
  using I = typename type::I;
  using F = typename type::F;
  using BilinearFormType = typename type::BilinearFormType;

  template <bool needs_sparsity_tag = !std::is_same<SparsityTag, void>::value, bool anything = true>
  struct addbind /*<true, ...>*/
  {
    static void leaf_factory(pybind11::module& m, const std::string& FactoryName)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      //      m.def(
      //          FactoryName.c_str(),
      //          [](GP& grid,
      //             const SS& source_space,
      //             const RS& range_space,
      //             const MatrixTag&,
      //             const SparsityTag&,
      //             const XT::LA::SparsityPatternDefault& pattern,
      //             const std::string& logging_prefix) {
      //            return new type(grid.leaf_view(), source_space, range_space, pattern, logging_prefix);
      //          },
      //          "grid"_a,
      //          "source_space"_a,
      //          "range_space"_a,
      //          "la_backend"_a,
      //          "sparsity_type"_a,
      //          "sparsity_pattern"_a,
      //          "logging_prefix"_a = "",
      //          py::keep_alive<0, 1>(),
      //          py::keep_alive<0, 2>(),
      //          py::keep_alive<0, 3>());
      //      m.def(
      //          FactoryName.c_str(),
      //          [](GP& grid,
      //             const SS& source_space,
      //             const RS& range_space,
      //             const MatrixTag&,
      //             const SparsityTag&,
      //             const std::string& logging_prefix) {
      //            return new type(grid.leaf_view(),
      //                            source_space,
      //                            range_space,
      //                            make_element_and_intersection_sparsity_pattern(range_space, source_space,
      //                            grid.leaf_view()), logging_prefix);
      //          },
      //          "grid"_a,
      //          "source_space"_a,
      //          "range_space"_a,
      //          "la_backend"_a,
      //          "sparsity_type"_a,
      //          "logging_prefix"_a = "",
      //          py::keep_alive<0, 1>(),
      //          py::keep_alive<0, 2>(),
      //          py::keep_alive<0, 3>());
    }
    static void coupling_factory(pybind11::module& m, const std::string& FactoryName)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;
    }
  }; // struct addbind<true, ...>

  template <bool anything>
  struct addbind<false, anything>
  {
    static void leaf_factory(pybind11::module& m, const std::string& FactoryName)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;

      m.def(
          FactoryName.c_str(),
          [](GP& grid,
             const SS& source_space,
             const RS& range_space,
             const MatrixTag&,
             const XT::LA::SparsityPatternDefault& pattern,
             const std::string& logging_prefix) {
            return new type(grid.leaf_view(),
                            source_space,
                            range_space,
                            new M(range_space.mapper().size(), source_space.mapper().size(), pattern),
                            logging_prefix);
          },
          "grid"_a,
          "source_space"_a,
          "range_space"_a,
          "la_backend"_a,
          "sparsity_pattern"_a,
          "logging_prefix"_a = "",
          py::keep_alive<0, 1>()
          //          py::keep_alive<0, 2>()
      );
      if (std::is_same<MatrixTag, XT::LA::bindings::Istl>::value) {
        m.def(
            FactoryName.c_str(),
            [](GP& grid,
               const SS& source_space,
               const RS& range_space,
               const XT::LA::SparsityPatternDefault& pattern,
               const MatrixTag&,
               const std::string& logging_prefix) {
              return new type(grid.leaf_view(),
                              source_space,
                              range_space,
                              new M(range_space.mapper().size(), source_space.mapper().size(), pattern),
                              logging_prefix);
            },
            "grid"_a,
            "source_space"_a,
            "range_space"_a,
            "sparsity_pattern"_a,
            "la_backend"_a = MatrixTag(),
            "logging_prefix"_a = "",
            py::keep_alive<0, 1>()
            //            py::keep_alive<0, 2>()
        );
        m.def(
            FactoryName.c_str(),
            [](GP& grid,
               const SS& source_space,
               const RS& range_space,
               const MatrixTag&,
               const std::string& logging_prefix) {
              return new type(
                  grid.leaf_view(),
                  source_space,
                  range_space,
                  new M(range_space.mapper().size(),
                        source_space.mapper().size(),
                        make_element_and_intersection_sparsity_pattern(range_space, source_space, grid.leaf_view())),
                  logging_prefix);
            },
            "grid"_a,
            "source_space"_a,
            "range_space"_a,
            "la_backend"_a = MatrixTag(),
            "logging_prefix"_a = "",
            py::keep_alive<0, 1>()
            //            py::keep_alive<0, 2>()
        );
      }
    }

    static void coupling_factory(pybind11::module& m, const std::string& FactoryName)
    {
      namespace py = pybind11;
      using namespace pybind11::literals;
      using CGP = XT::Grid::CouplingGridProvider<AGV>;
//      m.def(
//          FactoryName.c_str(),
//          [](CGP& grid,
//             const SS& source_space,
//             const RS& range_space,
//             const MatrixTag&,
//             const XT::LA::SparsityPatternDefault& pattern,
//             const std::string& logging_prefix) {
//            return new type(grid.coupling_view(),
//                            source_space,
//                            range_space,
//                            new M(range_space.mapper().size(), source_space.mapper().size(), pattern),
//                            logging_prefix);
//          },
//          "grid"_a,
//          "source_space"_a,
//          "range_space"_a,
//          "la_backend"_a,
//          "sparsity_pattern"_a,
//          "logging_prefix"_a = "",
//          py::keep_alive<0, 1>()
//          //          py::keep_alive<0, 2>()
//      );
      if (std::is_same<MatrixTag, XT::LA::bindings::Istl>::value) {
//        m.def(
//            FactoryName.c_str(),
//            [](CGP& grid,
//               const SS& source_space,
//               const RS& range_space,
//               const XT::LA::SparsityPatternDefault& pattern,
//               const MatrixTag&,
//               const std::string& logging_prefix) {
//              return new type(grid.coupling_view(),
//                              source_space,
//                              range_space,
//                              new M(range_space.mapper().size(), source_space.mapper().size(), pattern),
//                              logging_prefix);
//            },
//            "grid"_a,
//            "source_space"_a,
//            "range_space"_a,
//            "sparsity_pattern"_a,
//            "la_backend"_a = MatrixTag(),
//            "logging_prefix"_a = "",
//            py::keep_alive<0, 1>()
//            //            py::keep_alive<0, 2>()
//        );
        m.def(
            FactoryName.c_str(),
            [](CGP& grid,
               const SS& source_space,
               const RS& range_space,
               const MatrixTag&,
               const std::string& logging_prefix) {
              const auto cv = grid.coupling_view();
              /// which sparsity pattern for the coupling matrix??
              auto pattern = make_intersection_sparsity_pattern(
                  range_space, source_space, cv);
              return new type(cv,
                              source_space,
                              range_space,
                              new M(range_space.mapper().size(),
                                    source_space.mapper().size(),
                                    pattern),
                              logging_prefix);
            },
            "grid"_a,
            "source_space"_a,
            "range_space"_a,
            "la_backend"_a = MatrixTag(),
            "logging_prefix"_a = "",
            py::keep_alive<0, 1>()
            //            py::keep_alive<0, 2>()
        );
      }
    }
  }; // struct addbind<false, ...>

public:
  static bound_type bind_type(pybind11::module& m,
                              const std::string& matrix_id,
                              const std::string& grid_id,
                              const std::string& layer_id = "",
                              const std::string& class_id = "matrix_operator")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = XT::Common::to_camel_case(
        bindings::OperatorInterface<M, AGV, s, r, SGV, RGV>::class_name(matrix_id, grid_id, layer_id, class_id));
    bound_type c(m, ClassName.c_str(), ClassName.c_str());
    //    c.def(
    //        py::init(
    //            [](GP& grid, const SS& source_space, const RS& range_space, M& matrix, const std::string&
    //            logging_prefix) {
    //              return new type(grid.leaf_view(), source_space, range_space, matrix, logging_prefix);
    //            }),
    //        "grid"_a,
    //        "source_space"_a,
    //        "range_space"_a,
    //        "matrix"_a,
    //        "logging_prefix"_a = "",
    //        py::keep_alive<0, 1>(),
    //        py::keep_alive<1, 2>(),
    //        py::keep_alive<1, 3>(),
    //        py::keep_alive<1, 4>(),
    //        py::keep_alive<1, 5>());
    //    c.def(py::init([](GP& grid,
    //                      const SS& source_space,
    //                      const RS& range_space,
    //                      const XT::LA::SparsityPatternDefault& pattern,
    //                      const std::string& logging_prefix) {
    //            return new type(grid.leaf_view(),
    //                            source_space,
    //                            range_space,
    //                            new M(range_space.mapper().size(), source_space.mapper().size(), pattern),
    //                            logging_prefix);
    //          }),
    //          "grid"_a,
    //          "source_space"_a,
    //          "range_space"_a,
    //          "sparsity_pattern"_a,
    //          "logging_prefix"_a = "",
    //            py::keep_alive<0, 1>(),
    //          py::keep_alive<1, 2>(),
    //          py::keep_alive<1, 3>(),
    //          py::keep_alive<1, 4>());
    //    c.def(py::init([](GP& grid, const SS& source_space, const RS& range_space, const std::string& logging_prefix)
    //    {
    //            return new type(
    //                grid.leaf_view(),
    //                source_space,
    //                range_space,
    //                new M(range_space.mapper().size(),
    //                      source_space.mapper().size(),
    //                      make_element_and_intersection_sparsity_pattern(range_space, source_space,
    //                      grid.leaf_view())),
    //                logging_prefix);
    //          }),
    //          "grid"_a,
    //          "source_space"_a,
    //          "range_space"_a,
    //          "logging_prefix"_a = "",
    //            py::keep_alive<0, 1>(),
    //          py::keep_alive<1, 2>(),
    //          py::keep_alive<1, 3>(),
    //          py::keep_alive<1, 4>());

    // doing this so complicated to get an actual reference instead of a copy
    c.def_property("matrix",
                   /*getter=*/(const M& (type::*)() const) & type::matrix,
                   /*setter=*/[](type& self, const M& mat) {
                     DUNE_THROW_IF(mat.rows() != self.matrix().rows() || mat.cols() != self.matrix().cols(),
                                   XT::Common::Exceptions::shapes_do_not_match,
                                   "Cannot assign a matrix of size " << mat.rows() << "x" << mat.cols()
                                                                     << " to a matrix of size " << self.matrix().rows()
                                                                     << "x" << self.matrix().cols() << "!");
                     self.matrix() = mat;
                   });

    // methods from walker base, to allow for overloads
    //    XT::Grid::bindings::Walker<AGV>::addbind_methods(c);

    // methods from operator base, to allow for overloads
    bindings::OperatorInterface<M, AGV, s, r, SGV, RGV>::addbind_methods(c);

    // additional methods
    //    c.def("clear", [](type& self) { self.clear(); });
    c.def(
        "append",
        [](type& self, const BilinearFormType& bilinear_form, const XT::Common::Parameter& param) {
          self.append(bilinear_form, param);
        },
        "bilinear_form"_a,
        "param"_a = XT::Common::Parameter());
    c.def("__iadd__", // function ptr signature required for the right return type
          (type & (type::*)(const BilinearFormType&, const XT::Common::Parameter&)) & type::append,
          "bilinear_form"_a,
          "param"_a = XT::Common::Parameter(),
          py::is_operator());
    //    c.def("__iadd__", // function ptr signature required for the right return type
    //          (type
    //           & (type::*)(const std::tuple<const BilinearFormType&,
    //                                        const XT::Common::Parameter&,
    //                                        const XT::Grid::ElementFilter<AGV>&>&))
    //              & type::operator+=,
    //          "tuple_of_bilinearform_param_elementfilter"_a,
    //          py::is_operator());
    //    c.def(
    //        "append",
    //        [](type& self,
    //           const LocalIntersectionBilinearFormInterface<I, r, 1, F, F, s, 1, F>& local_intersection_bilinear_form,
    //           const XT::Common::Parameter& param) {
    //          self.append(local_intersection_bilinear_form, param);
    //        },
    //        "local_intersection_bilinear_form"_a,
    //        "param"_a = XT::Common::Parameter(),
    //        "intersection_filter"_a = XT::Grid::ApplyOn::AllIntersections<AGV>());
    //    c.def("__iadd__", // function ptr signature required for the right return type
    //          (type
    //           & (type::*)(const LocalIntersectionBilinearFormInterface<I, r, 1, F, F, s, 1, F>&,
    //                       const XT::Common::Parameter&,
    //                       const XT::Grid::IntersectionFilter<AGV>&))
    //              & type::append,
    //          "local_intersection_bilinear_form"_a,
    //          "param"_a = XT::Common::Parameter(),
    //          "intersection_filter"_a = XT::Grid::ApplyOn::AllIntersections<AGV>(),
    //          py::is_operator());
    //    c.def("__iadd__", // function ptr signature required for the right return type
    //          (type
    //           & (type::*)(const std::tuple<const LocalIntersectionBilinearFormInterface<I, r, 1, F, F, s, 1, F>&,
    //                                        const XT::Common::Parameter&,
    //                                        const XT::Grid::IntersectionFilter<AGV>&>&))
    //              & type::operator+=,
    //          "tuple_of_localintersectionbilinearform_param_elementfilter"_a,
    //          py::is_operator());
    //    c.def(
    //        "append",
    //        [](type& self,
    //           const LocalCouplingIntersectionBilinearFormInterface<I, r, 1, F, F, s, 1, F>&
    //               local_coupling_intersection_bilinear_form,
    //           const XT::Common::Parameter& param,
    //           const XT::Grid::IntersectionFilter<AGV>& intersection_filter) {
    //          self.append(local_coupling_intersection_bilinear_form, param, intersection_filter);
    //        },
    //        "local_coupling_intersection_bilinear_form"_a,
    //        "param"_a = XT::Common::Parameter(),
    //        "intersection_filter"_a = XT::Grid::ApplyOn::AllIntersections<AGV>());
    //    c.def("__iadd__", // function ptr signature required for the right return type
    //          (type
    //           & (type::*)(const LocalCouplingIntersectionBilinearFormInterface<I, r, 1, F, F, s, 1, F>&,
    //                       const XT::Common::Parameter&,
    //                       const XT::Grid::IntersectionFilter<AGV>&))
    //              & type::append,
    //          "local_coupling_intersection_bilinear_form"_a,
    //          "param"_a = XT::Common::Parameter(),
    //          "intersection_filter"_a = XT::Grid::ApplyOn::AllIntersections<AGV>(),
    //          py::is_operator());
    //    c.def("__iadd__", // function ptr signature required for the right return type
    //          (type
    //           & (type::*)(const std::tuple<const LocalCouplingIntersectionBilinearFormInterface<I, r, 1, F, F, s, 1,
    //           F>&,
    //                                        const XT::Common::Parameter&,
    //                                        const XT::Grid::IntersectionFilter<AGV>&>&))
    //              & type::operator+=,
    //          "tuple_of_localcouplingintersectionbilinearform_param_elementfilter"_a,
    //          py::is_operator());
    //    c.def(
    //        "assemble",
    //        [](type& self, const bool use_tbb) { self.assemble(use_tbb); },
    //        "parallel"_a = false,
    //        py::call_guard<py::gil_scoped_release>());
    return c;
  } // ... bind_type(...)

  //  static void
  //  bind_factory(pybind11::module& m, const std::string& matrix_id, const std::string& class_id = "matrix_operator")
  //  {
  //    namespace py = pybind11;
  //    using namespace pybind11::literals;

  //    //    // factories
  //    const auto FactoryName = XT::Common::to_camel_case(class_id);
  //    //    m.def(
  //    //        FactoryName.c_str(),
  //    //        [](GP& grid, const SS& source_space, const RS& range_space, M& matrix, const std::string&
  //    logging_prefix)
  //    //        {
  //    //          return new type(grid.leaf_view(), source_space, range_space, matrix, logging_prefix);
  //    //        },
  //    //        "grid"_a,
  //    //        "source_space"_a,
  //    //        "range_space"_a,
  //    //        "matrix"_a,
  //    //        "logging_prefix"_a = "",
  //    //        py::keep_alive<0, 1>(),
  //    //        py::keep_alive<0, 2>(),
  //    //        py::keep_alive<0, 3>(),
  //    //        py::keep_alive<0, 4>());
  //    addbind<>::leaf_factory(m, FactoryName);
  //    addbind<>::coupling_factory(m, FactoryName);
  //    //    m.def(
  //    //        XT::Common::to_camel_case(matrix_id + "_" + class_id).c_str(),
  //    //        [](GP& grid,
  //    //           const SS& source_space,
  //    //           const RS& range_space,
  //    //           const XT::LA::SparsityPatternDefault& pattern,
  //    //           const std::string& logging_prefix) {
  //    //          return new type(grid.leaf_view(),
  //    //                          source_space,
  //    //                          range_space,
  //    //                          new M(range_space.mapper().size(), source_space.mapper().size(), pattern),
  //    //                          logging_prefix);
  //    //        },
  //    //        "grid"_a,
  //    //        "source_space"_a,
  //    //        "range_space"_a,
  //    //        "sparsity_pattern"_a,
  //    //        "logging_prefix"_a = "",
  //    //        py::keep_alive<0, 1>(),
  //    //        py::keep_alive<0, 2>(),
  //    //        py::keep_alive<0, 3>());
  //  } // ... bind_factory(...)

  static bound_type bind(pybind11::module& m,
                         const std::string& matrix_id,
                         const std::string& grid_id,
                         const std::string& layer_id = "",
                         const std::string& class_id = "matrix_operator")
  {
    auto c = bind_type(m, matrix_id, grid_id, layer_id, class_id);
    return c;
  }

  static void
  bind_leaf_factory(pybind11::module& m, const std::string& matrix_id, const std::string& class_id = "matrix_operator")
  {
    const auto FactoryName = XT::Common::to_camel_case(class_id);
    addbind<>::leaf_factory(m, FactoryName);
  }

  static void bind_coupling_factory(pybind11::module& m,
                                    const std::string& matrix_id,
                                    const std::string& class_id = "matrix_operator")
  {
    const auto FactoryName = XT::Common::to_camel_case(class_id);
    addbind<>::coupling_factory(m, FactoryName);
  }


}; // class MatrixOperator


template <class M,
          class MatrixTag,
          class SparsityTag,
          class AGV,
          size_t s = 1,
          size_t r = s,
          class SGV = AGV,
          class RGV = AGV>
class ConstMatrixOperator
{
  using G = std::decay_t<XT::Grid::extract_grid_t<AGV>>;
  static const size_t d = G::dimension;
  using GP = XT::Grid::GridProvider<G>;

public:
  using type = GDT::ConstMatrixOperator<AGV, s, 1, r, 1, double, M, SGV, RGV>;
  using base_type = GDT::OperatorInterface<AGV, s, 1, r, 1, double, M, SGV, RGV>;
  using bound_type = pybind11::class_<type, base_type>;

private:
  using SS = typename type::SourceSpaceType;
  using RS = typename type::RangeSpaceType;

public:
  static bound_type bind_type(pybind11::module& m,
                              const std::string& matrix_id,
                              const std::string& grid_id,
                              const std::string& layer_id = "",
                              const std::string& class_id = "const_matrix_operator")
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    const auto ClassName = XT::Common::to_camel_case(
        bindings::OperatorInterface<M, AGV, s, r, SGV, RGV>::class_name(matrix_id, grid_id, layer_id, class_id));
    bound_type c(m, ClassName.c_str(), ClassName.c_str());

    return c;
  } // ... bind_type(...)


  static bound_type bind(pybind11::module& m,
                         const std::string& matrix_id,
                         const std::string& grid_id,
                         const std::string& layer_id = "",
                         const std::string& class_id = "const_matrix_operator")
  {
    auto c = bind_type(m, matrix_id, grid_id, layer_id, class_id);
    return c;
  }
}; // class ConstMatrixOperator

} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PYTHON_OPERATORS_MATRIX_BASED_HH
