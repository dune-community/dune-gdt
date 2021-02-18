// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef PYTHON_DUNE_GDT_OPERATORS_INTERFACES_BINDINGS_HH
#define PYTHON_DUNE_GDT_OPERATORS_INTERFACES_BINDINGS_HH

#include <dune/xt/grid/dd/glued.hh>
#include <dune/xt/grid/view/coupling.hh>

#include <python/dune/xt/grid/grids.bindings.hh>

#include "interfaces.hh"
#include "lincomb.hh"
#include "matrix-based.hh"

/**
 * We require ConstLincombOperator and LincombOperator for the numeric operators and MatrixOperator for the jacobian.
 */
template <class M, class MT, class ST, class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct OperatorInterface_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using LGV = typename G::LeafGridView;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m, const std::string& matrix_id)
  {
    using Dune::GDT::bindings::BilinearFormInterface;
    using Dune::GDT::bindings::ConstLincombOperator;
    using Dune::GDT::bindings::ConstMatrixOperator;
    using Dune::GDT::bindings::LincombOperator;
    using Dune::GDT::bindings::MatrixOperator;
    //    using Dune::GDT::bindings::ForwardOperatorInterface;
    using Dune::GDT::bindings::OperatorInterface;
    using Dune::XT::Grid::bindings::grid_name;


    OperatorInterface<M, LGV>::bind(m, matrix_id, grid_name<G>::value(), "leaf");
    //    BilinearFormInterface<GV>::bind(m, matrix_id, grid_name<G>::value(), "leaf");
    //    ForwardOperatorInterface<GV>::bind(m, matrix_id, grid_name<G>::value(), "leaf");
    //    ConstLincombOperator<M, GV>::bind(m, matrix_id, grid_name<G>::value(), "leaf");
    //    LincombOperator<M, GV>::bind(m, matrix_id, grid_name<G>::value(), "leaf");
    ConstMatrixOperator<M, MT, ST, LGV>::bind_type(m, matrix_id, grid_name<G>::value(), "leaf");
    MatrixOperator<M, MT, ST, LGV>::bind_type(m, matrix_id, grid_name<G>::value(), "leaf");
    if (d > 1) {
      OperatorInterface<M, LGV, d, 1>::bind(m, matrix_id, grid_name<G>::value(), "leaf");
      //      BilinearFormInterface<GV, d, 1>::bind(m, matrix_id, grid_name<G>::value(), "leaf");
      //      ForwardOperatorInterface<GV, d, 1>::bind(m, matrix_id, grid_name<G>::value(), "leaf");
      //      ConstLincombOperator<M, GV, d, 1>::bind(m, matrix_id, grid_name<G>::value(), "leaf");
      //      LincombOperator<M, GV, d, 1>::bind(m, matrix_id, grid_name<G>::value(), "leaf");
      ConstMatrixOperator<M, MT, ST, LGV, d, 1>::bind_type(m, matrix_id, grid_name<G>::value(), "leaf");
      MatrixOperator<M, MT, ST, LGV, d, 1>::bind_type(m, matrix_id, grid_name<G>::value(), "leaf");

      OperatorInterface<M, LGV, 1, d>::bind(m, matrix_id, grid_name<G>::value(), "leaf");
      //      BilinearFormInterface<GV, 1, d>::bind(m, matrix_id, grid_name<G>::value(), "leaf");
      //      ForwardOperatorInterface<GV, 1, d>::bind(m, matrix_id, grid_name<G>::value(), "leaf");
      //      ConstLincombOperator<M, GV, 1, d>::bind(m, matrix_id, grid_name<G>::value(), "leaf");
      //      LincombOperator<M, GV, 1, d>::bind(m, matrix_id, grid_name<G>::value(), "leaf");
      ConstMatrixOperator<M, MT, ST, LGV, 1, d>::bind_type(m, matrix_id, grid_name<G>::value(), "leaf");
      MatrixOperator<M, MT, ST, LGV, 1, d>::bind_type(m, matrix_id, grid_name<G>::value(), "leaf");

      OperatorInterface<M, LGV, d, d>::bind(m, matrix_id, grid_name<G>::value(), "leaf");
      //      BilinearFormInterface<GV, d, d>::bind(m, matrix_id, grid_name<G>::value(), "leaf");
      //      ForwardOperatorInterface<GV, d, d>::bind(m, matrix_id, grid_name<G>::value(), "leaf");
      //      ConstLincombOperator<M, GV, d, d>::bind(m, matrix_id, grid_name<G>::value(), "leaf");
      //      LincombOperator<M, GV, d, d>::bind(m, matrix_id, grid_name<G>::value(), "leaf");
      ConstMatrixOperator<M, MT, ST, LGV, d, d>::bind_type(m, matrix_id, grid_name<G>::value(), "leaf");
      MatrixOperator<M, MT, ST, LGV, d, d>::bind_type(m, matrix_id, grid_name<G>::value(), "leaf");
    }
    // add your extra dimensions here
    // ...

#if HAVE_DUNE_GRID_GLUE
    if constexpr (d == 2) {
      using GridGlueType = Dune::XT::Grid::DD::Glued<G, G, Dune::XT::Grid::Layers::leaf>;
      using CGV = Dune::XT::Grid::CouplingGridView<GridGlueType>;
        OperatorInterface<M, CGV, 1, 1, LGV, LGV>::bind(m, matrix_id, grid_name<G>::value(), "coupling");
        //    BilinearFormInterface<GV>::bind(m, matrix_id, grid_name<G>::value(), "coupling");
        //    ForwardOperatorInterface<GV>::bind(m, matrix_id, grid_name<G>::value(), "coupling");
        //    ConstLincombOperator<M, GV>::bind(m, matrix_id, grid_name<G>::value(), "coupling");
        //    LincombOperator<M, GV>::bind(m, matrix_id, grid_name<G>::value(), "coupling");
        ConstMatrixOperator<M, MT, ST, CGV, 1, 1, LGV, LGV>::bind_type(m, matrix_id, grid_name<G>::value(), "coupling");
        MatrixOperator<M, MT, ST, CGV, 1, 1, LGV, LGV>::bind_type(m, matrix_id, grid_name<G>::value(), "coupling");
        if (d > 1) {
          OperatorInterface<M, CGV, d, 1, LGV, LGV>::bind(m, matrix_id, grid_name<G>::value(), "coupling");
          //      BilinearFormInterface<GV, d, 1>::bind(m, matrix_id, grid_name<G>::value(), "coupling");
          //      ForwardOperatorInterface<GV, d, 1>::bind(m, matrix_id, grid_name<G>::value(), "coupling");
          //      ConstLincombOperator<M, GV, d, 1>::bind(m, matrix_id, grid_name<G>::value(), "coupling");
          //      LincombOperator<M, GV, d, 1>::bind(m, matrix_id, grid_name<G>::value(), "coupling");
          ConstMatrixOperator<M, MT, ST, CGV, d, 1, LGV, LGV>::bind_type(m, matrix_id, grid_name<G>::value(), "coupling");
          MatrixOperator<M, MT, ST, CGV, d, 1, LGV, LGV>::bind_type(m, matrix_id, grid_name<G>::value(), "coupling");

          OperatorInterface<M, CGV, 1, d, LGV, LGV>::bind(m, matrix_id, grid_name<G>::value(), "coupling");
          //      BilinearFormInterface<GV, 1, d>::bind(m, matrix_id, grid_name<G>::value(), "coupling");
          //      ForwardOperatorInterface<GV, 1, d>::bind(m, matrix_id, grid_name<G>::value(), "coupling");
          //      ConstLincombOperator<M, GV, 1, d>::bind(m, matrix_id, grid_name<G>::value(), "coupling");
          //      LincombOperator<M, GV, 1, d>::bind(m, matrix_id, grid_name<G>::value(), "coupling");
          ConstMatrixOperator<M, MT, ST, CGV, 1, d, LGV, LGV>::bind_type(m, matrix_id, grid_name<G>::value(), "coupling");
          MatrixOperator<M, MT, ST, CGV, 1, d, LGV, LGV>::bind_type(m, matrix_id, grid_name<G>::value(), "coupling");

          OperatorInterface<M, CGV, d, d, LGV, LGV>::bind(m, matrix_id, grid_name<G>::value(), "coupling");
          //      BilinearFormInterface<GV, d, d>::bind(m, matrix_id, grid_name<G>::value());
          //      ForwardOperatorInterface<GV, d, d>::bind(m, matrix_id, grid_name<G>::value());
          //      ConstLincombOperator<M, GV, d, d>::bind(m, matrix_id, grid_name<G>::value());
          //      LincombOperator<M, GV, d, d>::bind(m, matrix_id, grid_name<G>::value());
          ConstMatrixOperator<M, MT, ST, CGV, d, d, LGV, LGV>::bind_type(m, matrix_id, grid_name<G>::value(), "coupling");
          MatrixOperator<M, MT, ST, CGV, d, d, LGV, LGV>::bind_type(m, matrix_id, grid_name<G>::value(), "coupling");
        }
        // add your extra dimensions here
        // ...
    }
#endif

    OperatorInterface_for_all_grids<M, MT, ST, Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m, matrix_id);
  }
};

template <class M, class MT, class ST>
struct OperatorInterface_for_all_grids<M, MT, ST, Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/, const std::string& /*matrix_id*/) {}
};


#endif // PYTHON_DUNE_GDT_OPERATORS_INTERFACES_BINDINGS_HH
