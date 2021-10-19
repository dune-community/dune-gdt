// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014, 2016)
//   Rene Milk       (2014)
//   Tobias Leibner  (2016)

#include "config.h"

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/gdt/test/cellmodel/cellmodel.hh>

PYBIND11_MODULE(cellmodel, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;
  using namespace Dune;
  py::enum_<CellModelLinearSolverType>(m, "CellModelLinearSolverType")
      .value("direct", CellModelLinearSolverType::direct)
      .value("gmres", CellModelLinearSolverType::gmres)
      .value("fgmres_gmres", CellModelLinearSolverType::fgmres_gmres)
      .value("fgmres_bicgstab", CellModelLinearSolverType::fgmres_bicgstab)
      .value("schur_gmres", CellModelLinearSolverType::schur_gmres)
      .value("schur_fgmres_gmres", CellModelLinearSolverType::schur_fgmres_gmres)
      .value("schur_fgmres_bicgstab", CellModelLinearSolverType::schur_fgmres_bicgstab);

  py::enum_<StokesSolverType>(m, "StokesSolverType")
      .value("direct", StokesSolverType::direct)
      .value("schur_cg_A_direct", StokesSolverType::schur_cg_A_direct)
      .value("schur_cg_A_direct_prec_mass", StokesSolverType::schur_cg_A_direct_prec_mass)
      .value("schur_cg_A_direct_prec_masslumped", StokesSolverType::schur_cg_A_direct_prec_masslumped);

  py::enum_<CellModelMassMatrixSolverType>(m, "CellModelMassMatrixSolverType")
      .value("sparse_ldlt", CellModelMassMatrixSolverType::sparse_ldlt)
      .value("sparse_lu", CellModelMassMatrixSolverType::sparse_lu)
      .value("cg", CellModelMassMatrixSolverType::cg)
      .value("cg_incomplete_cholesky", CellModelMassMatrixSolverType::cg_incomplete_cholesky);

  py::class_<CellModelSolver>(m, "CellModelSolver")
      .def(py::init<const std::string /*testcase*/,
                    const double /*t_end*/,
                    const double /*dt*/,
                    const unsigned int /*num_elements_x*/,
                    const unsigned int /*num_elements_y*/,
                    const int /*pol_order*/,
                    const bool /*use_tbb*/,
                    const double /*Be*/,
                    const double /*Ca*/,
                    const double /*Pa*/,
                    const double /*Re*/,
                    const double /*Fa*/,
                    const double /*xi*/,
                    const double /*kappa*/,
                    const double /*c_1*/,
                    const double /*gamma*/,
                    const double /*epsilon*/,
                    const int /*overintegrate*/,
                    const CellModelLinearSolverType /*pfield_solver_type*/,
                    const CellModelMassMatrixSolverType /*pfield_mass_matrix_solver_type*/,
                    const CellModelLinearSolverType /*ofield_solver_type*/,
                    const CellModelMassMatrixSolverType /*ofield_mass_matrix_solver_type*/,
                    const StokesSolverType /*stokes_solver_type*/,
                    const double /*outer_reduction*/,
                    const int /*outer_restart*/,
                    const int /*outer_verbose*/,
                    const double /*inner_reduction*/,
                    const int /*inner_maxit*/,
                    const int /*inner_verbose*/,
                    const double /*bending*/,
                    const double /*conservative*/>(),
           "testcase"_a = "single_cell",
           "t_end"_a = 1.,
           "dt"_a = 0.01,
           "num_elements_x"_a = 50,
           "num_elements_y"_a = 50,
           "pol_order"_a = 1,
           "use_tbb"_a = true,
           "Be"_a = 0.3, // bending capillary number, ratio of viscous forces to bending forces
           "Ca"_a = 0.1, // capillary number, ratio of viscous forces to surface tension forces
           "Pa"_a = 1, // polarization elasticity number
           "Re"_a = 5e-13, // Reynolds number
           "Fa"_a = 1., // active force number
           "xi"_a = 1.1, // alignment of P with the flow, > 0 for rod-like cells and < 0 for oblate ones
           "kappa"_a = 1.65, // eta_rot/eta, scaling factor between rotational and dynamic viscosity
           "c_1"_a = 5., // double well shape parameter
           "gamma"_a = 0.025, // phase field mobility coefficient
           "epsilon"_a = 0.21, // phase field parameter
           "overintegrate"_a = 2,
           "pfield_solver_type"_a = CellModelLinearSolverType::gmres,
           "pfield_mass_matrix_solver_type"_a = CellModelMassMatrixSolverType::sparse_ldlt,
           "ofield_solver_type"_a = CellModelLinearSolverType::schur_gmres,
           "ofield_mass_matrix_solver_type"_a = CellModelMassMatrixSolverType::sparse_ldlt,
           "stokes_solver_type"_a = StokesSolverType::schur_cg_A_direct_prec_mass,
           "outer_reduction"_a = 1e-14,
           "outer_restart"_a = 100,
           "outer_verbose"_a = 0,
           "inner_reduction"_a = 1e-3,
           "inner_maxit"_a = 10,
           "inner_verbose"_a = 0,
           "bending"_a = true,
           "conservative"_a = true)
      .def("visualize",
           &CellModelSolver::visualize,
           "prefix"_a,
           "step"_a,
           "t"_a,
           "subsampling"_a = true,
           "vtu"_a = true,
           "txt"_a = true,
           "timings"_a = true)
      .def("visualize_pfield", &CellModelSolver::visualize_pfield)
      .def("visualize_ofield", &CellModelSolver::visualize_ofield)
      .def("visualize_stokes", &CellModelSolver::visualize_stokes)
      .def("prepare_pfield_operator", &CellModelSolver::prepare_pfield_operator, "cell"_a, "restricted"_a = false)
      .def("prepare_ofield_operator", &CellModelSolver::prepare_ofield_operator, "cell"_a, "restricted"_a = false)
      .def("prepare_stokes_operator", &CellModelSolver::prepare_stokes_operator, "restricted"_a = false)
      .def("apply_pfield_operator", &CellModelSolver::apply_pfield_operator, "y"_a, "cell"_a, "restricted"_a = false)
      .def("apply_ofield_operator", &CellModelSolver::apply_ofield_operator, "y"_a, "cell"_a, "restricted"_a = false)
      .def("apply_stokes_operator", &CellModelSolver::apply_stokes_operator, "y"_a, "restricted"_a = false)
      .def("set_pfield_vec", &CellModelSolver::set_pfield_vec)
      .def("set_ofield_vec", &CellModelSolver::set_ofield_vec)
      .def("set_stokes_vec", &CellModelSolver::set_stokes_vec)
      .def("set_pfield_vec_dofs",
           py::overload_cast<const size_t, const pybind11::array_t<double>&, const pybind11::list&>(
               &CellModelSolver::set_pfield_vec_dofs))
      .def("set_ofield_vec_dofs",
           py::overload_cast<const size_t, const pybind11::array_t<double>&, const pybind11::list&>(
               &CellModelSolver::set_ofield_vec_dofs))
      .def("set_stokes_vec_dofs",
           py::overload_cast<const pybind11::array_t<double>&, const pybind11::list&>(
               &CellModelSolver::set_stokes_vec_dofs))
      .def("solve", &CellModelSolver::solve)
      .def("next_n_timesteps", &CellModelSolver::next_n_timesteps)
      .def("apply_inverse_pfield_operator", &CellModelSolver::apply_inverse_pfield_operator)
      .def("apply_inverse_ofield_operator", &CellModelSolver::apply_inverse_ofield_operator)
      .def("apply_inverse_stokes_operator", &CellModelSolver::apply_inverse_stokes_operator)
      .def("apply_pfield_product_operator", &CellModelSolver::apply_pfield_product_operator)
      .def("apply_ofield_product_operator", &CellModelSolver::apply_ofield_product_operator)
      .def("apply_stokes_product_operator", &CellModelSolver::apply_stokes_product_operator)
      .def("set_pfield_jacobian_state",
           &CellModelSolver::set_pfield_jacobian_state,
           "source"_a,
           "cell"_a,
           "restricted"_a = false)
      .def("set_ofield_jacobian_state",
           &CellModelSolver::set_ofield_jacobian_state,
           "source"_a,
           "cell"_a,
           "restricted"_a = false)
      .def("set_pfield_jacobian_state_dofs", &CellModelSolver::set_pfield_jacobian_state_dofs)
      .def("set_ofield_jacobian_state_dofs", &CellModelSolver::set_ofield_jacobian_state_dofs)
      .def("apply_pfield_jacobian",
           &CellModelSolver::apply_pfield_jacobian,
           "source"_a,
           "cell"_a,
           "restricted"_a = false)
      .def("apply_ofield_jacobian",
           &CellModelSolver::apply_ofield_jacobian,
           "source"_a,
           "cell"_a,
           "restricted"_a = false)
      .def("apply_stokes_jacobian", &CellModelSolver::apply_stokes_jacobian, "source"_a, "restricted"_a = false)
      .def("apply_inverse_pfield_jacobian", &CellModelSolver::apply_inverse_pfield_jacobian)
      .def("apply_inverse_ofield_jacobian", &CellModelSolver::apply_inverse_ofield_jacobian)
      .def("apply_inverse_stokes_jacobian", &CellModelSolver::apply_inverse_stokes_jacobian)
      .def("update_pfield_parameters",
           &CellModelSolver::update_pfield_parameters,
           "Be"_a,
           "Ca"_a,
           "Pa"_a,
           "cell"_a,
           "restricted"_a = false)
      .def("update_ofield_parameters",
           &CellModelSolver::update_ofield_parameters,
           "Pa"_a,
           "cell"_a,
           "restricted"_a = false)
      .def("num_cells", &CellModelSolver::num_cells)
      .def("finished", &CellModelSolver::finished)
      .def("linear", &CellModelSolver::linear)
      .def("reset", &CellModelSolver::reset)
      .def("pfield_vec", &CellModelSolver::pfield_vec, py::return_value_policy::copy)
      .def("ofield_vec", &CellModelSolver::ofield_vec, py::return_value_policy::copy)
      .def("stokes_vec", &CellModelSolver::stokes_vec, py::return_value_policy::copy)
      .def("compute_restricted_pfield_dofs", &CellModelSolver::compute_restricted_pfield_dofs)
      .def("compute_restricted_ofield_dofs", &CellModelSolver::compute_restricted_ofield_dofs)
      .def("compute_restricted_stokes_dofs", &CellModelSolver::compute_restricted_stokes_dofs)
      .def("pfield_deim_source_dofs", &CellModelSolver::pfield_deim_source_dofs, py::return_value_policy::copy)
      .def("ofield_deim_source_dofs", &CellModelSolver::ofield_deim_source_dofs, py::return_value_policy::copy)
      .def("stokes_deim_source_dofs", &CellModelSolver::stokes_deim_source_dofs, py::return_value_policy::copy);
}
