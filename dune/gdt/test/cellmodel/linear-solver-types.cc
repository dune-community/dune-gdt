// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#include <string>

#include <dune/common/exceptions.hh>

#include "linear-solver-types.hh"

CellModelLinearSolverType string_to_solver_type(const std::string& type_string)
{
  if (type_string == "direct")
    return CellModelLinearSolverType::direct;
  else if (type_string == "gmres")
    return CellModelLinearSolverType::gmres;
  else if (type_string == "fgmres_gmres")
    return CellModelLinearSolverType::fgmres_gmres;
  else if (type_string == "fgmres_bicgstab")
    return CellModelLinearSolverType::fgmres_bicgstab;
  else if (type_string == "schur_gmres")
    return CellModelLinearSolverType::schur_gmres;
  else if (type_string == "schur_fgmres_gmres")
    return CellModelLinearSolverType::schur_fgmres_gmres;
  else if (type_string == "schur_fgmres_bicgstab")
    return CellModelLinearSolverType::schur_fgmres_bicgstab;
  else
    DUNE_THROW(Dune::InvalidStateException, "Solver type " + std::string(type_string) + " unknown!");
  return CellModelLinearSolverType::direct;
}

CellModelMassMatrixSolverType string_to_mass_matrix_solver_type(const std::string& type_string)
{
  if (type_string == "sparse_lu")
    return CellModelMassMatrixSolverType::sparse_lu;
  else if (type_string == "cg")
    return CellModelMassMatrixSolverType::cg;
  else if (type_string == "cg_incomplete_cholesky")
    return CellModelMassMatrixSolverType::cg_incomplete_cholesky;
  else
    DUNE_THROW(Dune::InvalidStateException, "Mass matrix solver type " + std::string(type_string) + " unknown!");
  return CellModelMassMatrixSolverType::sparse_lu;
}