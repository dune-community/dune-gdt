// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#include <string>

#include "config.h"

#include <dune/common/exceptions.hh>

#include "linear-solver-types.hh"

CellModelLinearSolverType string_to_solver_type(const std::string& type_string)
{
  if (type_string == "direct")
    return CellModelLinearSolverType::direct;
  if (type_string == "gmres")
    return CellModelLinearSolverType::gmres;
  if (type_string == "fgmres_gmres")
    return CellModelLinearSolverType::fgmres_gmres;
  if (type_string == "fgmres_bicgstab")
    return CellModelLinearSolverType::fgmres_bicgstab;
  if (type_string == "schur_gmres")
    return CellModelLinearSolverType::schur_gmres;
  if (type_string == "schur_fgmres_gmres")
    return CellModelLinearSolverType::schur_fgmres_gmres;
  if (type_string == "schur_fgmres_bicgstab")
    return CellModelLinearSolverType::schur_fgmres_bicgstab;
  DUNE_THROW(Dune::InvalidStateException, "Solver type " + std::string(type_string) + " unknown!");
  return CellModelLinearSolverType::direct;
}

StokesSolverType string_to_stokes_solver_type(const std::string& type_string)
{
  if (type_string == "direct")
    return StokesSolverType::direct;
  if (type_string == "schur_cg_A_direct")
    return StokesSolverType::schur_cg_A_direct;
  if (type_string == "schur_cg_A_direct_prec_mass")
    return StokesSolverType::schur_cg_A_direct_prec_mass;
  if (type_string == "schur_cg_A_direct_prec_masslumped")
    return StokesSolverType::schur_cg_A_direct_prec_masslumped;
  DUNE_THROW(Dune::InvalidStateException, "Stokes solver type " + std::string(type_string) + " unknown!");
  return StokesSolverType::direct;
}


CellModelMassMatrixSolverType string_to_mass_matrix_solver_type(const std::string& type_string)
{
  if (type_string == "sparse_ldlt")
    return CellModelMassMatrixSolverType::sparse_ldlt;
  if (type_string == "sparse_lu")
    return CellModelMassMatrixSolverType::sparse_lu;
  if (type_string == "cg")
    return CellModelMassMatrixSolverType::cg;
  if (type_string == "cg_incomplete_cholesky")
    return CellModelMassMatrixSolverType::cg_incomplete_cholesky;
  DUNE_THROW(Dune::InvalidStateException, "Mass matrix solver type " + std::string(type_string) + " unknown!");
  return CellModelMassMatrixSolverType::sparse_lu;
}
