# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2018)
#   René Fritze     (2018)
#   Tobias Leibner  (2018)
# ~~~

file(GLOB_RECURSE header "*.hh")

dune_pybindxi_add_module(_discretefunction_discretefunction EXCLUDE_FROM_ALL ${header} discretefunction/discretefunction.cc)
dune_pybindxi_add_module(_discretefunction_dof_vector EXCLUDE_FROM_ALL ${header} discretefunction/dof-vector.cc)
dune_pybindxi_add_module(_functionals_interfaces_common EXCLUDE_FROM_ALL ${header} functionals/interfaces_common.cc)
dune_pybindxi_add_module(_functionals_interfaces_eigen EXCLUDE_FROM_ALL ${header} functionals/interfaces_eigen.cc)
dune_pybindxi_add_module(_functionals_interfaces_istl EXCLUDE_FROM_ALL ${header} functionals/interfaces_istl.cc)
dune_pybindxi_add_module(_functionals_vector_based EXCLUDE_FROM_ALL ${header} functionals/vector-based.cc)
dune_pybindxi_add_module(_local_bilinear_forms_coupling_intersection_integrals EXCLUDE_FROM_ALL ${header} local/bilinear-forms/coupling_intersection_integrals.cc)
dune_pybindxi_add_module(_local_bilinear_forms_coupling_intersection_interface EXCLUDE_FROM_ALL ${header} local/bilinear-forms/coupling_intersection_interface.cc)
dune_pybindxi_add_module(_local_bilinear_forms_element_integrals EXCLUDE_FROM_ALL ${header} local/bilinear-forms/element_integrals.cc)
dune_pybindxi_add_module(_local_bilinear_forms_element_interface EXCLUDE_FROM_ALL ${header} local/bilinear-forms/element_interface.cc)
dune_pybindxi_add_module(_local_bilinear_forms_intersection_integrals EXCLUDE_FROM_ALL ${header} local/bilinear-forms/intersection_integrals.cc)
dune_pybindxi_add_module(_local_bilinear_forms_intersection_interface EXCLUDE_FROM_ALL ${header} local/bilinear-forms/intersection_interface.cc)
dune_pybindxi_add_module(_local_bilinear_forms_restricted_coupling_intersection_integrals EXCLUDE_FROM_ALL ${header} local/bilinear-forms/restricted_coupling_intersection_integrals.cc)
dune_pybindxi_add_module(_local_bilinear_forms_restricted_intersection_integrals EXCLUDE_FROM_ALL ${header} local/bilinear-forms/restricted_intersection_integrals.cc)
dune_pybindxi_add_module(_local_functionals_element_integrals EXCLUDE_FROM_ALL ${header} local/functionals/element_integrals.cc)
dune_pybindxi_add_module(_local_functionals_element_interface EXCLUDE_FROM_ALL ${header} local/functionals/element_interface.cc)
dune_pybindxi_add_module(_local_functionals_intersection_integrals EXCLUDE_FROM_ALL ${header} local/functionals/intersection_integrals.cc)
dune_pybindxi_add_module(_local_functionals_intersection_interface EXCLUDE_FROM_ALL ${header} local/functionals/intersection_interface.cc)
dune_pybindxi_add_module(_local_functionals_restricted_intersection_integrals EXCLUDE_FROM_ALL ${header} local/functionals/restricted_intersection_integrals.cc)
dune_pybindxi_add_module(_local_integrands_binary_element_interface EXCLUDE_FROM_ALL ${header} local/integrands/binary_element_interface.cc)
dune_pybindxi_add_module(_local_integrands_binary_intersection_interface EXCLUDE_FROM_ALL ${header} local/integrands/binary_intersection_interface.cc)
dune_pybindxi_add_module(_local_integrands_element_product EXCLUDE_FROM_ALL ${header} local/integrands/element_product.cc)
dune_pybindxi_add_module(_local_integrands_intersection_product EXCLUDE_FROM_ALL ${header} local/integrands/intersection_product.cc)
dune_pybindxi_add_module(_local_integrands_ipdg_boundary_penalty EXCLUDE_FROM_ALL ${header} local/integrands/ipdg_boundary_penalty.cc)
dune_pybindxi_add_module(_local_integrands_ipdg_inner_penalty EXCLUDE_FROM_ALL ${header} local/integrands/ipdg_inner_penalty.cc)
dune_pybindxi_add_module(_local_integrands_laplace EXCLUDE_FROM_ALL ${header} local/integrands/laplace.cc)
dune_pybindxi_add_module(_local_integrands_laplace_ipdg_dirichlet_coupling EXCLUDE_FROM_ALL ${header} local/integrands/laplace_ipdg_dirichlet_coupling.cc)
dune_pybindxi_add_module(_local_integrands_laplace_ipdg_inner_coupling EXCLUDE_FROM_ALL ${header} local/integrands/laplace_ipdg_inner_coupling.cc)
dune_pybindxi_add_module(_local_integrands_linear_advection EXCLUDE_FROM_ALL ${header} local/integrands/linear_advection.cc)
dune_pybindxi_add_module(_local_integrands_quaternary_intersection_interface EXCLUDE_FROM_ALL ${header} local/integrands/quaternary_intersection_interface.cc)
dune_pybindxi_add_module(_local_integrands_unary_element_interface EXCLUDE_FROM_ALL ${header} local/integrands/unary_element_interface.cc)
dune_pybindxi_add_module(_local_integrands_unary_intersection_interface EXCLUDE_FROM_ALL ${header} local/integrands/unary_intersection_interface.cc)
dune_pybindxi_add_module(_operators_interfaces_common EXCLUDE_FROM_ALL ${header} operators/interfaces_common.cc)
dune_pybindxi_add_module(_operators_interfaces_eigen EXCLUDE_FROM_ALL ${header} operators/interfaces_eigen.cc)
dune_pybindxi_add_module(_operators_interfaces_istl EXCLUDE_FROM_ALL ${header} operators/interfaces_istl.cc)
dune_pybindxi_add_module(_operators_matrix_based EXCLUDE_FROM_ALL ${header} operators/matrix-based.cc)
dune_pybindxi_add_module(_spaces_h1_continuous_lagrange EXCLUDE_FROM_ALL ${header} spaces/h1/continuous-lagrange.cc)
dune_pybindxi_add_module(_spaces_interface EXCLUDE_FROM_ALL ${header} spaces/interface.cc)
dune_pybindxi_add_module(_spaces_l2_discontinuous_lagrange EXCLUDE_FROM_ALL ${header} spaces/l2/discontinuous-lagrange.cc)
dune_pybindxi_add_module(_tools_dirichlet_constraints EXCLUDE_FROM_ALL ${header} tools/dirichlet-constraints.cc)
dune_pybindxi_add_module(_tools_grid_quality_estimates EXCLUDE_FROM_ALL ${header} tools/grid-quality-estimates.cc)
dune_pybindxi_add_module(_tools_sparsity_pattern EXCLUDE_FROM_ALL ${header} tools/sparsity-pattern.cc)
# dune_pybindxi_add_module(usercode EXCLUDE_FROM_ALL ${header} usercode.cc)
# dune_pybindxi_add_module(gamm_2019_talk_on_conservative_rb
#                          EXCLUDE_FROM_ALL
#                          ${header}
#                          gamm-2019-talk-on-conservative-rb.cc)

