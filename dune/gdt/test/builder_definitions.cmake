# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2017)
#   Rene Milk       (2016 - 2018)
# ~~~

set(DXT_BIN_COUNT "13" CACHE STRING "number of bins for test targets")
add_custom_target(test_binaries_builder_0 DEPENDS test_mpi_linearelliptic__swipdg_discretization)
set_tests_properties(test_mpi_linearelliptic__swipdg_discretization PROPERTIES LABELS "builder_0")
add_custom_target(test_binaries_builder_1
                  DEPENDS
                    headercheck__dune_gdt_assembler_wrapper.hh
                    headercheck__dune_gdt_local_fluxes_godunov.hh
                    headercheck__dune_gdt_local_fluxes_musta.hh
                    headercheck__dune_gdt_local_functionals_interfaces.hh
                    headercheck__dune_gdt_local_operators_fv.hh
                    headercheck__dune_gdt_operators_fv.hh
                    headercheck__dune_gdt_operators_fv_laxwendroff.hh
                    headercheck__dune_gdt_operators_fv_slopelimiters.hh
                    headercheck__dune_gdt_prolongations_l2.hh
                    headercheck__dune_gdt_spaces_dg_default.hh
                    headercheck__dune_gdt_test_hyperbolic_eocexpectations-fv-boltzmanncheckerboard-2dyaspgrid.hh
                    headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_fokkerplanck_fokkerplanckequation.hh
                    headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_fokkerplanck_planesource.hh
                    headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_kinetictransport_checkerboard.hh
                    headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_kinetictransport_twobeams.hh
                    headercheck__dune_gdt_test_instationary-testcase.hh
                    headercheck__dune_gdt_test_linearelliptic_discretizers_ipdg.hh
                    headercheck__dune_gdt_test_linearelliptic_eocexpectations_all.hh
                    headercheck__dune_gdt_test_linearelliptic_eocexpectations_block-swipdg-esv2007-2dyaspgrid.hh
                    headercheck__dune_gdt_test_linearelliptic_eocexpectations_cg-esv2007-2dyaspgrid.hh
                    headercheck__dune_gdt_test_linearelliptic_swipdg-estimators.hh
                    headercheck__python_dune_gdt_spaces_constraints.hh
                    test_linearelliptic__swipdg_discretization
                    test_operators__darcy)
set_tests_properties(test_linearelliptic__swipdg_discretization PROPERTIES LABELS "builder_1")
set_tests_properties(test_operators__darcy PROPERTIES LABELS "builder_1")
add_custom_target(test_binaries_builder_2
                  DEPENDS headercheck__dune_gdt_assembler_system.lib.hh
                          headercheck__dune_gdt_local_operators_dirichlet-projection.hh
                          headercheck__dune_gdt_operators_elliptic.hh
                          headercheck__dune_gdt_operators_interfaces.hh
                          headercheck__dune_gdt_projections_l2.hh
                          headercheck__dune_gdt_spaces.hh
                          headercheck__dune_gdt_spaces_constraints.hh
                          headercheck__dune_gdt_spaces_dg.hh
                          headercheck__dune_gdt_spaces_parallel_helper.hh
                          headercheck__dune_gdt_test_hyperbolic_eocexpectations-fv-burgers-2dyaspgrid.hh
                          headercheck__dune_gdt_test_hyperbolic_eocexpectations-fv-sourcebeam-1dyaspgrid.hh
                          headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_basisfunctions_hatfunctions.hh
                          headercheck__dune_gdt_test_hyperbolic_problems_sodshocktube.hh
                          headercheck__dune_gdt_test_linearelliptic_eocexpectations_cg-er2007-2dyaspgrid.hh
                          headercheck__dune_gdt_test_linearelliptic_eocexpectations_cg-mixedboundary-2dalugrid.hh
                          headercheck__dune_gdt_test_linearelliptic_problems_AO2013.hh
                          headercheck__dune_gdt_test_linearelliptic_problems_ER2007.hh
                          headercheck__dune_gdt_test_operators_l2.hh
                          headercheck__python_dune_gdt_operators_weighted-l2.hh
                          headercheck__python_dune_gdt_playground_operators_OS2015.hh
                          headercheck__python_dune_gdt_shared.hh
                          headercheck__python_dune_gdt_spaces_fv.hh
                          test_mpi_linearelliptic__block_swipdg_discretization
                          test_mpi_operators__l2__operator)
set_tests_properties(test_mpi_linearelliptic__block_swipdg_discretization PROPERTIES LABELS "builder_2")
set_tests_properties(test_mpi_operators__l2__operator PROPERTIES LABELS "builder_2")
add_custom_target(
  test_binaries_builder_3
  DEPENDS headercheck__dune_gdt_local_assembler.hh
          headercheck__dune_gdt_local_discretefunction.hh
          headercheck__dune_gdt_local_fluxes_laxfriedrichs.hh
          headercheck__dune_gdt_operators_base.hh
          headercheck__dune_gdt_operators_fluxreconstruction.hh
          headercheck__dune_gdt_operators_fv_datahandle.hh
          headercheck__dune_gdt_operators_fv_reconstructed_function.hh
          headercheck__dune_gdt_prolongations_lagrange.hh
          headercheck__dune_gdt_spaces_basefunctionset_fv.hh
          headercheck__dune_gdt_spaces_basefunctionset_interface.hh
          headercheck__dune_gdt_spaces_datahandles.hh
          headercheck__dune_gdt_spaces_interface.hh
          headercheck__dune_gdt_test_hyperbolic_eocexpectations-fv-transport-1dyaspgrid.hh
          headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_kinetictransport_kinetictransportequation.hh
          headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_kinetictransport_pointsource.hh
          headercheck__dune_gdt_test_linearelliptic_eocexpectations.hh
          headercheck__dune_gdt_test_linearelliptic_eocexpectations_base.hh
          headercheck__dune_gdt_test_linearelliptic_eocexpectations_cg-esv2007-2dalugrid.hh
          headercheck__dune_gdt_test_linearelliptic_problems.hh
          headercheck__dune_gdt_test_linearelliptic_problems_mixedboundary.hh
          headercheck__dune_gdt_test_operators_laplace.hh
          headercheck__dune_gdt_test_operators_weighted-l2.hh
          headercheck__dune_gdt_timestepper_matrix-exponential.hh
          headercheck__python_dune_gdt_functionals_elliptic-ipdg_bindings.hh
          headercheck__python_dune_gdt_spaces_dg.hh
          test_mpi_operators__l2__localizable_product
          test_mpi_prolongations_cg_dg)
set_tests_properties(test_mpi_operators__l2__localizable_product PROPERTIES LABELS "builder_3")
set_tests_properties(test_mpi_prolongations_cg_dg PROPERTIES LABELS "builder_3")
add_custom_target(test_binaries_builder_4
                  DEPENDS
                    headercheck__dune_gdt_discretefunction_reinterpret.hh
                    headercheck__dune_gdt_discretizations_default_fv-internal.hh
                    headercheck__dune_gdt_local_fluxes_interfaces.hh
                    headercheck__dune_gdt_local_integrands_elliptic.hh
                    headercheck__dune_gdt_operators_elliptic-ipdg.hh
                    headercheck__dune_gdt_spaces_rt_interface.hh
                    headercheck__dune_gdt_test_hyperbolic_fv-discretization.hh
                    headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_basisfunctions_piecewise_monomials.hh
                    headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_fokkerplanck_sourcebeam.hh
                    headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_fokkerplanck_twobeams.hh
                    headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_kineticequation.hh
                    headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_kinetictransport_sourcebeam.hh
                    headercheck__dune_gdt_test_hyperbolic_problems_shallowwater.hh
                    headercheck__dune_gdt_test_linearelliptic_eocexpectations_cg-ao2013-2dyaspgrid.hh
                    headercheck__dune_gdt_test_linearelliptic_eocexpectations_cg-mixedboundary-2dyaspgrid.hh
                    headercheck__dune_gdt_test_linearelliptic_eocexpectations_swipdg-esv2007-2dalugrid.hh
                    headercheck__dune_gdt_test_projections_l2-local.hh
                    headercheck__dune_gdt_timestepper_factory.hh
                    headercheck__python_dune_gdt_assembler_system.hh
                    headercheck__python_dune_gdt_operators_oswaldinterpolation.hh
                    headercheck__python_dune_gdt_spaces_cg.hh
                    test_linearelliptic__swipdg_estimators
                    test_mpi_operators__weighted_l2__localizable_product)
set_tests_properties(test_linearelliptic__swipdg_estimators PROPERTIES LABELS "builder_4")
set_tests_properties(test_mpi_operators__weighted_l2__localizable_product PROPERTIES LABELS "builder_4")
add_custom_target(test_binaries_builder_5
                  DEPENDS headercheck__dune_gdt_discretefunction_default.hh
                          headercheck__dune_gdt_local_integrands_sipdg.hh
                          headercheck__dune_gdt_operators_fv_force.hh
                          headercheck__dune_gdt_operators_fv_kinetic.hh
                          headercheck__dune_gdt_operators_fv_musta.hh
                          headercheck__dune_gdt_spaces_basefunctionset_default.hh
                          headercheck__dune_gdt_spaces_cg.lib.hh
                          headercheck__dune_gdt_spaces_cg_default.hh
                          headercheck__dune_gdt_spaces_mapper_interfaces.hh
                          headercheck__dune_gdt_test_hyperbolic_eocexpectations-fv-shocktube-1dyaspgrid.hh
                          headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_basisfunctions_base.hh
                          headercheck__dune_gdt_test_instationary-eocstudy.hh
                          headercheck__dune_gdt_test_linearelliptic_discretizers_base.hh
                          headercheck__dune_gdt_test_linearelliptic_eocexpectations_swipdg-er2007-2dyaspgrid.hh
                          headercheck__dune_gdt_test_linearelliptic_eocexpectations_swipdg-mixedboundary-2dalugrid.hh
                          headercheck__dune_gdt_test_linearelliptic_problems_ESV2007.hh
                          headercheck__dune_gdt_test_projections_base.hh
                          headercheck__dune_gdt_test_projections_lagrange.hh
                          headercheck__dune_gdt_test_prolongations_base.hh
                          headercheck__dune_gdt_timestepper_enums.hh
                          headercheck__dune_gdt_timestepper_fractional-step.hh
                          headercheck__python_dune_gdt_operators_base.hh
                          headercheck__python_dune_gdt_operators_fluxreconstruction.hh
                          headercheck__python_dune_gdt_spaces_bindings.hh
                          test_linearelliptic__cg_discretization
                          test_mpi_operators__elliptic__matrix_operator
                          test_mpi_operators__laplace__localizable_product)
set_tests_properties(test_linearelliptic__cg_discretization PROPERTIES LABELS "builder_5")
set_tests_properties(test_mpi_operators__elliptic__matrix_operator PROPERTIES LABELS "builder_5")
set_tests_properties(test_mpi_operators__laplace__localizable_product PROPERTIES LABELS "builder_5")
add_custom_target(test_binaries_builder_6
                  DEPENDS headercheck__dune_gdt_local_dof-vector.hh
                          headercheck__dune_gdt_local_fluxes_force.hh
                          headercheck__dune_gdt_local_fluxes_kinetic.hh
                          headercheck__dune_gdt_local_operators_integrals.hh
                          headercheck__dune_gdt_local_operators_lambda.hh
                          headercheck__dune_gdt_operators_fv_laxfriedrichs.hh
                          headercheck__dune_gdt_operators_fv_rhs.hh
                          headercheck__dune_gdt_prolongations_l2-global.hh
                          headercheck__dune_gdt_spaces_rt.hh
                          headercheck__dune_gdt_test_grids.hh
                          headercheck__dune_gdt_test_hyperbolic_eocexpectations-fv-shallowwater-1dyaspgrid.hh
                          headercheck__dune_gdt_test_hyperbolic_eocexpectations.hh
                          headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_basisfunctions.hh
                          headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_fokkerplanck_twopulses.hh
                          headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_kinetictransport_planesource.hh
                          headercheck__dune_gdt_test_linearelliptic_discretizers_cg.hh
                          headercheck__dune_gdt_test_linearelliptic_eocexpectations_swipdg-er2007-2dalugrid.hh
                          headercheck__dune_gdt_test_linearelliptic_eocstudy.hh
                          headercheck__dune_gdt_test_prolongations_l2.hh
                          headercheck__dune_gdt_test_spaces_rt.hh
                          headercheck__python_dune_gdt_discretefunction_bindings.hh
                          headercheck__python_dune_gdt_operators_l2.hh
                          test_mpi_operators__weighted_l2__operator
                          test_mpi_projections_part_1_a
                          test_spaces_dg)
set_tests_properties(test_mpi_operators__weighted_l2__operator PROPERTIES LABELS "builder_6")
set_tests_properties(test_mpi_projections_part_1_a PROPERTIES LABELS "builder_6")
set_tests_properties(test_spaces_dg PROPERTIES LABELS "builder_6")
add_custom_target(test_binaries_builder_7
                  DEPENDS
                    headercheck__dune_gdt_local_operators_interfaces.hh
                    headercheck__dune_gdt_local_operators_l2-projection.hh
                    headercheck__dune_gdt_projections_dirichlet.hh
                    headercheck__dune_gdt_projections_l2-local.hh
                    headercheck__dune_gdt_spaces_fv_product.hh
                    headercheck__dune_gdt_spaces_mapper_fv.hh
                    headercheck__dune_gdt_spaces_mapper_product.hh
                    headercheck__dune_gdt_test_hyperbolic_eocexpectations-fv-transport-2dyaspgrid.hh
                    headercheck__dune_gdt_test_hyperbolic_eocstudy.hh
                    headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_basisfunctions_spherical_harmonics.hh
                    headercheck__dune_gdt_test_linearelliptic_discretizers_block-ipdg.hh
                    headercheck__dune_gdt_test_linearelliptic_eocexpectations_cg-ao2013-2dalugrid.hh
                    headercheck__dune_gdt_test_linearelliptic_eocexpectations_cg-spe10-2dyaspgrid.hh
                    headercheck__dune_gdt_test_linearelliptic_problems_base.hh
                    headercheck__dune_gdt_test_operators_darcy.hh
                    headercheck__dune_gdt_test_prolongations_l2-global.hh
                    headercheck__dune_gdt_timestepper_adaptive-rungekutta.hh
                    headercheck__dune_gdt_timestepper_explicit-rungekutta.hh
                    headercheck__python_dune_gdt_local_elliptic-ipdg-operators.hh
                    headercheck__python_dune_gdt_operators_elliptic-ipdg_bindings.hh
                    headercheck__python_dune_gdt_playground_spaces_block.hh
                    test_empty
                    test_mpi_operators__elliptic__operator
                    test_mpi_projections_part_2_a
                    test_spaces_rt)
set_tests_properties(test_empty PROPERTIES LABELS "builder_7")
set_tests_properties(test_mpi_operators__elliptic__operator PROPERTIES LABELS "builder_7")
set_tests_properties(test_mpi_projections_part_2_a PROPERTIES LABELS "builder_7")
set_tests_properties(test_spaces_rt PROPERTIES LABELS "builder_7")
add_custom_target(test_binaries_builder_8
                  DEPENDS headercheck__dune_gdt_assembler_local-accumulators.hh
                          headercheck__dune_gdt_discretefunction_datahandle.hh
                          headercheck__dune_gdt_discretizations_default.hh
                          headercheck__dune_gdt_discretizations_interfaces.hh
                          headercheck__dune_gdt_exceptions.hh
                          headercheck__dune_gdt_functionals_base.hh
                          headercheck__dune_gdt_functionals_l2.hh
                          headercheck__dune_gdt_local_integrands_lambda.hh
                          headercheck__dune_gdt_operators_fv_base.hh
                          headercheck__dune_gdt_operators_fv_reconstruction.hh
                          headercheck__dune_gdt_operators_l2.hh
                          headercheck__dune_gdt_operators_oswaldinterpolation.hh
                          headercheck__dune_gdt_playground_spaces_mapper_restricted.hh
                          headercheck__dune_gdt_prolongations.hh
                          headercheck__dune_gdt_spaces_basefunctionset_product.hh
                          headercheck__dune_gdt_spaces_cg_interface.hh
                          headercheck__dune_gdt_spaces_dg.lib.hh
                          headercheck__dune_gdt_spaces_fv_default.hh
                          headercheck__dune_gdt_test_hyperbolic_problems_interface.hh
                          headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_basisfunctions_legendre.hh
                          headercheck__dune_gdt_test_hyperbolic_problems_transport.hh
                          headercheck__dune_gdt_test_prolongations_l2-local.hh
                          headercheck__python_dune_gdt_spaces_interface.hh
                          test_mpi_operators__elliptic__localizable_product
                          test_mpi_operators__laplace__matrix_operator
                          test_mpi_operators__laplace__operator
                          test_mpi_prolongations__local)
set_tests_properties(test_mpi_operators__elliptic__localizable_product PROPERTIES LABELS "builder_8")
set_tests_properties(test_mpi_operators__laplace__matrix_operator PROPERTIES LABELS "builder_8")
set_tests_properties(test_mpi_operators__laplace__operator PROPERTIES LABELS "builder_8")
set_tests_properties(test_mpi_prolongations__local PROPERTIES LABELS "builder_8")
add_custom_target(test_binaries_builder_9
                  DEPENDS headercheck__dune_gdt_assembler_codim0-matrix-datahandle.hh
                          headercheck__dune_gdt_assembler_local-assemblers.hh
                          headercheck__dune_gdt_assembler_system.hh
                          headercheck__dune_gdt_functionals_interfaces.hh
                          headercheck__dune_gdt_local_functionals_integrals.hh
                          headercheck__dune_gdt_local_integrands_product.hh
                          headercheck__dune_gdt_operators_darcy.hh
                          headercheck__dune_gdt_operators_fv_godunov.hh
                          headercheck__dune_gdt_playground_spaces_block.hh
                          headercheck__dune_gdt_playground_spaces_mapper_block.hh
                          headercheck__dune_gdt_projections_l2-global.hh
                          headercheck__dune_gdt_prolongations_l2-local.hh
                          headercheck__dune_gdt_spaces_fv.hh
                          headercheck__dune_gdt_spaces_localview.hh
                          headercheck__dune_gdt_test_hyperbolic_discretizers_base.hh
                          headercheck__dune_gdt_test_hyperbolic_eocexpectations-fv-burgers-1dyaspgrid.hh
                          headercheck__dune_gdt_test_hyperbolic_problems.hh
                          headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_fokkerplanck_onebeam.hh
                          headercheck__dune_gdt_test_linearelliptic_eocexpectations_cg-spe10-2dalugrid.hh
                          headercheck__dune_gdt_test_operators_elliptic.hh
                          headercheck__dune_gdt_type_traits.hh
                          headercheck__python_dune_gdt_functionals_base.hh
                          headercheck__python_dune_gdt_local_diffusive-flux-estimation-operator.hh
                          headercheck__python_dune_gdt_operators_elliptic_bindings.hh
                          headercheck__python_dune_gdt_playground_operators_RS2017.hh
                          test_mpi_operators__weighted_l2__matrix_operator
                          test_mpi_projections__lagrange_and_global
                          test_prolongations_rt_fv)
set_tests_properties(test_mpi_operators__weighted_l2__matrix_operator PROPERTIES LABELS "builder_9")
set_tests_properties(test_mpi_projections__lagrange_and_global PROPERTIES LABELS "builder_9")
set_tests_properties(test_prolongations_rt_fv PROPERTIES LABELS "builder_9")
add_custom_target(test_binaries_builder_10
                  DEPENDS headercheck__dune_gdt_local_integrands_elliptic-ipdg.hh
                          headercheck__dune_gdt_operators_fv_realizability.hh
                          headercheck__dune_gdt_operators_weighted-l2.hh
                          headercheck__dune_gdt_playground_spaces_restricted.hh
                          headercheck__dune_gdt_playground_timestepper_rosenbrock.hh
                          headercheck__dune_gdt_spaces_mapper_default.hh
                          headercheck__dune_gdt_spaces_rt_default.hh
                          headercheck__dune_gdt_test_hyperbolic_all_eocexpectations.hh
                          headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_fokkerplanck_rectangularic.hh
                          headercheck__dune_gdt_test_linearelliptic_eocexpectations_swipdg-ao2013-2dalugrid.hh
                          headercheck__dune_gdt_test_linearelliptic_eocexpectations_swipdg-spe10-2dyaspgrid.hh
                          headercheck__dune_gdt_test_linearelliptic_estimators_swipdg-fluxreconstruction.hh
                          headercheck__dune_gdt_test_linearelliptic_problems_interface.hh
                          headercheck__dune_gdt_test_projections_l2.hh
                          headercheck__dune_gdt_test_prolongations_lagrange.hh
                          headercheck__dune_gdt_test_prolongations_prolongations.hh
                          headercheck__dune_gdt_test_spaces_fv.hh
                          headercheck__dune_gdt_test_stationary-eocstudy.hh
                          headercheck__python_dune_gdt_functionals_l2_bindings.hh
                          headercheck__python_dune_gdt_projections_bindings.hh
                          headercheck__python_dune_gdt_projections_dirichlet.hh
                          test_mpi_operators__l2__matrix_operator
                          test_mpi_projections_part_0_a
                          test_mpi_prolongations__lagrange_and_global)
set_tests_properties(test_mpi_operators__l2__matrix_operator PROPERTIES LABELS "builder_10")
set_tests_properties(test_mpi_projections_part_0_a PROPERTIES LABELS "builder_10")
set_tests_properties(test_mpi_prolongations__lagrange_and_global PROPERTIES LABELS "builder_10")
add_custom_target(test_binaries_builder_11
                  DEPENDS headercheck__dune_gdt_discretizations_default_stationary-containerbased.hh
                          headercheck__dune_gdt_functionals_elliptic-ipdg.hh
                          headercheck__dune_gdt_local_integrands_fv.hh
                          headercheck__dune_gdt_local_integrands_interfaces.hh
                          headercheck__dune_gdt_local_operators_lagrange-projection.hh
                          headercheck__dune_gdt_operators_fv_enums.hh
                          headercheck__dune_gdt_operators_laplace.hh
                          headercheck__dune_gdt_playground_local_integrands_OS2014.hh
                          headercheck__dune_gdt_projections.hh
                          headercheck__dune_gdt_spaces_dg_interface.hh
                          headercheck__dune_gdt_spaces_fv_interface.hh
                          headercheck__dune_gdt_test_hyperbolic_problems_burgers.hh
                          headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_lebedevquadrature.hh
                          headercheck__dune_gdt_test_linearelliptic_eocexpectations_swipdg-ao2013-2dyaspgrid.hh
                          headercheck__dune_gdt_test_linearelliptic_eocexpectations_swipdg-mixedboundary-2dyaspgrid.hh
                          headercheck__dune_gdt_test_linearelliptic_eocexpectations_swipdg-spe10-2dalugrid.hh
                          headercheck__dune_gdt_test_linearelliptic_swipdg-estimator-expectations.hh
                          headercheck__dune_gdt_test_linearelliptic_swipdg-estimator-testcases.hh
                          headercheck__dune_gdt_test_linearelliptic_swipdg-testcases.hh
                          headercheck__dune_gdt_test_operators_base.hh
                          headercheck__dune_gdt_test_projections_projections.hh
                          headercheck__dune_gdt_test_spaces_cg.hh
                          headercheck__dune_gdt_test_stationary-testcase.hh
                          headercheck__dune_gdt_timestepper_interface.hh
                          test_mpi_projections_part_0_b
                          test_mpi_projections_part_1_b)
set_tests_properties(test_mpi_projections_part_0_b PROPERTIES LABELS "builder_11")
set_tests_properties(test_mpi_projections_part_1_b PROPERTIES LABELS "builder_11")
add_custom_target(test_binaries_builder_12
                  DEPENDS headercheck__dune_gdt_discretizations_default_fv.hh
                          headercheck__dune_gdt_local_fluxes_entropybased.hh
                          headercheck__dune_gdt_local_fluxes_laxwendroff.hh
                          headercheck__dune_gdt_local_integrands_ESV2007.hh
                          headercheck__dune_gdt_projections_lagrange.hh
                          headercheck__dune_gdt_spaces_cg.hh
                          headercheck__dune_gdt_spaces_parallel.hh
                          headercheck__dune_gdt_spaces_product.hh
                          headercheck__dune_gdt_test_hyperbolic_discretizers_fv.hh
                          headercheck__dune_gdt_test_hyperbolic_eocexpectations_base.hh
                          headercheck__dune_gdt_test_hyperbolic_problems_base.hh
                          headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_kinetictransport_linesource.hh
                          headercheck__dune_gdt_test_hyperbolic_problems_momentmodels_triangulation.hh
                          headercheck__dune_gdt_test_linearelliptic_eocexpectations_cg-er2007-2dalugrid.hh
                          headercheck__dune_gdt_test_linearelliptic_eocexpectations_swipdg-esv2007-2dyaspgrid.hh
                          headercheck__dune_gdt_test_linearelliptic_problems_spe10.hh
                          headercheck__dune_gdt_test_projections_l2-global.hh
                          headercheck__dune_gdt_test_spaces_base.hh
                          headercheck__dune_gdt_test_spaces_dg.hh
                          headercheck__dune_gdt_timestepper_implicit-rungekutta.hh
                          headercheck__python_dune_gdt_playground_operators_ESV2007.hh
                          test_mpi_projections__local
                          test_mpi_projections_part_2_b
                          test_spaces_cg)
set_tests_properties(test_mpi_projections__local PROPERTIES LABELS "builder_12")
set_tests_properties(test_mpi_projections_part_2_b PROPERTIES LABELS "builder_12")
set_tests_properties(test_spaces_cg PROPERTIES LABELS "builder_12")
