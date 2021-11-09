# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2017 - 2018)
#   René Fritze     (2016, 2018)
# ~~~

from tempfile import NamedTemporaryFile

from dune.xt import guarded_import
from dune.xt.common.config import config

from ._version import __version__

for mod_name in (     # order should not matter!
        '_discretefunction_bochner',
        '_discretefunction_discretefunction',
        '_discretefunction_dof_vector',
        '_functionals_interfaces_common',
        '_functionals_interfaces_eigen',
        '_functionals_interfaces_istl',
        '_functionals_vector_based',
        '_interpolations_boundary',
        '_interpolations_default',
        '_interpolations_oswald',
        '_local_bilinear_forms_coupling_intersection_integrals',
        '_local_bilinear_forms_coupling_intersection_interface',
        '_local_bilinear_forms_element_integrals',
        '_local_bilinear_forms_element_interface',
        '_local_bilinear_forms_intersection_integrals',
        '_local_bilinear_forms_intersection_interface',
        '_local_bilinear_forms_restricted_coupling_intersection_integrals',
        '_local_bilinear_forms_restricted_intersection_integrals',
        '_local_bilinear_forms_vectorized_element_integrals',
        '_local_functionals_element_integrals',
        '_local_functionals_element_interface',
        '_local_functionals_intersection_integrals',
        '_local_functionals_intersection_interface',
        '_local_functionals_restricted_intersection_integrals',
        '_local_functionals_vectorized_element_integrals',
        '_local_integrands_binary_element_interface',
        '_local_integrands_binary_intersection_interface',
        '_local_integrands_element_product',
        '_local_integrands_intersection_product',
        '_local_integrands_ipdg_boundary_penalty',
        '_local_integrands_ipdg_inner_penalty',
        '_local_integrands_jump_boundary',
        '_local_integrands_jump_inner',
        '_local_integrands_laplace',
        '_local_integrands_laplace_ipdg_dirichlet_coupling',
        '_local_integrands_laplace_ipdg_inner_coupling',
        '_local_integrands_linear_advection',
        '_local_integrands_linear_advection_upwind_dirichlet_coupling',
        '_local_integrands_linear_advection_upwind_inner_coupling',
        '_local_integrands_quaternary_intersection_interface',
        '_local_integrands_unary_element_interface',
        '_local_integrands_unary_intersection_interface',
        '_local_operators_coupling_intersection_indicator',
        '_local_operators_element_indicator',
        '_local_operators_element_interface',
        '_local_operators_intersection_indicator',
        '_local_operators_intersection_interface',
        '_operators_bilinear_form',
        '_operators_interfaces_common',
        '_operators_interfaces_eigen',
        '_operators_interfaces_istl_1d',
        '_operators_interfaces_istl_2d',
        '_operators_interfaces_istl_3d',
        '_operators_laplace_ipdg_flux_reconstruction',
        '_operators_matrix_based_factory',
        '_operators_operator',
        '_prolongations',
        '_spaces_bochner',
        '_spaces_h1_continuous_lagrange',
        '_spaces_hdiv_raviart_thomas',
        '_spaces_interface',
        '_spaces_l2_discontinuous_lagrange',
        '_spaces_l2_finite_volume',
        '_spaces_skeleton_finite_volume',
        '_tools_adaptation_helper',
        '_tools_dirichlet_constraints',
        '_tools_grid_quality_estimates',
        '_tools_sparsity_pattern',
):
    guarded_import(globals(), 'dune.gdt', mod_name)


if config.HAVE_K3D:
    from dune.xt.common.vtk.plot import plot

    def visualize_function(function, grid=None, subsampling=False):
        assert function.dim_domain <= 2, f'Not implemented yet for {function.dim_domain}-dimensional grids!'
        if function.dim_domain == 1:
            import numpy as np
            from matplotlib import pyplot as plt
            from dune.xt.functions import GridFunction
            from dune.gdt import ContinuousLagrangeSpace, default_interpolation, DiscreteFunction

            assert grid # not optimal
            P1_space = ContinuousLagrangeSpace(grid, order=1)
            interpolation_points = np.array(P1_space.interpolation_points(), copy=False)[:, 0]
            piecewise_linear_interpolant = default_interpolation(GridFunction(grid, function), P1_space)
            values = np.array(piecewise_linear_interpolant.dofs.vector, copy=False)

            plt.figure()
            plt.title(f'{function.name}')
            plt.plot(interpolation_points, values)

            return plt.gca()
        elif function.dim_domain == 2:
            assert function.dim_range == 1, f'Not implemented yet for {function.dim_domain}-dimensional functions!'
            tmpfile = NamedTemporaryFile(mode='wb', delete=False, suffix='.vtu').name
            failed = False
            try:     # discrete function
                function.visualize(filename=tmpfile[:-4])
                return plot(tmpfile, color_attribute_name=function.name)
            except TypeError:
                failed = True
            except AttributeError:
                failed = True

            if failed:
                from dune.xt.functions import visualize_function as visualize_xt_function

                assert grid
                return visualize_xt_function(function, grid, subsampling=subsampling)
