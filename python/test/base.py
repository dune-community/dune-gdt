# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Rene Milk (2018)
# ~~~

import pytest
from dune.xt.common.test import load_all_submodule
from dune.xt.grid.provider import make_cube_dd_subdomains_grid__2d_cube_yaspgrid
from dune.xt.grid.boundaryinfo import make_boundary_info_on_dd_subdomain_layer
from dune.xt.la import IstlRowMajorSparseMatrixDouble as Matrix
from dune.xt.functions import (
    make_constant_function_2x2,
    make_expression_function_1x1
)


def test_load_all():
    import dune.gdt as gdt
    load_all_submodule(gdt)


def make_grid(num_subdomains=[2, 2]):
    half_num_fine_elements_per_subdomain_and_dim = 3
    inner_boundary_segment_index = 18446744073709551573
    import dune.xt.grid.provider as prv
    import dune.gdt
    maker = getattr(prv, 'make_cube_dd_subdomains_grid__{}'.format(dune.gdt.GDT_BINDINGS_GRID))
    return maker(
            lower_left=[-1,-1],
            upper_right=[1,1],
            num_elements=[num_subdomains[0]*half_num_fine_elements_per_subdomain_and_dim,
                          num_subdomains[1]*half_num_fine_elements_per_subdomain_and_dim],
            num_refinements=2,
            num_partitions=num_subdomains,
            num_oversampling_layers=2*half_num_fine_elements_per_subdomain_and_dim,
            inner_boundary_segment_index=inner_boundary_segment_index)


def test_blockspace():
    from dune.gdt.spaces import make_block_dg_space
    from dune.gdt.__operators_elliptic_ipdg import \
        make_elliptic_swipdg_affine_factor_matrix_operator as make_elliptic_swipdg_matrix_operator
    grid = make_grid()
    block_space = make_block_dg_space(grid)
    diffusion = make_expression_function_1x1(grid, 'x', '1', order=2, name='lambda_0')
    kappa = make_constant_function_2x2(grid, [[1., 0.], [0., 1.]], name='kappa')
    local_all_neumann_boundary_info = make_boundary_info_on_dd_subdomain_layer(grid, {'type': 'xt.grid.boundaryinfo.allneumann'})

    local_matrices = [None] * grid.num_subdomains
    local_patterns = [block_space.local_space(ii).compute_pattern('face_and_volume')
                      for ii in range(block_space.num_blocks)]
    for ii in range(grid.num_subdomains):
        local_matrices[ii] = Matrix(block_space.local_space(ii).size(),
                                    block_space.local_space(ii).size(),
                                    local_patterns[ii])
        ipdg_operator = make_elliptic_swipdg_matrix_operator(diffusion_factor=diffusion, diffusion_tensor=kappa,
                                                             boundary_info=local_all_neumann_boundary_info,
                                                             matrix=local_matrices[ii],
                                                             space=block_space.local_space(ii), over_integrate=2)
        ipdg_operator.assemble()


def test_visualize():
    from dune.gdt.spaces import make_dg_space
    from dune.gdt.discretefunction import make_discrete_function
    grid = make_grid()
    space = make_dg_space(grid)
    df = make_discrete_function(space, 'test')
    df.visualize(filename='foo')


if __name__ == '__main__':
    pytest.main()