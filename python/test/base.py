
import pytest
from dune.xt.common.test import load_all_submodule
from dune.xt.grid import make_cube_dd_subdomains_grid__2d_simplex_aluconformgrid
from dune.xt.grid import make_boundary_info_on_dd_subdomain_layer
from dune.xt.la import IstlRowMajorSparseMatrixDouble as Matrix
from dune.gdt.__operators_elliptic_ipdg import make_elliptic_swipdg_affine_factor_matrix_operator as make_elliptic_swipdg_matrix_operator
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
    return make_cube_dd_subdomains_grid__2d_simplex_aluconformgrid(
            lower_left=[-1,-1],
            upper_right=[1,1],
            num_elements=[num_subdomains[0]*half_num_fine_elements_per_subdomain_and_dim,
                          num_subdomains[1]*half_num_fine_elements_per_subdomain_and_dim],
            num_refinements=2,
            num_partitions=num_subdomains,
            num_oversampling_layers=2*half_num_fine_elements_per_subdomain_and_dim,
            inner_boundary_segment_index=inner_boundary_segment_index)


def test_blockspace():
    from dune.gdt.__spaces_block import make_block_dg_dd_subdomain_view_to_1x1_gdt_p1_space as make_block_space
    grid = make_grid()
    block_space = make_block_space(grid)
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


if __name__ == '__main__':
    pytest.main()