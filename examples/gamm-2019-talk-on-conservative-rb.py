import numpy as np

from pymor.basic import *
from pymor.bindings.dunext import DuneXTMatrixOperator, DuneXTVectorSpace
from pymor.bindings.dunegdt import DuneGDTVisualizer

from dune.xt.la import IstlDenseVectorDouble
from dune.xt.functions import ConstantFunction__2d_to_1x1 as ConstantFunction
from dune.xt.functions import ExpressionFunction__2d_to_1x1 as ExpressionFunction

logger = getLogger('main.main');
set_log_levels({'main': 'INFO',
                'pymor.discretizations': 'WARN'})

diffusion = {'functions': (ExpressionFunction('x', ['1+cos(0.5*pi*x[0])*cos(0.5*pi*x[1])'], 3, 'lambda_1'),
                           ExpressionFunction('x', [ '-cos(0.5*pi*x[0])*cos(0.5*pi*x[1])'], 3, 'lambda_2')),
             'coefficients': (1., ProjectionParameterFunctional('switch', ()))}
diffusion_expression = '1+(1-{})*cos(0.5*pi*x[0])*cos(0.5*pi*x[1])'
f = ExpressionFunction('x', ['0.5*pi*pi*cos(0.5*pi*x[0])*cos(0.5*pi*x[1])'], 3, 'lambda_2')
zero = ConstantFunction(0.)

from dune.gdt.gamm_2019_talk_on_conservative_rb import (
        DiscontinuousLagrangeSpace,
        GridProvider,
        RaviartThomasSpace,
        assemble_DG_product_matrix,
        assemble_Hdiv_product_matrix,
        assemble_L2_vector,
        assemble_SWIPDG_matrix,
        compute_flux_reconstruction,
        compute_local_conservation_error,
        make_discrete_function,
        visualize,
        )

def make_marix_operator(mat, ID):
    return DuneXTMatrixOperator(mat, source_id=ID, range_id=ID)

grid = GridProvider([-1, -1], [1, 1], [4, 4]) # The ESV2007 test is [-1, 1]^2, 4x4 elements, ...
grid.refine(2) # 2 refinements with a conforming alugrid
dg_space = DiscontinuousLagrangeSpace(grid, 1)

lhs_op = LincombOperator([make_marix_operator(assemble_SWIPDG_matrix(dg_space, diff), 'PRESSURE')
                          for diff in diffusion['functions']],
                         diffusion['coefficients'])
rhs_func = VectorFunctional(lhs_op.range.make_array((assemble_L2_vector(dg_space, f),)))
dg_product = make_marix_operator(assemble_DG_product_matrix(dg_space), 'PRESSURE')

fom = StationaryDiscretization(lhs_op, rhs_func, products={'h1_penalty': dg_product},
                               visualizer=DuneGDTVisualizer(dg_space))
fom = fom.with_(parameter_space=CubicParameterSpace(fom.parameter_type, 0.1, 1.))

reductor = CoerciveRBReductor(fom,
        product=fom.h1_penalty_product,
        coercivity_estimator=ExpressionParameterFunctional('switch', fom.parameter_type))
training_set = fom.parameter_space.sample_uniformly(100)
greedy_data = greedy(fom, reductor, training_set, extension_params={'method': 'gram_schmidt'}, max_extensions=1)
rom = greedy_data['rd']

# build reduced H_div_0 basis, therefore ...
rtn_space = RaviartThomasSpace(grid, 0)
FluxVectorSpace = DuneXTVectorSpace(IstlDenseVectorDouble, rtn_space.num_DoFs, 'FLUX')
rtn_product = make_marix_operator(assemble_Hdiv_product_matrix(rtn_space), 'FLUX')
# - solve for the affine shift (by solving with diffusion one)
u_h_f = make_marix_operator(assemble_SWIPDG_matrix(dg_space, ConstantFunction(1.)),
                            'PRESSURE').apply_inverse(rhs_func.as_source_array())
u_h_f = u_h_f._list[0].impl
make_discrete_function(dg_space, u_h_f, 'u_h_f').visualize('u_h_f')
t_h_f = compute_flux_reconstruction(grid, dg_space, rtn_space, ConstantFunction(1.), u_h_f)
logger.info('conservation_error (t_h_f): {}'.format(compute_local_conservation_error(
    grid, make_discrete_function(rtn_space, t_h_f, 't_h_f'), f)))
del u_h_f
make_discrete_function(rtn_space, t_h_f, 't_h_f').visualize('t_h_f')
# - compute the reconstruction for each param
Hdiv_0_RB = FluxVectorSpace.empty()
for ii in range(len(greedy_data['max_err_mus'])):
    mu = greedy_data['max_err_mus'][ii]
    u_h = fom.solve(mu)._list[0].impl
    make_discrete_function(dg_space, u_h, 'u_h_f').visualize('u_h__mu_{}'.format(ii))
    diff = ExpressionFunction('x', [diffusion_expression.format(mu['switch'])], 3, 'diff')
    visualize(grid, diff, 'diffusion__mu_{}'.format(ii))
    t_h = compute_flux_reconstruction(grid, dg_space, rtn_space, diff, u_h)
    # logger.info('conservation_error (t_h, mu_{}): {}'.format(ii, compute_local_conservation_error(
        # grid, make_discrete_function(rtn_space, t_h, 't_h'), f)))
    make_discrete_function(rtn_space, t_h, 't_h').visualize('t_h__mu_{}'.format(ii))
    t_h_0 = t_h - t_h_f
    # logger.info('conservation_error (t_h_0, mu_{}): {}'.format(ii, compute_local_conservation_error(
        # grid, make_discrete_function(rtn_space, t_h_0, 't_h_0'), zero)))
    make_discrete_function(rtn_space, t_h_0, 't_h_0').visualize('t_h_0__mu_{}'.format(ii))
    Hdiv_0_RB.append(FluxVectorSpace.from_data([t_h_0,]))
Hdiv_0_RB = gram_schmidt(Hdiv_0_RB, product=rtn_product, copy=False)


logger.info('testing conservation property ...')
for mu in fom.parameter_space.sample_randomly(10):
    u_RB = rom.solve(mu)
    u_RB = reductor.reconstruct(u_RB)._list[0].impl
    diff = ExpressionFunction('x', [diffusion_expression.format(mu['switch'])], 3, 'diff')
    t_RB = compute_flux_reconstruction(grid, dg_space, rtn_space, diff, u_RB)
    logger.info('conservation_error (t_RB, mu_{}): {}'.format(mu['switch'], compute_local_conservation_error(
        grid, make_discrete_function(rtn_space, t_RB, 't_RB'), f)))
    # project onto Hdiv_0_RB
    t_RB_0 = rtn_product.apply2(Hdiv_0_RB, FluxVectorSpace.from_data([t_RB - t_h_f,]))
    t_RB_0 = Hdiv_0_RB.lincomb(t_RB_0[:,0])
    t_RB_f = t_RB_0._list[0].impl + t_h_f
    logger.info('conservation_error (t_RB_f, mu_{}): {}'.format(mu['switch'], compute_local_conservation_error(
        grid, make_discrete_function(rtn_space, t_RB_f, 't_RB_f'), f)))
    # compare to correct flux
    u_h = fom.solve(mu)._list[0].impl
    t_h = compute_flux_reconstruction(grid, dg_space, rtn_space, diff, u_h)
    def rtn_norm(t):
        t = FluxVectorSpace.from_data([t,])
        return np.sqrt(rtn_product.apply2(t, t)[0][0])
    logger.info('relative flux reconstruction error: {}'.format(rtn_norm(t_h - t_RB_f) / rtn_norm(t_h)))

