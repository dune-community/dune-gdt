import numpy as np

from pymor.basic import *
from pymor.bindings.dunext import DuneXTMatrixOperator, DuneXTVectorSpace
from pymor.bindings.dunegdt import DuneGDTVisualizer
from pymor.vectorarrays.list import ListVectorArray

from dune.xt.la import IstlDenseVectorDouble
from dune.xt.functions import ConstantFunction__2d_to_1x1 as ConstantFunction
from dune.xt.functions import ExpressionFunction__2d_to_1x1 as ExpressionFunction

logger = getLogger('main.main');
set_log_levels({'main': 'INFO',
                'pymor.discretizations': 'WARN'})

diffusion = {'functions': (ExpressionFunction('x', ['1+cos(0.5*pi*x[0])*cos(0.5*pi*x[1])'], 3, 'lambda_1'),
                           ExpressionFunction('x', [ '-cos(0.5*pi*x[0])*cos(0.5*pi*x[1])'], 3, 'lambda_2')),
                           'coefficients': (ExpressionParameterFunctional('1', {'switch': ()}),
                                            ProjectionParameterFunctional('switch', ()))}
diffusion_expression = '1+(1-{})*cos(0.5*pi*x[0])*cos(0.5*pi*x[1])'
mu_bar = {'switch': 0.1}
mu_hat = {'switch': 0.1}
diffusion_bar = ExpressionFunction('x', [diffusion_expression.format(mu_bar['switch'])], 3, 'diffusion_mu_bar')
diffusion_hat = ExpressionFunction('x', [diffusion_expression.format(mu_hat['switch'])], 3, 'diffusion_mu_hat')
f = ExpressionFunction('x', ['0.5*pi*pi*cos(0.5*pi*x[0])*cos(0.5*pi*x[1])'], 3, 'lambda_2')
zero = ConstantFunction(0.)

def alpha(mu, mu_bar):
    return np.min([theta.evaluate(mu) / theta.evaluate(mu_bar) for theta in diffusion['coefficients']])

def gamma(mu, mu_bar):
    return np.max([theta.evaluate(mu) / theta.evaluate(mu_bar) for theta in diffusion['coefficients']])

from dune.gdt.gamm_2019_talk_on_conservative_rb import (
        DiscontinuousLagrangeSpace,
        GridProvider,
        RaviartThomasSpace,
        assemble_energy_semi_product_matrix,
        assemble_DG_product_matrix,
        assemble_Hdiv_product_matrix,
        assemble_L2_vector,
        assemble_SWIPDG_matrix,
        compute_estimate,
        compute_flux_reconstruction,
        compute_local_conservation_error,
        make_discrete_function,
        prolong,
        visualize,
        )

def make_marix_operator(mat, ID):
    return DuneXTMatrixOperator(mat, source_id=ID, range_id=ID)

logger.info('discretizing ...')

def discretize(num_refinements):
    grid = GridProvider([-1, -1], [1, 1], [4, 4]) # The ESV2007 test is [-1, 1]^2, 4x4 elements, ...
    grid.refine(num_refinements)
    dg_space = DiscontinuousLagrangeSpace(grid, 1)

    lhs_op = LincombOperator([make_marix_operator(assemble_SWIPDG_matrix(dg_space, diff), 'PRESSURE')
                              for diff in diffusion['functions']],
                             diffusion['coefficients'])
    rhs_func = VectorFunctional(lhs_op.range.make_array((assemble_L2_vector(dg_space, f),)))
    dg_product = make_marix_operator(assemble_DG_product_matrix(dg_space), 'PRESSURE')

    fom = StationaryDiscretization(lhs_op, rhs_func, products={'h1_penalty': dg_product},
                                   visualizer=DuneGDTVisualizer(dg_space))
    fom = fom.with_(parameter_space=CubicParameterSpace(fom.parameter_type, 0.1, 1.))
    return grid, dg_space, dg_product, fom

grid, dg_space, dg_product, fom = discretize(2) # ... and 2 refinements with ALU_2D_SIMPLEX_CONFORMING
PressureVectorSpace = DuneXTVectorSpace(IstlDenseVectorDouble, dg_space.num_DoFs, 'PRESSURE')

logger.info('grid has {} elements'.format(grid.num_elements))
logger.info('space has {} DoFs'.format(dg_space.num_DoFs))
logger.info('')
logger.info('computing reference discretization ...')

reference_grid, reference_dg_space, _, reference_fom = discretize(2 + 3*2)
reference_energy_semi_product = make_marix_operator(
        assemble_energy_semi_product_matrix(reference_dg_space, diffusion_bar), 'PRESSURE')
ReferencePressureVectorSpace = DuneXTVectorSpace(IstlDenseVectorDouble, reference_dg_space.num_DoFs, 'PRESSURE')

def reference_dg_norm(u):
    if not isinstance(u, ListVectorArray):
        u = ReferencePressureVectorSpace.from_data([u,])
    return np.sqrt(reference_energy_semi_product.apply2(u, u)[0][0])

logger.info('')
logger.info('assembling Hdiv product ...')

rtn_space = RaviartThomasSpace(grid, 0)
FluxVectorSpace = DuneXTVectorSpace(IstlDenseVectorDouble, rtn_space.num_DoFs, 'FLUX')
rtn_product = make_marix_operator(assemble_Hdiv_product_matrix(rtn_space), 'FLUX')

def rtn_norm(t):
    if not isinstance(t, ListVectorArray):
        t = FluxVectorSpace.from_data([t,])
    return np.sqrt(rtn_product.apply2(t, t)[0][0])

logger.info('')
logger.info('computing ESV2007 pressure and flux ...')

u_h_f = make_marix_operator(assemble_SWIPDG_matrix(dg_space, ConstantFunction(1.)),
                            'PRESSURE').apply_inverse(fom.rhs.as_source_array())
t_h_f = FluxVectorSpace.from_data([
    compute_flux_reconstruction(grid, dg_space, rtn_space, ConstantFunction(1.), u_h_f._list[0].impl),])

logger.info('')
logger.info('computing [OS2015, table 1] estimates (should be 0.166, 0.723, 0.355) ...')

_, eta_NC, eta_R, eta_DF = compute_estimate(
        grid,
        make_discrete_function(dg_space, u_h_f._list[0].impl, 'u_h_f'),
        make_discrete_function(rtn_space, t_h_f._list[0].impl, 't_h_f'),
        f, ConstantFunction(1.), ConstantFunction(1.), ConstantFunction(1.),
        1, 1, 1)

logger.info('    are {}, {}, {}'.format(eta_NC, eta_R, eta_DF))
logger.info('')
logger.info('computing other OS2015 estimates (should be '
            '[table 3, eta_NC] 0.182, <= [table 1, eta_R] 0.166, [table 2, eta_DF] 0.316) ...')

mu = {'switch': 1}
diffusion_mu = ExpressionFunction('x', [diffusion_expression.format(mu['switch'])], 3, 'diffusion_mu')
u_h = fom.solve(mu)
t_h = FluxVectorSpace.from_data([
    compute_flux_reconstruction(grid, dg_space, rtn_space, diffusion_mu, u_h._list[0].impl),])
_, eta_NC, eta_R, eta_DF = compute_estimate(
        grid,
        make_discrete_function(dg_space, u_h._list[0].impl, 'u_h'),
        make_discrete_function(rtn_space, t_h._list[0].impl, 't_h'),
        f, diffusion_mu, diffusion_bar, diffusion_hat,
        alpha(mu, mu_bar), alpha(mu, mu_hat), gamma(mu, mu_bar))

logger.info('    are {}, {}, {}'.format(eta_NC, eta_R, eta_DF))
logger.info('')
max_extensions=1
logger.info('building pressure RB (simulating intermediate greedy step by max_extensions={}) ...'.format(
    max_extensions))
logger.info('')

reductor = CoerciveRBReductor(fom,
        product=fom.h1_penalty_product,
        coercivity_estimator=ExpressionParameterFunctional('switch', fom.parameter_type))
training_set = fom.parameter_space.sample_uniformly(100)
greedy_data = greedy(fom, reductor, training_set,
        extension_params={'method': 'gram_schmidt'}, max_extensions=max_extensions)
rom = greedy_data['rd']
greedy_mus = greedy_data['max_err_mus']

logger.info('')
logger.info('building Hdiv_0 RB ...')

Hdiv_0_RB = FluxVectorSpace.empty()
for ii in range(len(greedy_data['max_err_mus'])):
    mu = greedy_data['max_err_mus'][ii]
    u_h = fom.solve(mu)
    diffusion_mu = ExpressionFunction('x', [diffusion_expression.format(mu['switch'])], 3, 'diffusion_mu')
    t_h = FluxVectorSpace.from_data([
        compute_flux_reconstruction(grid, dg_space, rtn_space, diffusion_mu, u_h._list[0].impl),])
    Hdiv_0_RB.append(t_h - t_h_f)
Hdiv_0_RB = gram_schmidt(Hdiv_0_RB, product=rtn_product, copy=False)

logger.info('')
logger.info('testing for some parameters ...')

for mu in fom.parameter_space.sample_uniformly(3):
    logger.info('  mu = {}'.format(mu))
    u_RB = rom.solve(mu)
    u_RB = reductor.reconstruct(u_RB)
    diffusion_mu = ExpressionFunction('x', [diffusion_expression.format(mu['switch'])], 3, 'diffusion_mu')
    t_RB = FluxVectorSpace.from_data([
        compute_flux_reconstruction(grid, dg_space, rtn_space, diffusion_mu, u_RB._list[0].impl),])
    # logger.info('    conservation_error(t_RB): {}'.format(compute_local_conservation_error(
        # grid, make_discrete_function(rtn_space, t_RB._list[0].impl, 't_RB'), f)))
    # project onto Hdiv_0_RB (we know that Hdiv_0_R is orthonormalized wrt rtn_product)
    t_RB_0 = rtn_product.apply2(Hdiv_0_RB, FluxVectorSpace.from_data([t_RB._list[0].impl - t_h_f._list[0].impl,]))
    t_RB_0 = Hdiv_0_RB.lincomb(t_RB_0[:,0])
    t_RB_f = t_RB_0._list[0].impl + t_h_f._list[0].impl
    # logger.info('    conservation_error(t_RB_f): {}'.format(compute_local_conservation_error(
        # grid, make_discrete_function(rtn_space, t_RB_f, 't_RB_f'), f)))
    # compare to correct flux
    # u_h = fom.solve(mu)._list[0].impl
    # t_h = compute_flux_reconstruction(grid, dg_space, rtn_space, diffusion_mu, u_h)
    # logger.info('    relative flux reconstruction error: {}'.format(rtn_norm(t_h - t_RB_f) / rtn_norm(t_h)))
    # compute estimate
    logger.info('    estimating error ...')
    eta, eta_NC, eta_R, eta_DF = compute_estimate(
            grid,
            make_discrete_function(dg_space, u_RB._list[0].impl, 'u_RB'),
            make_discrete_function(rtn_space, t_RB_f, 't_RB_f'),
            f, diffusion_mu, diffusion_bar, diffusion_hat,
            alpha(mu, mu_bar), alpha(mu, mu_hat), gamma(mu, mu_bar))
    logger.info('      eta_NC = {}'.format(eta_NC))
    logger.info('      eta_R  = {}'.format(eta_R))
    logger.info('      eta_DF = {}'.format(eta_DF))
    logger.info('      eta    = {}'.format(eta))
    logger.info('    computing error by reference solution ...')
    u_RB_on_ref = prolong(dg_space, u_RB._list[0].impl, reference_dg_space)
    u_h_ref = reference_fom.solve(mu)
    error = reference_dg_norm(u_h_ref._list[0].impl - u_RB_on_ref)
    logger.info('      error = {}'.format(error))
    logger.info('   => efficiency = {}'.format(eta/error))

