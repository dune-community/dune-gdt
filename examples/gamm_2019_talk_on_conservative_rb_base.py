from pymor.basic import *

logger = getLogger('main.simulate_single_greedy_step');

def simulate_single_greedy_step(
        fom,
        dg_product,
        FluxVectorSpace,
        rtn_product,
        t_h_f,
        compute_flux_reconstruction,
        compute_estimate,
        compute_reference_error,
        max_extensions,
        num_samples):

    logger.info('building pressure RB ...')
    logger.info('')

    reductor = CoerciveRBReductor(fom,
            product=dg_product,
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
        t_h = FluxVectorSpace.from_data([compute_flux_reconstruction(mu, u_h._list[0].impl),])
        Hdiv_0_RB.append(t_h - t_h_f)
    Hdiv_0_RB = gram_schmidt(Hdiv_0_RB, product=rtn_product, copy=False)

    logger.info('')
    logger.info('testing for some parameters ...')

    mus = []
    etas = []
    eta_NCs = []
    eta_Rs = []
    eta_DFs = []
    errors = []
    efficiencies = []

    for mu in fom.parameter_space.sample_uniformly(num_samples):
        logger.info('  mu = {}'.format(mu))
        u_RB = rom.solve(mu)
        u_RB = reductor.reconstruct(u_RB)
        t_RB = FluxVectorSpace.from_data([compute_flux_reconstruction(mu, u_RB._list[0].impl),])
        # project onto Hdiv_0_RB (we know that Hdiv_0_R is orthonormalized wrt rtn_product)
        t_RB_0 = rtn_product.apply2(Hdiv_0_RB, FluxVectorSpace.from_data([t_RB._list[0].impl - t_h_f._list[0].impl,]))
        t_RB_0 = Hdiv_0_RB.lincomb(t_RB_0[:,0])
        t_RB_f = t_RB_0._list[0].impl + t_h_f._list[0].impl
        logger.info('    estimating error ...')
        eta, eta_NC, eta_R, eta_DF = compute_estimate(mu, u_RB._list[0].impl, t_RB_f)
        logger.info('      eta_NC = {}'.format(eta_NC))
        logger.info('      eta_R  = {}'.format(eta_R))
        logger.info('      eta_DF = {}'.format(eta_DF))
        logger.info('      eta    = {}'.format(eta))
        logger.info('    computing error by reference solution ...')
        error = compute_reference_error(mu, u_RB._list[0].impl)
        logger.info('      error  = {}'.format(error))
        logger.info('   => efficiency = {}'.format(eta/error))
        mus.append(mu)
        etas.append(eta)
        eta_NCs.append(eta_NC)
        eta_Rs.append(eta_R)
        eta_DFs.append(eta_DF)
        errors.append(error)
        efficiencies.append(eta/error)

    return mus, etas, eta_NCs, eta_Rs, eta_DFs, errors, efficiencies

