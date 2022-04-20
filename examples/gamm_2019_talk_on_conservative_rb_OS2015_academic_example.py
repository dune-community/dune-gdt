import numpy as np

from pymor.basic import *
from pymor.bindings.dunext import DuneXTMatrixOperator, DuneXTVectorSpace
from pymor.bindings.dunegdt import DuneGDTVisualizer
from pymor.vectorarrays.list import ListVectorArray

from dune.xt.la import IstlDenseVectorDouble
from dune.xt.functions import ConstantFunction__2d_to_1x1 as ConstantFunction
from dune.xt.functions import ExpressionFunction__2d_to_1x1 as ExpressionFunction

logger = getLogger("main.main")
set_log_levels(
    {"main": "INFO", "main.simulate_single_greedy_step": "WARN", "pymor": "WARN"}
)

diffusion = {
    "functions": (
        ExpressionFunction("x", ["1+cos(0.5*pi*x[0])*cos(0.5*pi*x[1])"], 3, "lambda_1"),
        ExpressionFunction("x", ["-cos(0.5*pi*x[0])*cos(0.5*pi*x[1])"], 3, "lambda_2"),
    ),
    "coefficients": (
        ExpressionParameterFunctional("1", {"switch": ()}),
        ProjectionParameterFunctional("switch", ()),
    ),
}
diffusion_expression = "1+(1-{})*cos(0.5*pi*x[0])*cos(0.5*pi*x[1])"
mu_bar = {"switch": 0.1}
mu_hat = {"switch": 0.1}
diffusion_bar = ExpressionFunction(
    "x", [diffusion_expression.format(mu_bar["switch"])], 3, "diffusion_mu_bar"
)
diffusion_hat = ExpressionFunction(
    "x", [diffusion_expression.format(mu_hat["switch"])], 3, "diffusion_mu_hat"
)
f = ExpressionFunction("x", ["0.5*pi*pi*cos(0.5*pi*x[0])*cos(0.5*pi*x[1])"], 3, "f")
zero = ConstantFunction(0.0)


def alpha(mu, mu_bar):
    return np.min(
        [
            theta.evaluate(mu) / theta.evaluate(mu_bar)
            for theta in diffusion["coefficients"]
        ]
    )


def gamma(mu, mu_bar):
    return np.max(
        [
            theta.evaluate(mu) / theta.evaluate(mu_bar)
            for theta in diffusion["coefficients"]
        ]
    )


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


logger.info("discretizing ...")


def discretize(num_refinements):
    grid = GridProvider(
        [-1, -1], [1, 1], [4, 4]
    )  # The ESV2007 test is [-1, 1]^2, 4x4 elements, ...
    grid.refine(num_refinements)
    dg_space = DiscontinuousLagrangeSpace(grid, 1)

    lhs_op = LincombOperator(
        [
            make_marix_operator(assemble_SWIPDG_matrix(dg_space, diff), "PRESSURE")
            for diff in diffusion["functions"]
        ],
        diffusion["coefficients"],
    )
    rhs_func = VectorFunctional(
        lhs_op.range.make_array((assemble_L2_vector(dg_space, f),))
    )
    dg_product = make_marix_operator(assemble_DG_product_matrix(dg_space), "PRESSURE")

    fom = StationaryDiscretization(
        lhs_op,
        rhs_func,
        products={"energy_penalty": dg_product},
        visualizer=DuneGDTVisualizer(dg_space),
    )
    fom = fom.with_(parameter_space=CubicParameterSpace(fom.parameter_type, 0.1, 1.0))
    fom.enable_caching("disk")
    return grid, dg_space, dg_product, fom


grid, dg_space, dg_product, fom = discretize(
    2
)  # ... and 2 refinements with ALU_2D_SIMPLEX_CONFORMING
PressureVectorSpace = DuneXTVectorSpace(
    IstlDenseVectorDouble, dg_space.num_DoFs, "PRESSURE"
)

logger.info("grid has {} elements".format(grid.num_elements))
logger.info("space has {} DoFs".format(dg_space.num_DoFs))
logger.info("computing reference discretization ...")

reference_grid, reference_dg_space, _, reference_fom = discretize(2 + 3 * 2)
reference_energy_semi_product = make_marix_operator(
    assemble_energy_semi_product_matrix(reference_dg_space, diffusion_bar), "PRESSURE"
)
ReferencePressureVectorSpace = DuneXTVectorSpace(
    IstlDenseVectorDouble, reference_dg_space.num_DoFs, "PRESSURE"
)


def reference_dg_norm(u):
    if not isinstance(u, ListVectorArray):
        u = ReferencePressureVectorSpace.from_data(
            [
                u,
            ]
        )
    return np.sqrt(reference_energy_semi_product.apply2(u, u)[0][0])


logger.info("reference grid has {} elements".format(reference_grid.num_elements))
logger.info("reference space has {} DoFs".format(reference_dg_space.num_DoFs))
logger.info("assembling Hdiv product ...")

rtn_space = RaviartThomasSpace(grid, 0)
FluxVectorSpace = DuneXTVectorSpace(IstlDenseVectorDouble, rtn_space.num_DoFs, "FLUX")
rtn_product = make_marix_operator(assemble_Hdiv_product_matrix(rtn_space), "FLUX")


def rtn_norm(t):
    if not isinstance(t, ListVectorArray):
        t = FluxVectorSpace.from_data(
            [
                t,
            ]
        )
    return np.sqrt(rtn_product.apply2(t, t)[0][0])


logger.info("computing ESV2007 pressure and flux ...")

u_h_f = make_marix_operator(
    assemble_SWIPDG_matrix(dg_space, ConstantFunction(1.0)), "PRESSURE"
).apply_inverse(fom.rhs.as_source_array())
t_h_f = FluxVectorSpace.from_data(
    [
        compute_flux_reconstruction(
            grid, dg_space, rtn_space, ConstantFunction(1.0), u_h_f._list[0].impl
        ),
    ]
)

logger.info("computing [OS2015, table 1] estimates (should be 0.166, 0.723, 0.355) ...")

_, eta_NC, eta_R, eta_DF = compute_estimate(
    grid,
    make_discrete_function(dg_space, u_h_f._list[0].impl, "u_h_f"),
    make_discrete_function(rtn_space, t_h_f._list[0].impl, "t_h_f"),
    f,
    ConstantFunction(1.0),
    ConstantFunction(1.0),
    ConstantFunction(1.0),
    1,
    1,
    1,
)

logger.info("    are {}, {}, {}".format(eta_NC, eta_R, eta_DF))
logger.info(
    "computing other OS2015 estimates (should be "
    "[table 3, eta_NC] 0.182, <= [table 1, eta_R] 0.166, [table 2, eta_DF] 0.316) ..."
)

mu = {"switch": 1}
diffusion_mu = ExpressionFunction(
    "x", [diffusion_expression.format(mu["switch"])], 3, "diffusion_mu"
)
u_h = fom.solve(mu)
t_h = FluxVectorSpace.from_data(
    [
        compute_flux_reconstruction(
            grid, dg_space, rtn_space, diffusion_mu, u_h._list[0].impl
        ),
    ]
)
_, eta_NC, eta_R, eta_DF = compute_estimate(
    grid,
    make_discrete_function(dg_space, u_h._list[0].impl, "u_h"),
    make_discrete_function(rtn_space, t_h._list[0].impl, "t_h"),
    f,
    diffusion_mu,
    diffusion_bar,
    diffusion_hat,
    alpha(mu, mu_bar),
    alpha(mu, mu_hat),
    gamma(mu, mu_bar),
)

logger.info("    are {}, {}, {}".format(eta_NC, eta_R, eta_DF))

from gamm_2019_talk_on_conservative_rb_base import simulate_single_greedy_step

RB_size = 0
for nn in range(1, 100):

    logger.info("simulating greedy step {} ...".format(nn))

    greedy_data, estimate_data = simulate_single_greedy_step(
        fom,
        dg_product=fom.energy_penalty_product,
        FluxVectorSpace=FluxVectorSpace,
        rtn_product=rtn_product,
        t_h_f=t_h_f,
        compute_flux_reconstruction=lambda mu, u_RB: compute_flux_reconstruction(
            grid,
            dg_space,
            rtn_space,
            ExpressionFunction(
                "x", [diffusion_expression.format(mu["switch"])], 3, "diffusion_mu"
            ),
            u_RB,
        ),
        compute_estimate=lambda mu, u_RB, t_RB_f: compute_estimate(
            grid,
            make_discrete_function(dg_space, u_RB, "u_RB"),
            make_discrete_function(rtn_space, t_RB_f, "t_RB_f"),
            f,
            ExpressionFunction(
                "x", [diffusion_expression.format(mu["switch"])], 3, "diffusion_mu"
            ),
            diffusion_bar,
            diffusion_hat,
            alpha(mu, mu_bar),
            alpha(mu, mu_hat),
            gamma(mu, mu_bar),
        ),
        compute_reference_error=lambda mu, u_RB: reference_dg_norm(
            reference_fom.solve(mu)._list[0].impl
            - prolong(dg_space, u_RB, reference_dg_space)
        ),
        max_extensions=nn,
        num_samples=10,
    )

    if greedy_data["extensions"] > RB_size:
        RB_size = greedy_data["extensions"]
    else:
        logger.info("  finished")
        break

    logger.info("  max greedy error: {}".format(greedy_data["max_errs"][-1]))
    logger.info("  worst error:      {}".format(np.max(estimate_data["errors"])))
    logger.info("  worst estimate:   {}".format(np.max(estimate_data["etas"])))
    logger.info("  worst efficiency: {}".format(np.max(estimate_data["efficiencies"])))
