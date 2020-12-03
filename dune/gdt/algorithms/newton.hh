// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef DUNE_GDT_ALGORITHMS_NEWTON_HH
#define DUNE_GDT_ALGORITHMS_NEWTON_HH

#include <dune/xt/common/parameter.hh>
#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/la/container/vector-interface.hh>
#include <dune/xt/la/solver.hh>

#include <dune/gdt/exceptions.hh>

namespace Dune {
namespace GDT {


// forward, no include
template <class AGV, size_t s_r, size_t s_rC, size_t r_r, size_t r_rC, class F, class M, class SGV, class RGV>
class OperatorInterface;


static inline XT::Common::Configuration default_newton_solve_options()
{
  return {{{"precision", "1e-7"}, {"max_iter", "100"}, {"max_dampening_iter", "1000"}}};
}


/**
 * \brief Computes the inverse action of the operator.
 *
 * Currently implemented is the dampened Newton from [DF2015, Sec. 8.4.4.1].
 *
 * \todo Allow to pass jacobian options as subcfg in newton.
 **/
template <class AGV, size_t s_r, size_t s_rC, size_t r_r, size_t r_rC, class F, class M, class SGV, class RGV, class V>
void newton_solve(const OperatorInterface<AGV, s_r, s_rC, r_r, r_rC, F, M, SGV, RGV>& lhs_operator,
                  const XT::LA::VectorInterface<V>& rhs_vector,
                  XT::LA::VectorInterface<V>& initial_guess_vector,
                  const XT::Common::Parameter& param = {},
                  const XT::Common::Configuration& opts = default_newton_solve_options())
{
  XT::Common::DefaultLogger logger("newton_solve");
  LOG(debug) << "(lhs_operator=" << &lhs_operator << ", rhs_vector.sup_norm()=" << rhs_vector.sup_norm()
             << ", initial_guess_vector.sup_norm()=" << initial_guess_vector.sup_norm() << ", param=" << print(param)
             << ", opts " << print(opts, {{"oneline", "true"}}) << ")" << std::endl;
  // some preparations
  const auto& rhs = rhs_vector.as_imp();
  auto& initial_guess = initial_guess_vector.as_imp();
  auto residual_op = lhs_operator - rhs;
  auto residual = rhs.copy();
  auto update = initial_guess.copy();
  auto candidate = initial_guess.copy();
  const auto default_opts = default_newton_solve_options();
  const auto precision = opts.get("precision", default_opts.get<double>("precision"));
  const auto max_iter = opts.get("max_iter", default_opts.get<size_t>("max_iter"));
  const auto max_dampening_iter = opts.get("max_dampening_iter", default_opts.get<size_t>("max_dampening_iter"));
  // one matrix for all jacobians
  auto jacobian_op = residual_op.empty_jacobian_op();
  auto jacobian_solver = XT::LA::make_solver(jacobian_op.matrix());
  LOG(info) << "solving system by dampened newton:" << std::endl;
  size_t l = 0;
  Timer timer;
  while (true) {
    const auto prefix = "      " + XT::Common::whitespaceify(l);
    timer.reset();
    LOG(info) << "l = " << l << ": computing residual ... " << std::flush;
    residual_op.apply(initial_guess, residual, param);
    auto res = residual.l2_norm();
    LOG(info) << "took " << timer.elapsed() << "s, |residual|_l2 = " << res << std::endl;
    if (res < precision) {
      LOG(info) << prefix << "residual below tolerance, succeeded!" << std::endl;
      break;
    }
    DUNE_THROW_IF(l >= max_iter,
                  Exceptions::newton_error,
                  "max iterations reached!\n|residual|_l2 = " << res << "\nopts:\n"
                                                              << opts);
    LOG(info) << prefix << "computing jacobi matrix ... " << std::flush;
    timer.reset();
    jacobian_op.matrix() *= 0.;
    residual_op.jacobian(initial_guess, jacobian_op, {{"type", residual_op.jacobian_options().at(0)}}, param);
    jacobian_op.assemble(/*use_tbb=*/true);
    LOG(info) << "took " << timer.elapsed() << "s\n" << prefix << "solving for defect ... " << std::flush;
    timer.reset();
    residual *= -1.;
    update = initial_guess; // <- initial guess for the linear solver
    bool linear_solve_succeeded = false;
    std::vector<std::string> tried_linear_solvers;
    for (const auto& linear_solver_type : jacobian_solver.types()) {
      try {
        tried_linear_solvers.push_back(linear_solver_type);
        jacobian_solver.apply(
            residual, update, {{"type", linear_solver_type}, {"precision", XT::Common::to_string(0.1 * precision)}});
        linear_solve_succeeded = true;
        break;
      } catch (const XT::LA::Exceptions::linear_solver_failed&) {
      }
    }
    DUNE_THROW_IF(!linear_solve_succeeded,
                  Exceptions::newton_error,
                  "could not solve linear system for defect!\nTried the following linear solvers: "
                      << tried_linear_solvers << "\nopts:\n"
                      << opts);
    LOG(info) << "took " << timer.elapsed() << "s";
    if (tried_linear_solvers.size() > 1) {
      LOG(info) << " (and " << tried_linear_solvers.size() << " attempts with different linear solvers)";
    }
    LOG(info) << "\n" << prefix << "computing update ... " << std::flush;
    timer.reset();
    // try the automatic dampening strategy proposed in [DF2015, Sec. 8.4.4.1, p. 432]
    size_t k = 0;
    auto candidate_res = 2 * res; // any number such that we enter the while loop at least once
    double lambda = 1;
    while (!(candidate_res / res < 1)) {
      DUNE_THROW_IF(k >= max_dampening_iter,
                    Exceptions::newton_error,
                    "max iterations reached when trying to compute automatic dampening!\n|residual|_l2 = "
                        << res << "\nl = " << l << "\nopts:\n"
                        << opts);
      candidate = initial_guess + update * lambda;
      residual_op.apply(candidate, residual, param);
      candidate_res = residual.l2_norm();
      lambda /= 2;
      k += 1;
    }
    initial_guess = candidate;
    LOG(info) << "took " << timer.elapsed() << "s and a dampening of " << 2 * lambda << std::endl;
    l += 1;
  }
} // ... newton (...)


} // namespace GDT
} // namespace Dune


#endif // DUNE_GDT_ALGORITHMS_NEWTON_HH
