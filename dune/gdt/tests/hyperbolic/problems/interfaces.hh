// This file is part of the dune-hdd project:
//   http://users.dune-project.org/projects/dune-hdd
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_INTERFACES_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_INTERFACES_HH

#include <ostream>

#include <dune/grid/yaspgrid.hh>

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/functions/default.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/functions/checkerboard.hh>

namespace Dune {
namespace GDT {
namespace Hyperbolic {


/* Interface for problem of the form delta_t u + div f(u) = q(u) where u: R^d \to R.
 * The flux f is a function f: R \to R^d, and q: R \to R is a source. */
template< class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim >
class ProblemInterface
{
  typedef ProblemInterface< EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim > ThisType;
public:
  typedef EntityImp         EntityType;
  typedef DomainFieldImp    DomainFieldType;
  static const size_t dimDomain = domainDim;
  typedef RangeFieldImp     RangeFieldType;
  static const size_t dimRange = rangeDim;

//  we do not have a grid for the range space, but still need an EntityType with dimension dimRange for FluxType and
//  SourceType, so just choose an arbitrary one
  typedef typename Dune::YaspGrid< dimRange >::template Codim< 0 >::Entity           FluxSourceEntityType;
  typedef Dune::Stuff::GlobalFunctionInterface< FluxSourceEntityType,
                                                RangeFieldType, dimRange,
                                                RangeFieldType, dimRange, dimDomain >                   FluxType;
  typedef Dune::Stuff::GlobalFunctionValuedFunctionInterface< EntityType, DomainFieldType, dimDomain,
                                                              FluxSourceEntityType, RangeFieldType, dimRange,
                                                              RangeFieldType, dimRange, 1 >             SourceType;
  typedef Dune::Stuff::GlobalFunctionValuedFunctionInterface< EntityType, DomainFieldType, dimDomain,
                                                              EntityType, DomainFieldType, dimDomain,
                                                              RangeFieldType, dimRange, 1 >             FunctionType;
  typedef typename Dune::Stuff::Functions::TimeDependentExpression
                < EntityImp, DomainFieldImp, dimDomain, RangeFieldImp, dimRange, 1, double >            BoundaryValueType;
  typedef Dune::Stuff::Common::Configuration                                                            ConfigType;
  typedef DS::TimeDependentFunctionInterface
                  < typename DS::LocalizableFunctionInterface < EntityType,
                                                                DomainFieldType, dimDomain,
                                                                RangeFieldType, dimRange, 1 > >         SolutionType;

  virtual ~ProblemInterface() {}

  static std::string static_id()
  {
    return "hdd.hyperbolic.problem";
  }

  virtual std::string type() const
  {
    return "hdd.hyperbolic.problem";
  }

  virtual const std::shared_ptr< const FluxType >& flux() const = 0;

  virtual const std::shared_ptr< const SourceType >& source() const = 0;

  virtual const std::shared_ptr< const FunctionType >& initial_values() const = 0;

  virtual const ConfigType grid_config() const = 0;

  virtual const ConfigType boundary_info() const = 0;

  virtual const std::shared_ptr< const BoundaryValueType >& boundary_values() const = 0;

  virtual double t_end() const = 0;

  virtual double CFL() const = 0;

  virtual bool is_linear() const
  {
    return false;
  }

  virtual void report(std::ostream& out, std::string prefix = "") const
  {
    out << prefix << "problem '" << type() << "':" << std::endl;
    // TODO: implement
  } // ... report(...)

private:
  template< class T >
  friend std::ostream& operator<<(std::ostream& /*out*/, const ThisType& /*problem*/);
}; // ProblemInterface


template< class E, class D, int d, class R, int r >
std::ostream& operator<<(std::ostream& out, const ProblemInterface< E, D, d, R, r >& problem)
{
  problem.report(out);
  return out;
} // ... operator<<(...)


} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_INTERFACES_HH
