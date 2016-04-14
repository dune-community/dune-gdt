// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_LOCAL_FLUXES_ANALYTICAL_HH
#define DUNE_GDT_LOCAL_FLUXES_ANALYTICAL_HH

#include <dune/stuff/common/configuration.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


template <class GlobalFunctionImp, class E, class D, size_t d, class R, size_t r, size_t rC>
class GlobalFunctionBasedAnalyticalFlux : public AutonomousAnalyticalFluxInterface<E, D, d, R, r, rC>
{
  typedef AutonomousAnalyticalFluxInterface<E, D, d, R, r, rC> BaseType;
  typedef GlobalFunctionBasedAnalyticalFlux<GlobalFunctionImp, E, D, d, R, r, rC> ThisType;

public:
  typedef GlobalFunctionImp GlobalFunctionType;
  using typename BaseType::FluxRangeType;
  using typename BaseType::FluxJacobianRangeType;
  using typename BaseType::RangeType;
  static_assert(std::is_same<RangeType, typename GlobalFunctionType::DomainType>::value, "");
  static_assert(std::is_same<FluxRangeType, typename GlobalFunctionType::RangeType>::value, "");
  static_assert(std::is_same<FluxJacobianRangeType, typename GlobalFunctionType::JacobianRangeType>::value, "");

  GlobalFunctionBasedAnalyticalFlux(const GlobalFunctionType& global_function)
    : global_function_(global_function)
  {
  }

  virtual FluxRangeType evaluate(const RangeType& u) const
  {
    return global_function_.evaluate(u);
  }

  virtual FluxJacobianRangeType jacobian(const RangeType& u) const
  {
    return global_function_.jacobian(u);
  }

  static std::string static_id()
  {
    return "gdt.GlobalFunctionBasedAnalyticalFluxflux";
  }

  static std::unique_ptr<ThisType> create(const DSC::Configuration global_function_config,
                                          const std::string sub_name = static_id())
  {
    return DSC::make_unique<ThisType>(*(GlobalFunctionType::create(global_function_config, sub_name)));
  } // ... create(...)
private:
  const GlobalFunctionType global_function_;
}; // class GlobalFunctionBasedAnalyticalFlux


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FLUXES_ANALYTICAL_HH
