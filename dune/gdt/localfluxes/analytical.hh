#ifndef DUNE_GDT_LOCALFLUXES_ANALYTICAL_HH
#define DUNE_GDT_LOCALFLUXES_ANALYTICAL_HH

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
    return "gdt.globalfunctionbasedanalyticalflux";
  }

  static std::unique_ptr<ThisType> create(const DSC::Configuration global_function_config,
                                          const std::string sub_name = static_id())
  {
    return DSC::make_unique<ThisType>(*(GlobalFunctionType::create(global_function_config, sub_name)));
  } // ... create(...)
private:
  const GlobalFunctionType global_function_;
};

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALFLUXES_ANALYTICAL_HH
