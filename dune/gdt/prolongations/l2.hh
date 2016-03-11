#ifndef DUNE_GDT_PROLONGATIONS_L2_HH
#define DUNE_GDT_PROLONGATIONS_L2_HH

#include <dune/gdt/discretefunction/default.hh>

#include "l2-local.hh"
#include "l2-global.hh"

namespace Dune {
namespace GDT {

// forward
template <class GridViewImp, class FieldImp = double>
class L2ProlongationOperator;


namespace internal {


template <class GridViewType, class SourceType, class RangeType>
class L2ProlongationLocalizableOperatorTraits
{
  static_assert(is_const_discrete_function<SourceType>::value, "");
  static_assert(is_discrete_function<RangeType>::value, "");

  template <class G, class S, class R, bool c = true>
  struct Helper
  {
    typedef L2GlobalProlongationLocalizableOperator<G, S, R> type;
  };

  template <class G, class S, class R>
  struct Helper<G, S, R, false>
  {
    typedef L2LocalProlongationLocalizableOperator<G, S, R> type;
  };

public:
  typedef typename Helper<GridViewType, SourceType, RangeType, RangeType::SpaceType::continuous>::type BaseType;
}; // class L2ProlongationLocalizableOperatorTraits


} // namespace internal


template <class GridViewImp, class SourceImp, class RangeImp>
class L2ProlongationLocalizableOperator
    : public internal::L2ProlongationLocalizableOperatorTraits<GridViewImp, SourceImp, RangeImp>::BaseType
{
  typedef
      typename internal::L2ProlongationLocalizableOperatorTraits<GridViewImp, SourceImp, RangeImp>::BaseType BaseType;

public:
  template <class... Args>
  explicit L2ProlongationLocalizableOperator(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
};


template <class GridViewType, class SS, class SV, class RS, class RV>
typename std::enable_if<Stuff::Grid::is_grid_layer<GridViewType>::value,
                        std::unique_ptr<L2ProlongationLocalizableOperator<GridViewType, ConstDiscreteFunction<SS, SV>,
                                                                          DiscreteFunction<RS, RV>>>>::type
make_l2_prolongation_localizable_operator(const GridViewType& grid_view, const ConstDiscreteFunction<SS, SV>& source,
                                          DiscreteFunction<RS, RV>& range, const size_t over_integrate = 0)
{
  return DSC::make_unique<L2ProlongationLocalizableOperator<GridViewType,
                                                            ConstDiscreteFunction<SS, SV>,
                                                            DiscreteFunction<RS, RV>>>(
      over_integrate, grid_view, source, range);
} // ... make_l2_prolongation_localizable_operator(...)

template <class SS, class SV, class RS, class RV>
std::unique_ptr<L2ProlongationLocalizableOperator<typename RS::GridViewType, ConstDiscreteFunction<SS, SV>,
                                                  DiscreteFunction<RS, RV>>>
make_l2_prolongation_localizable_operator(const ConstDiscreteFunction<SS, SV>& source, DiscreteFunction<RS, RV>& range,
                                          const size_t over_integrate = 0)
{
  return DSC::make_unique<L2ProlongationLocalizableOperator<typename RS::GridViewType,
                                                            ConstDiscreteFunction<SS, SV>,
                                                            DiscreteFunction<RS, RV>>>(
      over_integrate, range.space().grid_view(), source, range);
} // ... make_l2_prolongation_localizable_operator(...)


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PROLONGATIONS_L2_HH
