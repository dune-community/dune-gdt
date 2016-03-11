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


template <class GridViewType, class FieldImp = double>
class L2ProlongationOperatorTraits
{
public:
  typedef L2ProlongationOperator<GridViewType, FieldImp> derived_type;
  typedef FieldImp FieldType;
};


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


template <class GridViewImp, class FieldImp>
class L2ProlongationOperator : public OperatorInterface<internal::L2ProlongationOperatorTraits<GridViewImp, FieldImp>>
{
  typedef OperatorInterface<internal::L2ProlongationOperatorTraits<GridViewImp, FieldImp>> BaseType;

public:
  typedef internal::L2ProlongationOperatorTraits<GridViewImp, FieldImp> Traits;
  typedef GridViewImp GridViewType;
  using typename BaseType::FieldType;

private:
  typedef typename Stuff::Grid::Entity<GridViewType>::Type E;
  typedef typename GridViewType::ctype D;
  static const size_t d = GridViewType::dimension;

public:
  L2ProlongationOperator(const size_t over_integrate, GridViewType grid_view)
    : grid_view_(grid_view)
    , over_integrate_(over_integrate)
  {
  }

  L2ProlongationOperator(GridViewType grid_view)
    : grid_view_(grid_view)
    , over_integrate_(0)
  {
  }

  template <class SS, class SV, class RS, class RV>
  void apply(const ConstDiscreteFunction<SS, SV>& source, DiscreteFunction<RS, RV>& range) const
  {
    redirect<RS::continuous>::apply(grid_view_, source, range, over_integrate_);
  }

  template <class RangeType, class SourceType>
  FieldType apply2(const RangeType& /*range*/, const SourceType& /*source*/) const
  {
    DUNE_THROW(NotImplemented, "Go ahead if you think this makes sense!");
  }

  template <class RangeType, class SourceType>
  void apply_inverse(const RangeType& /*range*/, SourceType& /*source*/,
                     const Stuff::Common::Configuration& /*opts*/) const
  {
    DUNE_THROW(NotImplemented, "Go ahead if you think this makes sense!");
  }

  std::vector<std::string> invert_options() const
  {
    DUNE_THROW(NotImplemented, "Go ahead if you think this makes sense!");
  }

  Stuff::Common::Configuration invert_options(const std::string& /*type*/) const
  {
    DUNE_THROW(NotImplemented, "Go ahead if you think this makes sense!");
  }

private:
  template <bool continuous = true, bool anything = true>
  struct redirect
  {
    template <class SourceType, class RangeType>
    static void apply(const GridViewType& grd_vw, const SourceType& src, RangeType& rng, const size_t over_integrate)
    {
      L2GlobalProlongationLocalizableOperator<GridViewType, SourceType, RangeType>(over_integrate, grd_vw, src, rng)
          .apply();
    }
  };

  template <bool anything>
  struct redirect<false, anything>
  {
    template <class SourceType, class RangeType>
    static void apply(const GridViewType& grd_vw, const SourceType& src, RangeType& rng, const size_t over_integrate)
    {
      L2LocalProlongationLocalizableOperator<GridViewType, SourceType, RangeType>(over_integrate, grd_vw, src, rng)
          .apply();
    }
  };

  GridViewType grid_view_;
  const size_t over_integrate_;
}; // class L2ProlongationOperator


template <class GridViewType>
typename std::enable_if<Stuff::Grid::is_grid_layer<GridViewType>::value,
                        std::unique_ptr<L2ProlongationOperator<GridViewType>>>::type
make_l2_prolongation_operator(const GridViewType& grid_view, const size_t over_integrate = 0)
{
  return DSC::make_unique<L2ProlongationOperator<GridViewType>>(over_integrate, grid_view);
}


template <class GridViewType, class SS, class SV, class RS, class RV>
typename std::enable_if<Stuff::Grid::is_grid_layer<GridViewType>::value, void>::type
prolong_l2(const GridViewType& grid_view, const ConstDiscreteFunction<SS, SV>& source, DiscreteFunction<RS, RV>& range,
           const size_t over_integrate = 0)
{
  make_l2_prolongation_operator(grid_view, over_integrate)->apply(source, range);
}

template <class SS, class SV, class RS, class RV>
void prolong_l2(const ConstDiscreteFunction<SS, SV>& source, DiscreteFunction<RS, RV>& range,
                const size_t over_integrate = 0)
{
  make_l2_prolongation_operator(range.space().grid_view(), over_integrate)->apply(source, range);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PROLONGATIONS_L2_HH
