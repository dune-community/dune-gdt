// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_DEFAULT_HH
#define DUNE_GDT_OPERATORS_DEFAULT_HH

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/grid/walker/apply-on.hh>
#include <dune/stuff/grid/walker.hh>
#include <dune/stuff/la/container/pattern.hh>

#include <dune/gdt/assembler/local/codim0.hh>
#include <dune/gdt/assembler/wrapper.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/localoperator/interfaces.hh>
#include <dune/gdt/spaces/interface.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * \todo Check parallel case: there is probably/definitely communication missing in apply2!
 */
template< class GridViewImp, class RangeImp, class SourceImp = RangeImp, class FieldImp = typename RangeImp::RangeFieldType >
class LocalizableProductDefault
  : public Stuff::Grid::Walker< GridViewImp >
{
  typedef Stuff::Grid::Walker< GridViewImp > BaseType;
public:
  using typename BaseType::GridViewType;
  using typename BaseType::EntityType;
  typedef RangeImp         RangeType;
  typedef SourceImp        SourceType;
  typedef FieldImp         FieldType;
private:
  static_assert(Stuff::is_localizable_function< SourceType >::value,
                "SourceType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(Stuff::is_localizable_function< RangeType >::value,
                "RangeType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(std::is_same< typename SourceType::EntityType, EntityType >::value,
                "The EntityType of SourceType and GridViewType have to match!");
  static_assert(std::is_same< typename RangeType::EntityType, EntityType >::value,
                "The EntityType of RangeType and GridViewType have to match!");
  static_assert(std::is_same< typename SourceType::DomainFieldType, typename GridViewType::ctype >::value,
                "The DomainFieldType of SourceType and GridViewType have to match!");
  static_assert(std::is_same< typename RangeType::DomainFieldType, typename GridViewType::ctype >::value,
                "The DomainFieldType of RangeType and GridViewType have to match!");
  static_assert(SourceType::dimDomain == GridViewType::dimension,
                "The dimDomain of SourceType and GridViewType have to match!");
  static_assert(RangeType::dimDomain  == GridViewType::dimension,
                "The dimDomain of RangeType and GridViewType have to match!");

public:
  LocalizableProductDefault(GridViewType grid_view, const RangeType& range, const SourceType& source)
    : BaseType(grid_view)
    , range_(range)
    , source_(source)
    , walked_(false)
  {}

  LocalizableProductDefault(GridViewType grid_view, const RangeType& range)
    : BaseType(grid_view)
    , range_(range)
    , source_(range)
    , walked_(false)
  {}

  const SourceType& source() const
  {
    return source_;
  }

  const RangeType& range() const
  {
    return range_;
  }

  template< class V >
  void add(const LocalVolumeTwoFormInterface< V >& local_volume_twoform,
           const DSG::ApplyOn::WhichEntity< GridViewType >* where = new DSG::ApplyOn::AllEntities< GridViewType >())
  {
    typedef LocalVolumeTwoFormAccumulator
        < GridViewType, typename LocalVolumeTwoFormInterface< V >::derived_type, RangeType, SourceType, FieldType > AccumulateFunctor;
    local_volume_twoforms_.emplace_back(
          new AccumulateFunctor(grid_view_, local_volume_twoform.as_imp(), range_, source_, *where));
    BaseType::add(*local_volume_twoforms_.back(), where);
  }

  FieldType compute_locally(const EntityType& entity) const
  {
    FieldType local_result = 0.;
    for (const auto& local_volume_twoform : local_volume_twoforms_)
      local_result += local_volume_twoform->compute_locally(entity);
    return local_result;
  }

  FieldType apply2()
  {
    if (!walked_) {
      this->walk();
      walked_ = true;
    }
    FieldType result = 0.;
    for (const auto& local_volume_twoform : local_volume_twoforms_)
      result += local_volume_twoform->result();
    return result;
  }

private:
  using BaseType::grid_view_;

  const RangeType& range_;
  const SourceType& source_;
  std::vector< std::unique_ptr< DSG::internal::Codim0ReturnObject< GridViewType, FieldType > > > local_volume_twoforms_;
  bool walked_;
}; // class LocalizableProductDefault


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_DEFAULT_HH
