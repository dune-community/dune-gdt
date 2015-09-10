// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_PROJECTIONS_L2_HH
#define DUNE_GDT_OPERATORS_PROJECTIONS_L2_HH

#include <dune/stuff/common/timedlogging.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/grid/entity.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/localoperator/l2-projection.hh>
#include <dune/gdt/operators/default.hh>
#include <dune/gdt/operators/interfaces.hh>
#include <dune/gdt/spaces/dg/interface.hh>

#include "../default.hh"

namespace Dune {
namespace GDT {


// forward
template< class GridViewImp, class FieldImp = double >
class L2LocalProjectionOperator;


namespace internal {


template< class GridViewImp, class FieldImp >
class L2LocalProjectionOperatorTraits
{
public:
  typedef L2LocalProjectionOperator< GridViewImp, FieldImp > derived_type;
  typedef FieldImp                                           FieldType;
};


} // namespace internal


template< class GridViewImp, class SourceImp, class RangeImp >
class L2LocalProjectionLocalizableOperator
  : public LocalizableOperatorDefault< GridViewImp, SourceImp, RangeImp >
{
  typedef LocalizableOperatorDefault< GridViewImp, SourceImp, RangeImp > BaseType;
  typedef LocalL2ProjectionOperator LocalOperatorType;
public:
  using typename BaseType::GridViewType;
  using typename BaseType::SourceType;
  using typename BaseType::RangeType;

  template< class ...Args >
  explicit L2LocalProjectionLocalizableOperator(const size_t over_integrate, Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
    , local_operator_(over_integrate)
  {
    this->add(local_operator_);
    issue_warning(this->range().space());
  }

  template< class ...Args >
  explicit L2LocalProjectionLocalizableOperator(Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
    , local_operator_()
  {
    this->add(local_operator_);
    issue_warning(this->range().space());
  }

private:
  template< class S, bool warn >
  struct Warning { static inline void issue() {} };

  template< class S >
  struct Warning< S, true >
  {
    static inline void issue()
    {
      DSC::TimedLogger().get("gdt.l2localprojectionlocalizableoperator").warn()
          << "You are using this operator to project onto a continuous discrete function space (see below)!\n"
          << "You can disable this warning by defining "
          << "DUNE_GDT_OPERATORS_PROJECTIONS_L2_LOCAL_LOCALIZABLE_DISABLE_WARNING\n"
          << "at compile time or by disabling the Dune::Stuff::Common::TimedLogger() instance at runtime.\n"
          << "The type of the range space is: " << DSC::Typename< S >::value() << std::endl;
    } // ... issue_warning(...)
  };

  template< class S >
  static inline void issue_warning(const S&)
  {
#ifndef DUNE_GDT_OPERATORS_PROJECTIONS_L2_LOCAL_DISABLE_WARNING
    Warning< S, S::continuous >::issue();
#endif
  }

  const LocalOperatorType local_operator_;
}; // class L2LocalProjectionLocalizableOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_PROJECTIONS_L2_HH
