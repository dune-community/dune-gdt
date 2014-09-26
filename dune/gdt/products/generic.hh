// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_GENERIC_HH
#define DUNE_GDT_PRODUCTS_GENERIC_HH

#include <type_traits>

#include <dune/stuff/common/disable_warnings.hh>
#include <dune/grid/common/gridview.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/common/crtp.hh>

#include <dune/gdt/products/l2_internal.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/products/base.hh>

namespace Dune {
namespace GDT {
namespace Products {

template <class GridViewImp, class RangeImp, class SourceImp,
          template <class, class, class, class, template <class> class> class TraitsTemplate,
          template <class> class LocalEvaluationTemplate>
class GenericLocalizable
    : public LocalizableBase<TraitsTemplate<GridViewImp, RangeImp, SourceImp,
                                            GenericLocalizable<GridViewImp, RangeImp, SourceImp, TraitsTemplate,
                                                               LocalEvaluationTemplate>,
                                            LocalEvaluationTemplate>>,
      public internal::L2Base<GridViewImp, typename RangeImp::RangeFieldType, LocalEvaluationTemplate>
{
  typedef GenericLocalizable<GridViewImp, RangeImp, SourceImp, TraitsTemplate, LocalEvaluationTemplate> ThisType;

public:
  typedef TraitsTemplate<GridViewImp, RangeImp, SourceImp, ThisType, LocalEvaluationTemplate> Traits;

private:
  typedef Products::LocalizableBase<Traits> LocalizableBaseType;
  typedef internal::L2Base<GridViewImp, typename RangeImp::RangeFieldType, LocalEvaluationTemplate> L2BaseType;

public:
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::SourceType SourceType;

private:
  typedef typename Traits::LocalOperatorType LocalOperatorType;

public:
  GenericLocalizable(const GridViewType& grd_vw, const RangeType& rng, const SourceType& src,
                     const size_t over_integrate = 0)
    : LocalizableBaseType(grd_vw, rng, src)
    , L2BaseType(over_integrate)
  {
  }

  GenericLocalizable(const GridViewType& grd_vw, const RangeType& rng, const size_t over_integrate = 0)
    : LocalizableBaseType(grd_vw, rng, rng)
    , L2BaseType(over_integrate)
  {
  }

private:
  virtual const LocalOperatorType& local_operator() const DS_OVERRIDE DS_FINAL
  {
    return this->local_operator_;
  }
}; // class L2Localizable

template <class MatrixImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp,
          template <class, class, class, class, class, template <class> class> class TraitsTemplate,
          template <class> class LocalEvaluationTemplate>
class GenericAssemblable
    : public AssemblableBase<TraitsTemplate<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp,
                                            GenericAssemblable<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp,
                                                               TraitsTemplate, LocalEvaluationTemplate>,
                                            LocalEvaluationTemplate>>,
      public internal::L2Base<GridViewImp, typename RangeSpaceImp::RangeFieldType, LocalEvaluationTemplate>
{
  typedef GenericAssemblable<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp, TraitsTemplate,
                             LocalEvaluationTemplate> ThisType;

public:
  typedef TraitsTemplate<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp, ThisType, LocalEvaluationTemplate>
      Traits;

private:
  typedef Products::AssemblableBase<Traits> AssemblableBaseType;
  typedef internal::L2Base<GridViewImp, typename RangeSpaceImp::RangeFieldType, LocalEvaluationTemplate> L2BaseType;

public:
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::RangeSpaceType RangeSpaceType;
  typedef typename Traits::SourceSpaceType SourceSpaceType;
  typedef typename Traits::MatrixType MatrixType;

private:
  typedef typename Traits::LocalOperatorType LocalOperatorType;

public:
  using AssemblableBaseType::pattern;

  static Stuff::LA::SparsityPatternDefault pattern(const RangeSpaceType& range_space,
                                                   const SourceSpaceType& source_space, const GridViewType& grid_view)
  {
    return range_space.compute_volume_pattern(grid_view, source_space);
  }

  GenericAssemblable(MatrixType& mtrx, const RangeSpaceType& rng_scp, const GridViewType& grd_vw,
                     const SourceSpaceType& src_scp, const size_t over_integrate = 0)
    : AssemblableBaseType(mtrx, rng_scp, grd_vw, src_scp)
    , L2BaseType(over_integrate)
  {
  }

  GenericAssemblable(MatrixType& mtrx, const RangeSpaceType& rng_scp, const GridViewType& grd_vw,
                     const size_t over_integrate = 0)
    : AssemblableBaseType(mtrx, rng_scp, grd_vw, rng_scp)
    , L2BaseType(over_integrate)
  {
  }

  GenericAssemblable(MatrixType& matrix, const RangeSpaceType& range_space, const size_t over_integrate = 0)
    : AssemblableBaseType(matrix, range_space, *(range_space.grid_view()), range_space)
    , L2BaseType(over_integrate)
  {
  }

private:
  virtual const LocalOperatorType& local_operator() const DS_OVERRIDE DS_FINAL
  {
    return this->local_operator_;
  }
}; // class GenericAssemblable


template <class GridViewImp, class FieldImp, template <class, class, class> class OperatorTemplate>
class GenericProduct
    : public ProductInterface<internal::DerivedType<GridViewImp, FieldImp,
                                                    GenericProduct<GridViewImp, FieldImp, OperatorTemplate>>>
{
public:
  typedef internal::DerivedType<GridViewImp, FieldImp, GenericProduct<GridViewImp, FieldImp, OperatorTemplate>> Traits;
  typedef GridViewImp GridViewType;
  typedef FieldImp FieldType;

  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;

  GenericProduct(const GridViewType& grd_vw, const size_t over_integrate = 0)
    : grid_view_(grd_vw)
    , over_integrate_(over_integrate)
  {
  }

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  template <class RR, int rRR, int rCR, class RS, int rRS, int rCS>
  FieldType apply2(
      const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RR, rRR, rCR>& /*range*/,
      const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RS, rRS, rCS>& /*source*/) const
  {
    static_assert((Dune::AlwaysFalse<RR>::value), "Not implemented for this combination!");
  }

  template <int dimRangeRows, int dimRangeCols>
  FieldType apply2(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType,
                                                             dimRangeRows, dimRangeCols>& range,
                   const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType,
                                                             dimRangeRows, dimRangeCols>& source) const
  {
    typedef Stuff::
        LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType, dimRangeRows, dimRangeCols>
            RangeType;
    OperatorTemplate<GridViewType, RangeType, RangeType> product_operator(grid_view_, range, source, over_integrate_);
    return product_operator.apply2();
  } // ... apply2(...)

private:
  const GridViewType& grid_view_;
  const size_t over_integrate_;
}; // class L2


} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_GENERIC_HH
