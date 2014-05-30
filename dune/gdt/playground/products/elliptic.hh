// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PLAYGROUND_PRODUCTS_ELLIPTIC_HH
#define DUNE_GDT_PLAYGROUND_PRODUCTS_ELLIPTIC_HH

#include <dune/gdt/playground/localevaluation/elliptic.hh>
#include <dune/gdt/products/elliptic.hh>

namespace Dune {
namespace GDT {
namespace Products {


namespace internal {


template <class DiffusionFactorType, class GridViewImp, class FieldImp, class DiffusionTensorType>
class EllipticBaseTraits
{
public:
  typedef GridViewImp GridViewType;
  typedef FieldImp FieldType;
  typedef LocalOperator::Codim0Integral<LocalEvaluation::Elliptic<DiffusionFactorType, DiffusionTensorType>>
      LocalOperatorType;
}; // class EllipticBaseTraits


template <class DiffusionFactorType, class GridViewImp, class FieldImp, class DiffusionTensorType>
class EllipticBase
{
  typedef EllipticBaseTraits<DiffusionFactorType, GridViewImp, FieldImp, DiffusionTensorType> Traits;
  typedef typename Traits::LocalOperatorType LocalOperatorType;

public:
  EllipticBase(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
               const size_t over_integrate = 0)
    : diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , local_operator_(over_integrate, diffusion_factor_, diffusion_tensor_)
  {
  }

private:
  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;

protected:
  const LocalOperatorType local_operator_;
}; // class EllipticBase


template <class DiffusionFactorType, class GridViewImp, class RangeImp, class SourceImp, class FieldImp,
          class DiffusionTensorType>
class EllipticLocalizableTraits
    : public EllipticBaseTraits<DiffusionFactorType, GridViewImp, FieldImp, DiffusionTensorType>
{
public:
  typedef EllipticLocalizable<DiffusionFactorType, GridViewImp, RangeImp, SourceImp, FieldImp, DiffusionTensorType>
      derived_type;
  typedef RangeImp RangeType;
  typedef SourceImp SourceType;
  typedef FieldImp FieldType;
};


template <class DiffusionFactorType, class MatrixImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp,
          class DiffusionTensorType>
class EllipticAssemblableTraits
    : public EllipticBaseTraits<DiffusionFactorType, GridViewImp, typename MatrixImp::ScalarType, DiffusionTensorType>
{
public:
  typedef EllipticAssemblable<DiffusionFactorType, MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp,
                              DiffusionTensorType> derived_type;
  typedef RangeSpaceImp RangeSpaceType;
  typedef SourceSpaceImp SourceSpaceType;
  typedef MatrixImp MatrixType;
}; // class EllipticAssemblableTraits


template <class DiffusionFactorType, class GridViewImp, class FieldImp, class DiffusionTensorType>
class EllipticTraits
{
public:
  typedef Elliptic<DiffusionFactorType, GridViewImp, FieldImp, DiffusionTensorType> derived_type;
  typedef GridViewImp GridViewType;
  typedef FieldImp FieldType;
};


} // namespace internal


template <class DiffusionFactorType, class GridViewImp, class RangeImp, class SourceImp, class FieldImp,
          class DiffusionTensorType>
class EllipticLocalizable
    : public Products::LocalizableBase<internal::EllipticLocalizableTraits<DiffusionFactorType, GridViewImp, RangeImp,
                                                                           SourceImp, FieldImp, DiffusionTensorType>>,
      public internal::EllipticBase<DiffusionFactorType, GridViewImp, FieldImp, DiffusionTensorType>
{
  typedef Products::LocalizableBase<internal::EllipticLocalizableTraits<DiffusionFactorType, GridViewImp, RangeImp,
                                                                        SourceImp, FieldImp, DiffusionTensorType>>
      LocalizableBaseType;
  typedef internal::EllipticBase<DiffusionFactorType, GridViewImp, FieldImp, DiffusionTensorType> EllipticBaseType;

public:
  typedef internal::EllipticLocalizableTraits<DiffusionFactorType, GridViewImp, RangeImp, SourceImp, FieldImp,
                                              DiffusionTensorType> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::SourceType SourceType;

private:
  typedef typename Traits::LocalOperatorType LocalOperatorType;

public:
  EllipticLocalizable(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
                      const GridViewType& grd_vw, const RangeType& rng, const SourceType& src,
                      const size_t over_integrate = 0)
    : LocalizableBaseType(grd_vw, rng, src)
    , EllipticBaseType(diffusion_factor, diffusion_tensor, over_integrate)
  {
  }

  EllipticLocalizable(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
                      const GridViewType& grd_vw, const RangeType& rng, const size_t over_integrate = 0)
    : LocalizableBaseType(grd_vw, rng)
    , EllipticBaseType(diffusion_factor, diffusion_tensor, over_integrate)
  {
  }

private:
  virtual const LocalOperatorType& local_operator() const DS_OVERRIDE DS_FINAL
  {
    return this->local_operator_;
  }
}; // class EllipticLocalizable


template <class DiffusionFactorType, class MatrixImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp,
          class DiffusionTensorType>
class EllipticAssemblable
    : public Products::AssemblableBase<internal::EllipticAssemblableTraits<DiffusionFactorType, MatrixImp,
                                                                           RangeSpaceImp, GridViewImp, SourceSpaceImp,
                                                                           DiffusionTensorType>>,
      public internal::EllipticBase<DiffusionFactorType, GridViewImp, typename MatrixImp::ScalarType,
                                    DiffusionTensorType>
{
  typedef Products::AssemblableBase<internal::EllipticAssemblableTraits<DiffusionFactorType, MatrixImp, RangeSpaceImp,
                                                                        GridViewImp, SourceSpaceImp,
                                                                        DiffusionTensorType>> AssemblableBaseType;
  typedef internal::EllipticBase<DiffusionFactorType, GridViewImp, typename MatrixImp::ScalarType, DiffusionTensorType>
      EllipticBaseType;

public:
  typedef internal::EllipticAssemblableTraits<DiffusionFactorType, MatrixImp, RangeSpaceImp, GridViewImp,
                                              SourceSpaceImp, DiffusionTensorType> Traits;
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

  EllipticAssemblable(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
                      MatrixType& matrix, const RangeSpaceType& range_space, const GridViewType& grid_view,
                      const SourceSpaceType& source_space)
    : AssemblableBaseType(matrix, range_space, grid_view, source_space)
    , EllipticBaseType(diffusion_factor, diffusion_tensor)
  {
  }

  EllipticAssemblable(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
                      MatrixType& matrix, const RangeSpaceType& range_space, const GridViewType& grid_view)
    : AssemblableBaseType(matrix, range_space, grid_view, range_space)
    , EllipticBaseType(diffusion_factor, diffusion_tensor)
  {
  }

  EllipticAssemblable(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
                      MatrixType& matrix, const RangeSpaceType& range_space)
    : AssemblableBaseType(matrix, range_space, *(range_space.grid_view()), range_space)
    , EllipticBaseType(diffusion_factor, diffusion_tensor)
  {
  }

private:
  virtual const LocalOperatorType& local_operator() const DS_OVERRIDE DS_FINAL
  {
    return this->local_operator_;
  }
}; // class EllipticAssemblable


template <class DiffusionFactorType, class GridViewImp, class FieldImp, class DiffusionTensorType>
class Elliptic
    : public ProductInterface<internal::EllipticTraits<DiffusionFactorType, GridViewImp, FieldImp, DiffusionTensorType>>
{
public:
  typedef internal::EllipticTraits<DiffusionFactorType, GridViewImp, FieldImp, DiffusionTensorType> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::FieldType FieldType;

  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;

  Elliptic(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
           const GridViewType& grd_vw)
    : diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , grid_view_(grd_vw)
  {
  }

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  template <class RR, int rRR, int rCR, class RS, int rRS, int rCS>
  FieldType
  apply2(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RR, rRR, rCR>& /*range*/,
         const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RS, rRS, rCS>& /*source*/,
         const size_t /*over_integrate*/ = 0) const
  {
    static_assert((Dune::AlwaysFalse<RR>::value), "Not implemented for this combination!");
  }

  template <int dimRangeRows, int dimRangeCols>
  FieldType apply2(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType,
                                                             dimRangeRows, dimRangeCols>& range,
                   const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType,
                                                             dimRangeRows, dimRangeCols>& source,
                   const size_t over_integrate = 0) const
  {
    typedef Stuff::
        LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType, dimRangeRows, dimRangeCols>
            FunctionType;
    EllipticLocalizable<DiffusionFactorType, GridViewType, FunctionType, FunctionType, FieldType, DiffusionTensorType>
        product_operator(diffusion_factor_, diffusion_tensor_, grid_view_, range, source, over_integrate);
    return product_operator.apply2();
  } // ... apply2(...)

private:
  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
  const GridViewType& grid_view_;
}; // class Elliptic


} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_PRODUCTS_ELLIPTIC_HH
