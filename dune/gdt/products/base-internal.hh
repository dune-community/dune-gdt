// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_BASE_INTERNAL_HH
#define DUNE_GDT_PRODUCTS_BASE_INTERNAL_HH

#include <type_traits>

#include <dune/grid/common/gridview.hh>

#include <dune/stuff/functions/interfaces.hh>

#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {
namespace Products {


// forwards
template <class LocalOperatorProvider, class RangeImp, class SourceImp>
class LocalizableBase;


template <class LocalOperatorProvider, class MatrixImp, class RangeSpaceImp, class SourceSpaceImp>
class AssemblableBase;


template <class LocalOperatorProvider>
class GenericBase;


namespace internal {


/**
 * \note not of interest to the average user
 */
template <class LocalOperatorProvider, class RangeImp, class SourceImp>
class LocalizableBaseTraits
{
public:
  typedef LocalizableBase<LocalOperatorProvider, RangeImp, SourceImp> derived_type;
  typedef typename LocalOperatorProvider::GridViewType GridViewType;
  typedef typename LocalOperatorProvider::FieldType FieldType;
  typedef RangeImp RangeType;
  typedef SourceImp SourceType;

private:
  static_assert(std::is_base_of<Dune::GridView<typename GridViewType::Traits>, GridViewType>::value,
                "GridViewType has to be derived from GridView!");
  static_assert(std::is_base_of<Stuff::Tags::LocalizableFunction, RangeType>::value,
                "RangeType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(std::is_base_of<Stuff::Tags::LocalizableFunction, SourceType>::value,
                "SourceType has to be derived from Stuff::LocalizableFunctionInterface!");
}; // class LocalizableBaseTraits


/**
 * \note not of interest to the average user
 */
template <class LocalOperatorProvider, class MatrixImp, class RangeSpaceImp, class SourceSpaceImp>
class AssemblableBaseTraits
{
public:
  typedef AssemblableBase<LocalOperatorProvider, MatrixImp, RangeSpaceImp, SourceSpaceImp> derived_type;
  typedef typename LocalOperatorProvider::GridViewType GridViewType;
  typedef RangeSpaceImp RangeSpaceType;
  typedef SourceSpaceImp SourceSpaceType;
  typedef MatrixImp MatrixType;

private:
  static_assert(std::is_base_of<Dune::GridView<typename GridViewType::Traits>, GridViewType>::value,
                "GridViewType has to be derived from GridView!");
  static_assert(std::is_base_of<SpaceInterface<typename RangeSpaceType::Traits>, RangeSpaceType>::value,
                "RangeSpaceType has to be derived from SpaceInterface!");
  static_assert(std::is_base_of<SpaceInterface<typename SourceSpaceType::Traits>, SourceSpaceType>::value,
                "SourceSpaceType has to be derived from SpaceInterface!");
  static_assert(std::is_base_of<Stuff::LA::MatrixInterface<typename MatrixType::Traits>, MatrixType>::value,
                "MatrixType has to be derived from Stuff::LA::MatrixInterface!");
}; // class AssemblableBaseTraits


/**
 * \note not of interest to the average user
 */
template <class LocalOperatorProvider>
class GenericBaseTraits
{
public:
  typedef GenericBase<LocalOperatorProvider> derived_type;
  typedef typename LocalOperatorProvider::FieldType FieldType;
  typedef typename LocalOperatorProvider::GridViewType GridViewType;

private:
  static_assert(std::is_base_of<Dune::GridView<typename GridViewType::Traits>, GridViewType>::value,
                "GridViewType has to be derived from GridView!");
}; // class GenericBaseTraits


} // namespace internal
} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_BASE_INTERNAL_HH
