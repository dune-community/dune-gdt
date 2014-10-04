// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_BASE_INTERNAL_HH
#define DUNE_GDT_PRODUCTS_BASE_INTERNAL_HH

#include <type_traits>

#include <dune/grid/common/gridview.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/walker.hh>

#include <dune/gdt/assembler/local/codim0.hh>
#include <dune/gdt/assembler/local/codim1.hh>
#include <dune/gdt/localoperator/interface.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {
namespace Products {


// forward
template <class LocalOperatorProvider, class RangeImp, class SourceImp>
class LocalizableBase;


// forward
template <class LocalOperatorProvider, class MatrixImp, class RangeSpaceImp, class SourceSpaceImp>
class AssemblableBase;


// forward
template <class LocalOperatorProvider>
class GenericBase;


/**
 * \brief Base class for all local operator providers.
 *
 *        You have to derive a custom class LocalOperatorProvider from this class in order to use \sa LocalizableBase!
 *        All classes derived from this class have to provide the following types:
 *        - `GridViewType`
 *        - `FieldType`
 *        Depending on the type of local operator you want to use for the product (aka provide to LocalizableBase) you
 *        have to additionally implement (you can have zero or one local operator per type each):
 *        - if a local operator is derived from LocalOperator::Codim0Interface you have to provide the public type
 *          - `VolumeOperatorType`
 *          and you have to provide the public members
 *          - `static const bool has_volume_operator = true;`
 *          - `const VolumeOperatorType volume_operator_;
 *          If you want to restrict the entities the local operator will be applied on you have to additionally provide
 *          a method:
 *          - `DSG::ApplyOn::WhichEntity< GridViewType >* entities() const`
 *          See \sa Products::internal::WeightedL2Base in weightedl2-internal.hh for an example of a
 *          LocalOperatorProvider based on a local codim 0 operator.
 *        - if a local operator is derived from LocalOperator::Codim1BoundaryInterface you have to provide the public
 *          type
 *          - `BoundaryOperatorType`
 *          and you have to provide the public members
 *          - `static const bool has_boundary_operator = true;`
 *          - `const BoundaryOperatorType boundary_operator_;`
 *          If you want to restrict the intersections your the local operator will be applied on you have to
 *          additionally provide a method
 *          - `DSG::ApplyOn::WhichIntersections< GridViewType >* boundary_intersections() const`
 *          See \sa Products::internal::BoundaryL2Base in boundaryl2-internal.hh for an example of a
 *          LocalOperatorProvider based on a local codim 1 boundary operator.
 */
template <class GridViewType>
class LocalOperatorProviderBase
{
public:
  static const bool has_volume_operator   = false;
  static const bool has_boundary_operator = false;

  DSG::ApplyOn::AllEntities<GridViewType>* entities() const
  {
    return new DSG::ApplyOn::AllEntities<GridViewType>();
  }

  DSG::ApplyOn::BoundaryIntersections<GridViewType>* boundary_intersections() const
  {
    return new DSG::ApplyOn::BoundaryIntersections<GridViewType>();
  }
}; // class LocalOperatorProviderBase


namespace internal {


/**
 * \brief Helper class for \sa Products::LocalizableBase
 *
 *        This class manages the creation of the appropriate functors needed to handle the local operators provided
 *        by any class derived from \sa LocalOperatorProviderBase.
 * \note  This class is usually not of interest to the average user.
 */
template <class LocalOperatorProvider, class RangeType, class SourceType>
class LocalizableBaseHelper
{
  static_assert(std::is_base_of<LocalOperatorProviderBase<typename LocalOperatorProvider::GridViewType>,
                                LocalOperatorProvider>::value,
                "LocalOperatorProvider has to be derived from LocalOperatorProviderBase!");

  typedef typename LocalOperatorProvider::GridViewType GridViewType;
  typedef typename LocalOperatorProvider::FieldType FieldType;
  typedef Stuff::Grid::Walker<GridViewType> WalkerType;

  template <class LO, bool anthing = false>
  struct Volume
  {
    Volume(const GridViewType&, const LocalOperatorProvider&, const RangeType&, const SourceType&)
    {
    }

    void add(WalkerType&)
    {
    }

    FieldType result() const
    {
      return 0.0;
    }
  }; // struct Volume< ..., false >

  template <class LO>
  struct Volume<LO, true>
  {
    // if you get an error here you have defined has_volume_operator to true but either do not provide
    // VolumeOperatorType or your VolumeOperatorType is not derived from LocalOperator::Codim0Interface
    typedef typename LocalOperatorProvider::VolumeOperatorType LocalOperatorType;
    typedef LocalAssembler::Codim0OperatorAccumulateFunctor<GridViewType, LocalOperatorType, RangeType, SourceType,
                                                            FieldType> FunctorType;

    Volume(const GridViewType& grid_view, const LocalOperatorProvider& local_operators, const RangeType& range,
           const SourceType& source)
      : local_operators_(local_operators)
      , functor_(grid_view, local_operators_.volume_operator_, range, source) // <- if you get an error here you have
    {
    } //    defined has_volume_operator to true
    //    but do not provide volume_operator_

    void add(WalkerType& grid_walker)
    {
      grid_walker.add(functor_, local_operators_.entities()); // <- if you get an error here you have defined
    } //    has_volume_operator to true but implemented the
    //    wrong entities()

    FieldType result() const
    {
      return functor_.result();
    }

    const LocalOperatorProvider& local_operators_;
    FunctorType functor_;
  }; // struct Volume< ..., true >

  template <class LO, bool anthing = false>
  struct Boundary
  {
    Boundary(const GridViewType&, const LocalOperatorProvider&, const RangeType&, const SourceType&)
    {
    }

    void add(WalkerType&)
    {
    }

    FieldType result() const
    {
      return 0.0;
    }
  }; // struct Boundary< ..., false >

  template <class LO>
  struct Boundary<LO, true>
  {
    // if you get an error here you have defined has_boundary_operator to true but either do not provide
    // BoundaryOperatorType or your BoundaryOperatorType is not derived from LocalOperator::Codim1BoundaryInterface
    typedef typename LocalOperatorProvider::BoundaryOperatorType LocalOperatorType;
    typedef LocalAssembler::Codim1BoundaryOperatorAccumulateFunctor<GridViewType, LocalOperatorType, RangeType,
                                                                    SourceType, FieldType> FunctorType;

    Boundary(const GridViewType& grid_view, const LocalOperatorProvider& local_operators, const RangeType& range,
             const SourceType& source)
      : local_operators_(local_operators)
      , functor_(grid_view, local_operators_.boundary_operator_, range, source) // <- if you get an error here you have
    {
    } //    defined has_boundary_operator to
    //    true but do not provide
    //    boundary_operator_

    void add(WalkerType& grid_walker)
    {
      grid_walker.add(functor_, local_operators_.boundary_intersections()); // <- if you get an error here you have
    } //    defined has_boundary_operator to true
    //    but implemented the wrong
    //    boundary_intersections()

    FieldType result() const
    {
      return functor_.result();
    }

    const LocalOperatorProvider& local_operators_;
    FunctorType functor_;
  }; // struct Boundary< ..., true >

public:
  LocalizableBaseHelper(WalkerType& walker, const LocalOperatorProvider& local_operators, const RangeType& range,
                        const SourceType& source)
    : volume_helper_(walker.grid_view(), local_operators, range, source)
    , boundary_helper_(walker.grid_view(), local_operators, range, source)
  {
    volume_helper_.add(walker);
    boundary_helper_.add(walker);
  }

  FieldType result() const
  {
    return volume_helper_.result() + boundary_helper_.result();
  }

private:
  Volume<LocalOperatorProvider, LocalOperatorProvider::has_volume_operator> volume_helper_;
  Boundary<LocalOperatorProvider, LocalOperatorProvider::has_boundary_operator> boundary_helper_;
}; // class LocalizableBaseHelper


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
