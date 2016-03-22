// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_BASE_HH
#define DUNE_GDT_PRODUCTS_BASE_HH

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/tmp-storage.hh>
#include <dune/stuff/grid/walker.hh>
#include <dune/stuff/la/container/pattern.hh>

#include <dune/gdt/assembler/system.hh>

#include "interfaces.hh"
#include "base-internal.hh"

namespace Dune {
namespace GDT {
namespace Products {


/**
 * \brief Base class for all localizable products.
 *
 *        The purpose of this class is to facilitate the implementation of localizable products (that are based on
 *        local operators) by implementing as much as possible. All you have to do is to implement a class
 *        LocalOperatorProvider, derived from LocalOperatorProviderBase. See the documentation of
 *        \sa LocalOperatorProviderBase for the requirements of your class, depending on the type of local operator you
 *        want to use. We do not use the CRTP machinery here but trust you to implement what is necessarry. Any
 *        additional ctor arguments given to LocalizableBase are forwarded to your implementation of
 *        LocalOperatorProvider. Static checks of GridViewType, RangeType and SourceType are performed in
 *        internal::LocalizableBaseTraits. You can finally implement the product by deriving from this class and
 *        providing the appropriate LocalOperatorProvider. To get an idea see \sa Products::internal::WeightedL2Base in
 *        weightedl2-internal.hh for an example of a LocalOperatorProvider and Products::WeightedL2Localizable in
 *        weightedl2.hh for an example of the final product.
 * \note  If you want to know about the internal workings of this class, take a look at \sa
 *        internal::LocalizableBaseHelper.
 */
template< class LocalOperatorProvider, class RangeImp, class SourceImp >
class LocalizableBase
  : public LocalizableProductInterface< internal::LocalizableBaseTraits< LocalOperatorProvider, RangeImp, SourceImp > >
  , public Stuff::Grid::Walker< typename LocalOperatorProvider::GridViewType >
{
  typedef LocalizableBase< LocalOperatorProvider, RangeImp, SourceImp >                 ThisType;
  typedef LocalizableProductInterface
      < internal::LocalizableBaseTraits< LocalOperatorProvider, RangeImp, SourceImp > > ProductBaseType;
  typedef Stuff::Grid::Walker< typename LocalOperatorProvider::GridViewType >           WalkerBaseType;
public:
  typedef internal::LocalizableBaseTraits< LocalOperatorProvider, RangeImp, SourceImp > Traits;

  typedef typename WalkerBaseType::GridViewType     GridViewType;
  typedef typename WalkerBaseType::EntityType       EntityType;
  typedef typename WalkerBaseType::IntersectionType IntersectionType;
  typedef typename ProductBaseType::RangeType  RangeType;
  typedef typename ProductBaseType::SourceType SourceType;
  typedef typename ProductBaseType::FieldType  FieldType;
private:
  typedef internal::LocalizableBaseHelper< LocalOperatorProvider, RangeType, SourceType > HelperType;

public:
  template< class... Args >
  LocalizableBase(const GridViewType& grd_vw, const RangeType& rng, const SourceType& src, Args&& ...args)
    : WalkerBaseType(grd_vw)
    , range_(rng)
    , source_(src)
    , local_operators_(std::forward< Args >(args)...)
    , helper_(*this, local_operators_, range_, source_)
    , walked_(false)
  {}

  template< class... Args >
  LocalizableBase(const GridViewType& grd_vw, const RangeType& rng, Args&& ...args)
    : WalkerBaseType(grd_vw)
    , range_(rng)
    , source_(rng)
    , local_operators_(std::forward< Args >(args)...)
    , helper_(*this, local_operators_, range_, source_)
    , walked_(false)
  {}

  LocalizableBase(const ThisType& other)
    : WalkerBaseType(other.grid_view())
    , range_(other.range_)
    , source_(other.source_)
    , local_operators_(other.local_operators_)
    , helper_(*this, local_operators_, range_, source_)
    , walked_(false)
  {}

  using WalkerBaseType::grid_view;

  const RangeType& range() const
  {
    return range_;
  }

  const SourceType& source() const
  {
    return source_;
  }

  /**
   *       This method can be used to compute a local semi product on one entity, e.g. in the context of error
   *       estimation.
   * \note Calling this method should not alter the result obtained by apply2.
   */
  FieldType compute_locally(const EntityType& entity)
  {
    return helper_.compute_locally(entity);
  }

  /**
   *       This method can be used to compute a local semi product on one intersection, e.g. in the context of error
   *       estimation.
   * \note Calling this method should not alter the result obtained by apply2.
   */
  FieldType compute_locally(const IntersectionType& intersection,
                            const EntityType& inside_entity,
                            const EntityType& outside_entity)
  {
    return helper_.compute_locally(intersection, inside_entity, outside_entity);
  }

  virtual void finalize() override
  {
    walked_ = true;
    WalkerBaseType::finalize();
  }

  FieldType apply2()
  {
    if (!walked_) {
      this->walk();
      walked_ = true;
    }
    return helper_.result();
  } // ... apply2(...)

private:
  const RangeType& range_;
  const SourceType& source_;
  const LocalOperatorProvider local_operators_;
  HelperType helper_;
  bool walked_;
}; // class LocalizableBase


/**
 * \brief Base class for all assembable products.
 *
 *        The purpose of this class is to facilitate the implementation of assembable products (that are based on
 *        local operators) by implementing as much as possible. All you have to do is to implement a class
 *        LocalOperatorProvider, derived from LocalOperatorProviderBase. See the documentation of
 *        \sa LocalOperatorProviderBase for the requirements of your class, depending on the type of local operator you
 *        want to use. We do not use the CRTP machinery here but trust you to implement what is necessarry. Any
 *        additional ctor arguments given to AssemblableBase are forwarded to your implementation of
 *        LocalOperatorProvider. Static checks of MatrixType, RangeSpaceType and SourceSpaceType are performed in
 *        internal::AssemblableBaseTraits. You can finally implement the product by deriving from this class and
 *        providing the appropriate LocalOperatorProvider. To get an idea see \sa Products::internal::WeightedL2Base in
 *        weightedl2-internal.hh for an example of a LocalOperatorProvider and Products::WeightedL2Assemblable in
 *        weightedl2.hh for an example of the final product.
 * \note  If you want to know about the internal workings of this class, take a look at \sa
 *        internal::AssemblableBaseHelper.
 * \note  This class proves an automatic creation of the matrix that is being assembled into. Thus each kind of ctor
 *        exists in two variants (one taking a matrix reference, one not). If you provide an external matrix you are
 *        responsible for the well being of the matrix:
 *        - If this is a sparse matrix it needs to have the correct pattern (can be obtained from this class).
 *        - During assemble values are added to the matrix (not set). Thus it should be empty beforehand or you have to
 *          know what you are doing!
 */
template< class LocalOperatorProvider, class MatrixImp, class RangeSpaceImp, class SourceSpaceImp >
class AssemblableBase
  : DSC::StorageProvider< MatrixImp >
  , public AssemblableProductInterface< internal::AssemblableBaseTraits< LocalOperatorProvider,
                                                                         MatrixImp,
                                                                         RangeSpaceImp,
                                                                         SourceSpaceImp > >
  , public SystemAssembler< RangeSpaceImp, typename LocalOperatorProvider::GridViewType, SourceSpaceImp >
{
  typedef DSC::StorageProvider< MatrixImp >                                                 MatrixProvider;
  typedef AssemblableProductInterface< internal::AssemblableBaseTraits< LocalOperatorProvider,
                                                                        MatrixImp,
                                                                        RangeSpaceImp,
                                                                        SourceSpaceImp > >  ProductBaseType;
  typedef SystemAssembler
      < RangeSpaceImp, typename LocalOperatorProvider::GridViewType, SourceSpaceImp >       AssemblerBaseType;
  typedef AssemblableBase< LocalOperatorProvider, MatrixImp, RangeSpaceImp, SourceSpaceImp > ThisType;
public:
  typedef internal::AssemblableBaseTraits
      < LocalOperatorProvider, MatrixImp, RangeSpaceImp, SourceSpaceImp > Traits;
  typedef typename ProductBaseType::GridViewType    GridViewType;
  typedef typename ProductBaseType::RangeSpaceType  RangeSpaceType;
  typedef typename ProductBaseType::SourceSpaceType SourceSpaceType;
  typedef typename ProductBaseType::MatrixType      MatrixType;
  typedef typename ProductBaseType::FieldType       FieldType;
private:
  typedef internal::AssemblableBaseHelper< ThisType, LocalOperatorProvider > HelperType;

public:
  using ProductBaseType::pattern;

  static Stuff::LA::SparsityPatternDefault pattern(const RangeSpaceType& range_space,
                                                   const SourceSpaceType& source_space,
                                                   const GridViewType& grid_view)
  {
    if (LocalOperatorProvider::has_coupling_operator)
      return range_space.compute_face_and_volume_pattern(grid_view, source_space);
    else if (LocalOperatorProvider::has_volume_operator || LocalOperatorProvider::has_boundary_operator)
      return range_space.compute_volume_pattern(grid_view, source_space);
    else
      return Stuff::LA::SparsityPatternDefault();
  }

  template< class... Args >
  AssemblableBase(MatrixType& mtrx,
                  const RangeSpaceType& rng_spc,
                  const GridViewType& grd_vw,
                  const SourceSpaceType& src_spc,
                  Args&& ...args)
    : MatrixProvider(mtrx)
    , AssemblerBaseType(rng_spc, src_spc, grd_vw)
    , local_operators_(std::forward< Args >(args)...)
    , helper_(*this, MatrixProvider::storage_access(), local_operators_)
    , assembled_(false)
  {}

  template< class... Args >
  AssemblableBase(const RangeSpaceType& rng_spc,
                  const GridViewType& grd_vw,
                  const SourceSpaceType& src_spc,
                  Args&& ...args)
    : MatrixProvider(new MatrixType(rng_spc.mapper().size(),
                                     src_spc.mapper().size(),
                                     pattern(rng_spc, src_spc, grd_vw)))
    , AssemblerBaseType(rng_spc, src_spc, grd_vw)
    , local_operators_(std::forward< Args >(args)...)
    , helper_(*this, MatrixProvider::storage_access(), local_operators_)
    , assembled_(false)
  {}

  template< class... Args >
  AssemblableBase(MatrixType& mtrx, const RangeSpaceType& rng_spc, const GridViewType& grd_vw, Args&& ...args)
    : MatrixProvider(mtrx)
    , AssemblerBaseType(rng_spc, grd_vw)
    , local_operators_(std::forward< Args >(args)...)
    , helper_(*this, MatrixProvider::storage_access(), local_operators_)
    , assembled_(false)
  {}

  template< class... Args >
  AssemblableBase(const RangeSpaceType& rng_spc, const GridViewType& grd_vw, Args&& ...args)
    : MatrixProvider(new MatrixType(rng_spc.mapper().size(), rng_spc.mapper().size(), pattern(rng_spc, grd_vw)))
    , AssemblerBaseType(rng_spc, grd_vw)
    , local_operators_(std::forward< Args >(args)...)
    , helper_(*this, MatrixProvider::storage_access(), local_operators_)
    , assembled_(false)
  {}

  template< class... Args >
  AssemblableBase(MatrixType& mtrx, const RangeSpaceType& rng_spc, Args&& ...args)
    : MatrixProvider(mtrx)
    , AssemblerBaseType(rng_spc)
    , local_operators_(std::forward< Args >(args)...)
    , helper_(*this, MatrixProvider::storage_access(), local_operators_)
    , assembled_(false)
  {}

  template< class... Args >
  AssemblableBase(const RangeSpaceType& rng_spc, Args&& ...args)
    : MatrixProvider(new MatrixType(rng_spc.mapper().size(), rng_spc.mapper().size(), pattern(rng_spc)))
    , AssemblerBaseType(rng_spc)
    , local_operators_(std::forward< Args >(args)...)
    , helper_(*this, MatrixProvider::storage_access(), local_operators_)
    , assembled_(false)
  {}

  using AssemblerBaseType::grid_view;

  const RangeSpaceType& range_space() const
  {
    return AssemblerBaseType::test_space();
  }

  const SourceSpaceType& source_space() const
  {
    return AssemblerBaseType::ansatz_space();
  }

  MatrixType& matrix()
  {
    return MatrixProvider::storage_access();
  }

  const MatrixType& matrix() const
  {
    return MatrixProvider::storage_access();
  }

  void assemble(const bool use_tbb = false)
  {
    if (!assembled_) {
      AssemblerBaseType::assemble(use_tbb);
      assembled_ = true;
    }
  } // ... assemble()

private:
  const LocalOperatorProvider local_operators_;
  HelperType helper_;
  bool assembled_;
}; // class AssemblableBase


/**
 * \brief Base class for all generic products.
 *
 *        The purpose of this class is to facilitate the implementation of products that are based on a local operator
 *        that is derived from LocalOperator::Codim0Interface by implementing as much as possible. All you have to do is
 *        to implement a class LocalOperatorProvider that provides the following in addition to the requirements of \sa
 *        LocalizableBase:
 *        - it has to be copyable
 *        GenericBase derives from that provided class and forwards any additional ctor arguments to its ctor. A static
 *        check of GridViewType is performed in internal::GenericBaseTraits. You can finally implement the product by
 *        deriving from this class and providing the appropriate LocalOperatorProvider. To get an idea see \sa
 *        Products::internal::WeightedL2Base in l2-internal.hh for an example of a LocalOperatorProvider and
 *        Products::WeightedL2 in l2.hh for an example of the final product. This product is implemented by a creating a
 *        LocalizableBase product of appropriate type for the given source and range.
 */
template< class LocalOperatorProvider >
class GenericBase
  : LocalOperatorProvider
  , public ProductInterface< internal::GenericBaseTraits< LocalOperatorProvider > >
{
  typedef ProductInterface< internal::GenericBaseTraits< LocalOperatorProvider > > BaseType;
public:
  typedef internal::GenericBaseTraits< LocalOperatorProvider > Traits;
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::FieldType    FieldType;

  typedef typename GridViewType::template Codim< 0 >::Entity  EntityType;
  typedef typename GridViewType::ctype                        DomainFieldType;
  static const size_t                                         dimDomain = GridViewType::dimension;

  template< class... Args >
  GenericBase(const GridViewType& grd_vw, Args&& ...args)
    : LocalOperatorProvider(std::forward< Args >(args)...)
    , grid_view_(grd_vw)
  {}

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  template< class RR, size_t rRR, size_t rCR, class RS, size_t rRS, size_t rCS >
  FieldType apply2(const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, RR, rRR, rCR >& range,
                   const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, RS, rRS, rCS >& source) const
  {
    typedef Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, RR, rRR, rCR > RangeType;
    typedef Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, RS, rRS, rCS > SourceType;
    LocalizableBase< LocalOperatorProvider, RangeType, SourceType > product(grid_view_,
                                                                            range,
                                                                            source,
                                                                            static_cast< const LocalOperatorProvider& >(*this));
    return product.apply2();
  } // ... apply2(...)

private:
  const GridViewType& grid_view_;
}; // class GenericBase


} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_BASE_HH
