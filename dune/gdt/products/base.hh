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

#include <dune/gdt/assembler/local/codim0.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/localoperator/interface.hh>

#include "interfaces.hh"
#include "base-internal.hh"

namespace Dune {
namespace GDT {
namespace Products {


/**
 * \brief Base class for all localizable products.
 *
 *        The purpose of this class is to facilitate the implementation of localizable products that are based on a
 *        local operator that is derived from LocalOperator::Codim0Interface by implementing as much as possible. All
 *        you have to do is to implement a class LocalOperatorProvider that provides the following:
 *        - a protected member local_operator_
 *        - a typedef GridViewType
 *        - a typedef FieldType
 *        LocalizableBase derives from that provided class and forwards any additional ctor arguments to its ctor.
 *        Static checks of GridViewType, RangeType and SourceType are performed in internal::LocalizableBaseTraits. You
 *        can finally implement the product by deriving from this class and providing the appropriate
 *        LocalOperatorProvider. To get an idea see \sa Products::internal::WeightedL2Base in l2-internal.hh for an
 *        example of a LocalOperatorProvider and Products::WeightedL2Localizable in l2.hh for an example of the final
 *        product.
 */
template <class LocalOperatorProvider, class RangeImp, class SourceImp>
class LocalizableBase
    : LocalOperatorProvider,
      public LocalizableProductInterface<internal::LocalizableBaseTraits<LocalOperatorProvider, RangeImp, SourceImp>>,
      public Stuff::Grid::Functor::Codim0<typename LocalOperatorProvider::GridViewType>
{
  typedef LocalizableProductInterface<internal::LocalizableBaseTraits<LocalOperatorProvider, RangeImp, SourceImp>>
      BaseType;

public:
  typedef internal::LocalizableBaseTraits<LocalOperatorProvider, RangeImp, SourceImp> Traits;

  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::SourceType SourceType;
  typedef typename BaseType::FieldType FieldType;
  typedef typename BaseType::EntityType EntityType;

private:
  typedef DSC::TmpMatricesStorage<FieldType> TmpMatricesProviderType;

public:
  template <class... Args>
  LocalizableBase(const GridViewType& grd_vw, const RangeType& rng, const SourceType& src, Args&&... args)
    : LocalOperatorProvider(std::forward<Args>(args)...)
    , grid_view_(grd_vw)
    , range_(rng)
    , source_(src)
    , tmp_storage_(nullptr)
    , prepared_(false)
    , finalized_(false)
    , result_(0)
    , finalized_result_(0)
  {
  }

  template <class... Args>
  LocalizableBase(const GridViewType& grd_vw, const RangeType& rng, Args&&... args)
    : LocalOperatorProvider(std::forward<Args>(args)...)
    , grid_view_(grd_vw)
    , range_(rng)
    , source_(range_)
    , tmp_storage_(nullptr)
    , prepared_(false)
    , finalized_(false)
    , result_(0)
    , finalized_result_(0)
  {
  }

  virtual ~LocalizableBase()
  {
  }

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  const RangeType& range() const
  {
    return range_;
  }

  const SourceType& source() const
  {
    return source_;
  }

  virtual void prepare() DS_OVERRIDE
  {
    if (!prepared_) {
      tmp_storage_ = std::unique_ptr<TmpMatricesProviderType>(
          new TmpMatricesProviderType({1, this->local_operator_.numTmpObjectsRequired()}, 1, 1));
      result_   = FieldType(0.0);
      prepared_ = true;
    }
  } // ... prepare()

  FieldType compute_locally(const EntityType& entity) const
  {
    assert(prepared_);
    assert(tmp_storage_);
    auto& tmp_storage = tmp_storage_->matrices();
    assert(tmp_storage.size() >= 2);
    assert(tmp_storage[0].size() >= 1);
    auto& local_operator_result = tmp_storage[0][0];
    auto& tmp_matrices          = tmp_storage[1];
    // get the local functions
    const auto local_source_ptr = this->source().local_function(entity);
    const auto local_range_ptr  = this->range().local_function(entity);
    // apply local operator
    this->local_operator_.apply(*local_range_ptr, *local_source_ptr, local_operator_result, tmp_matrices);
    assert(local_operator_result.rows() == 1);
    assert(local_operator_result.cols() == 1);
    return local_operator_result[0][0];
  } // ... compute_locally(...)

  virtual void apply_local(const EntityType& entity) DS_OVERRIDE
  {
    *result_ += compute_locally(entity);
  }

  virtual void finalize() DS_OVERRIDE
  {
    if (!finalized_) {
      finalized_result_ = result_.sum();
      finalized_result_ = grid_view_.comm().sum(finalized_result_);
      finalized_        = true;
    }
  }

  FieldType apply2()
  {
    if (!finalized_) {
      Stuff::Grid::Walker<GridViewType> grid_walker(grid_view_);
      grid_walker.add(*this);
      grid_walker.walk();
    }
    return finalized_result_;
  }

private:
  const GridViewType& grid_view_;
  const RangeType& range_;
  const SourceType& source_;
  std::unique_ptr<TmpMatricesProviderType> tmp_storage_;
  bool prepared_;
  bool finalized_;
  DS::PerThreadValue<FieldType> result_;
  FieldType finalized_result_;
}; // class LocalizableBase


/**
 * \brief Base class for all assembable products.
 *
 *        The purpose of this class is to facilitate the implementation of assembable products that are based on a
 *        local operator that is derived from LocalOperator::Codim0Interface by implementing as much as possible. All
 *        you have to do is to implement a class LocalOperatorProvider that provides the following:
 *        - a protected member local_operator_
 *        - a typedef GridViewType
 *        AssemblableBase derives from that provided class and forwards any additional ctor arguments to its ctor.
 *        Static checks of MatrixType, RangeSpaceType and SourceSpaceType are performed in
 *        internal::AssemblableBaseTraits. You can finally implement the product by deriving from this class and
 *        providing the appropriate LocalOperatorProvider. To get an idea see \sa Products::internal::WeightedL2Base in
 *        l2-internal.hh for an example of a LocalOperatorProvider and Products::WeightedL2Assemblable in l2.hh for an
 *        example of the final product.
 * \note  This class proves an automatic creation of the matrix that is being assembled into. Thus each kind of ctor
 *        exists in two variants (one taking a matrix reference, one not). If you provide an external matrix you are
 *        responsible for the well being of the matrix:
 *        - If this is a sparse matrix it needs to have the correct pattern (can be obtained from this class).
 *        - During assemble values are added to the matrix (not set). Thus it should be empty beforehand or you have to
 *          know what you are doing!
 */
template <class LocalOperatorProvider, class MatrixImp, class RangeSpaceImp, class SourceSpaceImp>
class AssemblableBase
    : LocalOperatorProvider,
      DSC::StorageProvider<MatrixImp>,
      public AssemblableProductInterface<internal::AssemblableBaseTraits<LocalOperatorProvider, MatrixImp,
                                                                         RangeSpaceImp, SourceSpaceImp>>,
      public SystemAssembler<RangeSpaceImp, typename LocalOperatorProvider::GridViewType, SourceSpaceImp>
{
  typedef DSC::StorageProvider<MatrixImp> StorageBaseType;
  typedef AssemblableProductInterface<internal::AssemblableBaseTraits<LocalOperatorProvider, MatrixImp, RangeSpaceImp,
                                                                      SourceSpaceImp>> ProductBaseType;
  typedef SystemAssembler<RangeSpaceImp, typename LocalOperatorProvider::GridViewType, SourceSpaceImp>
      AssemblerBaseType;
  typedef LocalAssembler::Codim0Matrix<typename LocalOperatorProvider::LocalOperatorType> LocalAssemblerType;

public:
  typedef internal::AssemblableBaseTraits<LocalOperatorProvider, MatrixImp, RangeSpaceImp, SourceSpaceImp> Traits;
  typedef typename ProductBaseType::GridViewType GridViewType;
  typedef typename ProductBaseType::RangeSpaceType RangeSpaceType;
  typedef typename ProductBaseType::SourceSpaceType SourceSpaceType;
  typedef typename ProductBaseType::MatrixType MatrixType;
  typedef typename ProductBaseType::FieldType FieldType;

  using ProductBaseType::pattern;

  static Stuff::LA::SparsityPatternDefault pattern(const RangeSpaceType& range_space,
                                                   const SourceSpaceType& source_space, const GridViewType& grid_view)
  {
    return range_space.compute_volume_pattern(grid_view, source_space);
  }

  template <class... Args>
  AssemblableBase(MatrixType& mtrx, const RangeSpaceType& rng_spc, const GridViewType& grd_vw,
                  const SourceSpaceType& src_spc, Args&&... args)
    : LocalOperatorProvider(std::forward<Args>(args)...)
    , StorageBaseType(mtrx)
    , AssemblerBaseType(rng_spc, src_spc, grd_vw)
    , local_assembler_(this->local_operator_)
    , assembled_(false)
  {
    setup();
  }

  template <class... Args>
  AssemblableBase(const RangeSpaceType& rng_spc, const GridViewType& grd_vw, const SourceSpaceType& src_spc,
                  Args&&... args)
    : LocalOperatorProvider(std::forward<Args>(args)...)
    , StorageBaseType(
          new MatrixType(rng_spc.mapper().size(), src_spc.mapper().size(), pattern(rng_spc, src_spc, grd_vw)))
    , AssemblerBaseType(rng_spc, src_spc, grd_vw)
    , local_assembler_(this->local_operator_)
    , assembled_(false)
  {
    setup();
  }

  template <class... Args>
  AssemblableBase(MatrixType& mtrx, const RangeSpaceType& rng_spc, const GridViewType& grd_vw, Args&&... args)
    : LocalOperatorProvider(std::forward<Args>(args)...)
    , StorageBaseType(mtrx)
    , AssemblerBaseType(rng_spc, grd_vw)
    , local_assembler_(this->local_operator_)
    , assembled_(false)
  {
    setup();
  }

  template <class... Args>
  AssemblableBase(const RangeSpaceType& rng_spc, const GridViewType& grd_vw, Args&&... args)
    : LocalOperatorProvider(std::forward<Args>(args)...)
    , StorageBaseType(new MatrixType(rng_spc.mapper().size(), rng_spc.mapper().size(), pattern(rng_spc, grd_vw)))
    , AssemblerBaseType(rng_spc, grd_vw)
    , local_assembler_(this->local_operator_)
    , assembled_(false)
  {
    setup();
  }

  template <class... Args>
  AssemblableBase(MatrixType& mtrx, const RangeSpaceType& rng_spc, Args&&... args)
    : LocalOperatorProvider(std::forward<Args>(args)...)
    , StorageBaseType(mtrx)
    , AssemblerBaseType(rng_spc)
    , local_assembler_(this->local_operator_)
    , assembled_(false)
  {
    setup();
  }

  template <class... Args>
  AssemblableBase(const RangeSpaceType& rng_spc, Args&&... args)
    : LocalOperatorProvider(std::forward<Args>(args)...)
    , StorageBaseType(new MatrixType(rng_spc.mapper().size(), rng_spc.mapper().size(), pattern(rng_spc)))
    , AssemblerBaseType(rng_spc)
    , local_assembler_(this->local_operator_)
    , assembled_(false)
  {
    setup();
  }

  const GridViewType& grid_view() const
  {
    return AssemblerBaseType::grid_view();
  }

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
    return StorageBaseType::storage_access();
  }

  const MatrixType& matrix() const
  {
    return StorageBaseType::storage_access();
  }

  void assemble()
  {
    if (!assembled_) {
      AssemblerBaseType::assemble();
      assembled_ = true;
    }
  } // ... assemble()

private:
  void setup()
  {
    this->add(local_assembler_, matrix());
  }

  const LocalAssemblerType local_assembler_;
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
template <class LocalOperatorProvider>
class GenericBase : LocalOperatorProvider, public ProductInterface<internal::GenericBaseTraits<LocalOperatorProvider>>
{
  typedef ProductInterface<internal::GenericBaseTraits<LocalOperatorProvider>> BaseType;

public:
  typedef internal::GenericBaseTraits<LocalOperatorProvider> Traits;
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::FieldType FieldType;

  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;

  template <class... Args>
  GenericBase(const GridViewType& grd_vw, Args&&... args)
    : LocalOperatorProvider(std::forward<Args>(args)...)
    , grid_view_(grd_vw)
  {
  }

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  template <class RR, int rRR, int rCR, class RS, int rRS, int rCS>
  FieldType
  apply2(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RR, rRR, rCR>& range,
         const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RS, rRS, rCS>& source) const
  {
    typedef Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RR, rRR, rCR> RangeType;
    typedef Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RS, rRS, rCS> SourceType;
    LocalizableBase<LocalOperatorProvider, RangeType, SourceType> product(
        grid_view_, range, source, static_cast<const LocalOperatorProvider&>(*this));
    return product.apply2();
  } // ... apply2(...)

private:
  const GridViewType& grid_view_;
}; // class GenericBase


} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_BASE_HH
