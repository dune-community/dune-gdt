// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_BASE_HH
#define DUNE_GDT_PRODUCTS_BASE_HH

#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/assembler/tmp-storage.hh>
#include <dune/gdt/assembler/local/codim0.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace Products {


template <class Traits>
class LocalizableBase : public LocalizableProductInterface<Traits>,
                        public Functor::Codim0<typename Traits::GridViewType>
{
  typedef LocalizableProductInterface<Traits> InterfaceType;

public:
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::SourceType SourceType;
  typedef typename Traits::FieldType FieldType;

private:
  typedef typename Traits::LocalOperatorType LocalOperatorType;
  typedef TmpStorageProvider::Matrices<FieldType> TmpMatricesProviderType;

public:
  using typename InterfaceType::EntityType;

public:
  LocalizableBase(const GridViewType& grid_view, const RangeType& range, const SourceType& source)
    : grid_view_(grid_view)
    , range_(range)
    , source_(source)
    , tmp_storage_(nullptr)
    , prepared_(false)
    , finalized_(false)
    , result_(0)
  {
  }

  LocalizableBase(const GridViewType& grid_view, const RangeType& range)
    : grid_view_(grid_view)
    , range_(range)
    , source_(range_)
    , tmp_storage_(nullptr)
    , prepared_(false)
    , finalized_(false)
    , result_(0)
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

private:
  virtual const LocalOperatorType& local_operator() const = 0;

public:
  virtual void prepare() DS_OVERRIDE
  {
    if (!prepared_) {
      tmp_storage_ = std::unique_ptr<TmpMatricesProviderType>(
          new TmpMatricesProviderType({1, local_operator().numTmpObjectsRequired()}, 1, 1));
      result_ *= 0.0;
      prepared_ = true;
    }
  } // ... prepare()

  virtual void apply_local(const EntityType& entity) DS_OVERRIDE
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
    local_operator().apply(*local_range_ptr, *local_source_ptr, local_operator_result, tmp_matrices);
    assert(local_operator_result.rows() == 1);
    assert(local_operator_result.cols() == 1);
    result_ += local_operator_result[0][0];
  } // ... apply_local(...)

  virtual void finalize() DS_OVERRIDE
  {
    if (!finalized_) {
      result_    = grid_view_.comm().sum(result_);
      finalized_ = true;
    }
  }

  FieldType apply2()
  {
    if (!finalized_) {
      GridWalker<GridViewType> grid_walker(grid_view_);
      grid_walker.add(*this);
      grid_walker.walk();
    }
    return result_;
  }

private:
  const GridViewType& grid_view_;
  const RangeType& range_;
  const SourceType& source_;
  std::unique_ptr<TmpMatricesProviderType> tmp_storage_;
  bool prepared_;
  bool finalized_;
  FieldType result_;
}; // class LocalizableBase


template <class Traits>
class AssemblableBase : public AssemblableProductInterface<Traits>,
                        public Functor::Codim0<typename Traits::GridViewType>
{
  typedef AssemblableProductInterface<Traits> InterfaceType;

public:
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::RangeSpaceType RangeSpaceType;
  typedef typename Traits::SourceSpaceType SourceSpaceType;
  typedef typename Traits::MatrixType MatrixType;

  typedef typename MatrixType::ScalarType FieldType;

private:
  typedef TmpStorageProvider::Matrices<FieldType> TmpMatricesProviderType;
  typedef typename Traits::LocalOperatorType LocalOperatorType;
  typedef LocalAssembler::Codim0Matrix<LocalOperatorType> LocalAssemblerType;

public:
  using typename InterfaceType::EntityType;

  using InterfaceType::pattern;

  AssemblableBase(MatrixType& matrix, const RangeSpaceType& range_space, const GridViewType& grid_view,
                  const SourceSpaceType& source_space)
    : matrix_(matrix)
    , range_space_(range_space)
    , grid_view_(grid_view)
    , source_space_(source_space)
    , local_assembler_(nullptr)
    , tmp_storage_(nullptr)
    , prepared_(false)
    , assembled_(false)
  {
  }

  AssemblableBase(MatrixType& matrix, const RangeSpaceType& range_space, const GridViewType& grid_view)
    : matrix_(matrix)
    , range_space_(range_space)
    , grid_view_(grid_view)
    , source_space_(range_space_)
    , local_assembler_(nullptr)
    , tmp_storage_(nullptr)
    , prepared_(false)
    , assembled_(false)
  {
  }

  AssemblableBase(MatrixType& matrix, const RangeSpaceType& range_space)
    : matrix_(matrix)
    , range_space_(range_space)
    , grid_view_(*(range_space_->grid_view()))
    , source_space_(range_space_)
    , local_assembler_(nullptr)
    , tmp_storage_(nullptr)
    , prepared_(false)
    , assembled_(false)
  {
  }

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  const RangeSpaceType& range_space() const
  {
    return range_space_;
  }

  const SourceSpaceType& source_space() const
  {
    return source_space_;
  }

  MatrixType& matrix()
  {
    return matrix_;
  }

  const MatrixType& matrix() const
  {
    return matrix_;
  }

private:
  virtual const LocalOperatorType& local_operator() const = 0;

public:
  virtual void prepare() DS_OVERRIDE
  {
    if (!assembled_ && !prepared_) {
      local_assembler_ = std::unique_ptr<LocalAssemblerType>(new LocalAssemblerType(local_operator()));
      tmp_storage_ = std::unique_ptr<TmpMatricesProviderType>(
          new TmpMatricesProviderType(local_assembler_->numTmpObjectsRequired(),
                                      range_space_.mapper().maxNumDofs(),
                                      source_space_.mapper().maxNumDofs()));
      prepared_ = true;
    }
  } // ... prepare()

  virtual void apply_local(const EntityType& entity) DS_OVERRIDE
  {
    assert(prepared_);
    assert(local_assembler_);
    assert(tmp_storage_);
    local_assembler_->assembleLocal(
        range_space_, source_space_, entity, matrix_, tmp_storage_->matrices(), tmp_storage_->indices());
  } // ... apply_local(...)

  void assemble()
  {
    if (!assembled_) {
      GridWalker<GridViewType> grid_walker(grid_view_);
      grid_walker.add(*this);
      grid_walker.walk();
      assembled_ = true;
    }
  } // ... assemble()

  using InterfaceType::apply2;

private:
  MatrixType& matrix_;
  const RangeSpaceType& range_space_;
  const GridViewType& grid_view_;
  const SourceSpaceType& source_space_;
  std::unique_ptr<LocalAssemblerType> local_assembler_;
  std::unique_ptr<TmpMatricesProviderType> tmp_storage_;
  bool prepared_;
  bool assembled_;
}; // class AssemblableBase


} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_BASE_HH
