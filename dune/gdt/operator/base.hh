// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATOR_BASE_HH
#define DUNE_GDT_OPERATOR_BASE_HH

#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/assembler/local/codim0.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace Operator {


template< class Traits >
class AssemblableVolumeBase
    : public AssemblableOperatorInterface< Traits >
    , public Functor::Codim0< typename Traits::GridViewType >
{
  typedef AssemblableOperatorInterface< Traits > InterfaceType;
public:
  typedef typename Traits::GridViewType     GridViewType;
  typedef typename Traits::SourceSpaceType  SourceSpaceType;
  typedef typename Traits::RangeSpaceType   RangeSpaceType;
  typedef typename Traits::MatrixType       MatrixType;

  typedef typename MatrixType::ScalarType   FieldType;
private:
  typedef TmpStorageProvider::Matrices< FieldType >         TmpMatricesProviderType;
  typedef typename Traits::LocalOperatorType                LocalOperatorType;
  typedef LocalAssembler::Codim0Matrix< LocalOperatorType > LocalAssemblerType;
  typedef Stuff::LA::Solver< MatrixType >                   LinearSolverType;

public:
  using typename InterfaceType::EntityType;

  AssemblableVolumeBase(MatrixType& matrix,
                        const SourceSpaceType& source_space,
                        const RangeSpaceType& range_space,
                        const GridViewType& grid_view)
    : matrix_(matrix)
    , source_space_(source_space)
    , range_space_(range_space)
    , grid_view_(grid_view)
    , local_assembler_(nullptr)
    , tmp_storage_(nullptr)
    , linear_solver_(nullptr)
    , prepared_(false)
    , assembled_(false)
  {}

  AssemblableVolumeBase(MatrixType& matrix,
                        const SourceSpaceType& source_space,
                        const RangeSpaceType& range_space)
    : matrix_(matrix)
    , source_space_(source_space)
    , range_space_(range_space)
    , grid_view_(*(source_space_.grid_view()))
    , local_assembler_(nullptr)
    , tmp_storage_(nullptr)
    , linear_solver_(nullptr)
    , prepared_(false)
    , assembled_(false)
  {}

  AssemblableVolumeBase(MatrixType& matrix,
                        const SourceSpaceType& source_space)
    : matrix_(matrix)
    , source_space_(source_space)
    , range_space_(source_space_)
    , grid_view_(*(source_space_.grid_view()))
    , local_assembler_(nullptr)
    , tmp_storage_(nullptr)
    , linear_solver_(nullptr)
    , prepared_(false)
    , assembled_(false)
  {}

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
  virtual void prepare()
  {
    if (!assembled_ && !prepared_) {
      local_assembler_ = std::unique_ptr< LocalAssemblerType >(new LocalAssemblerType(local_operator()));
      tmp_storage_
          = std::unique_ptr< TmpMatricesProviderType >(new TmpMatricesProviderType(local_assembler_->numTmpObjectsRequired(),
                                                                                   range_space_.mapper().maxNumDofs(),
                                                                                   source_space_.mapper().maxNumDofs()));
      prepared_ = true;
    }
  } // ... prepare()

  virtual void apply_local(const EntityType& entity)
  {
    assert(prepared_);
    assert(local_assembler_);
    assert(tmp_storage_);
    local_assembler_->assembleLocal(range_space_, source_space_,
                                    entity,
                                    matrix(),
                                    tmp_storage_->matrices(), tmp_storage_->indices());
  } // ... apply_local(...)

  virtual void finalize()
  {
    if (!linear_solver_)
      linear_solver_ = std::unique_ptr< LinearSolverType >(new LinearSolverType(matrix_));
  }

  void assemble()
  {
    if (!assembled_) {
      GridWalker< GridViewType > grid_walker(grid_view_);
      grid_walker.add(*this);
      grid_walker.walk();
      assembled_ = true;
    }
  } // ... assemble()

  template< class S, class R >
  void apply(const Stuff::LA::VectorInterface< S >& source, Stuff::LA::VectorInterface< R >& range)
  {
    typedef typename S::derived_type SourceType;
    typedef typename R::derived_type RangeType;
    assemble();
    matrix_.mv(static_cast< const SourceType& >(source), static_cast< RangeType& >(range));
  }

  static std::vector< std::string > invert_options()
  {
    return LinearSolverType::options();
  }

  static Stuff::Common::ConfigTree invert_options(const std::string& type)
  {
    return LinearSolverType::options(type);
  }

  template< class R, class S >
  void apply_inverse(const Stuff::LA::VectorInterface< R >& range,
                     Stuff::LA::VectorInterface< S >& source,
                     const Stuff::Common::ConfigTree& opts)
  {
    typedef typename S::derived_type SourceType;
    typedef typename R::derived_type RangeType;
    assemble();
    assert(linear_solver_);
    linear_solver_->apply(static_cast< const RangeType& >(range), static_cast< SourceType& >(source), opts);
  }

private:
  MatrixType& matrix_;
  const SourceSpaceType& source_space_;
  const RangeSpaceType& range_space_;
  const GridViewType& grid_view_;
  std::unique_ptr< LocalAssemblerType > local_assembler_;
  std::unique_ptr< TmpMatricesProviderType > tmp_storage_;
  std::unique_ptr< LinearSolverType > linear_solver_;
  bool prepared_;
  bool assembled_;
}; // class AssemblableVolumeBase


} // namespace Operator
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATOR_BASE_HH
