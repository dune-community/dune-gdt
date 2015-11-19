// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_VIRTUAL_REFINED_HH
#define DUNE_GDT_OPERATORS_VIRTUAL_REFINED_HH

#include <type_traits>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container/interfaces.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localoperator/virtual_codim0.hh>
#include <dune/gdt/assembler/local/codim0.hh>
#include <dune/gdt/assembler/local/virtual_codim0.hh>
#include <dune/gdt/assembler/system.hh>

#include "base.hh"

namespace Dune
{
namespace GDT
{
namespace Operators
{


// forward
template <class RealEllipticOperatorImp>
class VirtualRefinedEllipticCG;


namespace internal
{


template <class RealEllipticOperatorImp>
class VirtualRefinedEllipticCGTraits
{
public:
  typedef VirtualRefinedEllipticCG<RealEllipticOperatorImp> derived_type;
  typedef typename RealEllipticOperatorImp::MatrixType MatrixType;
  typedef typename RealEllipticOperatorImp::SourceSpaceType SourceSpaceType;
  typedef typename RealEllipticOperatorImp::RangeSpaceType RangeSpaceType;
  typedef typename RealEllipticOperatorImp::GridViewType GridViewType;
  typedef typename RealEllipticOperatorImp::Traits::DiffusionFactorType DiffusionFactorType;
  typedef SystemAssembler<RangeSpaceType, GridViewType, SourceSpaceType> SystemAssemblerType;
}; // class VirtualRefinedEllipticCGTraits


} // namespace internal


template <class RealEllipticOperatorImp>
class VirtualRefinedEllipticCG
    : Stuff::Common::StorageProvider<
          typename internal::VirtualRefinedEllipticCGTraits<RealEllipticOperatorImp>::MatrixType>,
      public Operators::MatrixBased<internal::VirtualRefinedEllipticCGTraits<RealEllipticOperatorImp>>,
      public internal::VirtualRefinedEllipticCGTraits<RealEllipticOperatorImp>::SystemAssemblerType
{

public:
  typedef internal::VirtualRefinedEllipticCGTraits<RealEllipticOperatorImp> Traits;
  using DiffusionType = typename Traits::DiffusionFactorType;

private:
  typedef
      typename internal::VirtualRefinedEllipticCGTraits<RealEllipticOperatorImp>::SystemAssemblerType AssemblerBaseType;
  typedef Operators::MatrixBased<internal::VirtualRefinedEllipticCGTraits<RealEllipticOperatorImp>> OperatorBaseType;
  
//  typedef LocalOperator::Codim0Integral<LocalEvaluation::Elliptic<DiffusionType>> RealLocalOperatorType;
//    typedef LocalAssembler::Codim0Matrix<RealLocalOperatorType> RealLocalAssemblerType;
  
  typedef LocalOperator::VirtualRefinedCodim0Integral<LocalEvaluation::Elliptic<DiffusionType>> VirtualRefinedLocalOperatorType;
  typedef LocalAssembler::VirtualRefinedCodim0Matrix<VirtualRefinedLocalOperatorType> VirtualRefinedLocalAssemblerType;
  
  typedef Stuff::Common::StorageProvider<
      typename internal::VirtualRefinedEllipticCGTraits<RealEllipticOperatorImp>::MatrixType> StorageProvider;

public:
  typedef typename Traits::MatrixType MatrixType;
  typedef typename Traits::SourceSpaceType SourceSpaceType;
  typedef typename Traits::RangeSpaceType RangeSpaceType;
  typedef typename Traits::GridViewType GridViewType;


  VirtualRefinedEllipticCG(const DiffusionType& diffusion, MatrixType& mtrx, const SourceSpaceType& src_spc)
    : StorageProvider(mtrx)
    , OperatorBaseType(this->storage_access(), src_spc)
    , AssemblerBaseType(src_spc)
    , diffusion_(diffusion)
    , local_operator_(diffusion_)
//    , local_assembler_(local_operator_)
    , real_local_assembler_(local_operator_)
    , assembled_(false)
  {
    this->add(real_local_assembler_, this->matrix());
  }
  virtual void assemble() override final
  {
    if (!assembled_) {
      AssemblerBaseType::assemble(true);
      assembled_ = true;
    }
  } // ... assemble(...)
private:
  const DiffusionType& diffusion_;
  const VirtualRefinedLocalOperatorType local_operator_;
//  const RealLocalAssemblerType local_assembler_;
  VirtualRefinedLocalAssemblerType real_local_assembler_;
  bool assembled_;
}; // class VirtualRefinedEllipticCG

} // namespace Operators
} // namespace GDT
} // namespace Dune


#endif // VIRTUAL_REFINED_HH
