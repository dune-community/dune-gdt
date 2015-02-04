// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_ELLIPTIC_INTERNAL_HH
#define DUNE_GDT_PRODUCTS_ELLIPTIC_INTERNAL_HH

#include <type_traits>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/interfaces.hh>

#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localevaluation/elliptic.hh>

#include "base-internal.hh"

namespace Dune {
namespace GDT {
namespace Products {
namespace internal {


/**
 * \brief Base class for all elliptic products.
 * \note  Most likely you do not want to use this class directly, but Products::EllipticLocalizable,
 *        Products::EllipticAssemblable or Products::Elliptic instead!
 */
template <class DiffusionFactorType, class GV, class FieldImp, class DiffusionTensorType = void>
class EllipticBase : public LocalOperatorProviderBase<GV>
{
  static_assert(Stuff::is_localizable_function<DiffusionFactorType>::value,
                "DiffusionFactorType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(Stuff::is_localizable_function<DiffusionTensorType>::value,
                "DiffusionTensorType has to be derived from Stuff::LocalizableFunctionInterface!");
  typedef EllipticBase<DiffusionFactorType, GV, FieldImp, DiffusionTensorType> ThisType;

public:
  typedef GV GridViewType;
  typedef FieldImp FieldType;
  typedef LocalOperator::Codim0Integral<LocalEvaluation::Elliptic<DiffusionFactorType, DiffusionTensorType>>
      VolumeOperatorType;

  static const bool has_volume_operator = true;

  EllipticBase(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
               const size_t over_integrate = 0)
    : diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , volume_operator_(over_integrate, diffusion_factor_, diffusion_tensor_)
  {
  }

  EllipticBase(const ThisType& other) = default;

private:
  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;

public:
  const VolumeOperatorType volume_operator_;
}; // class EllipticBase


/**
 * \brief Base class for all simplified elliptic products.
 *
 *        In contrast to the generic variant (s.a.) only a diffusion has to be provided.
 *
 * \note  One could argue to drop this specialization since the L\sa ocalEvaluation::Elliptic already implements an
 *        appropriate specialization and we could just forward to it. We keep this class however because it allows us to
 *        manually have different ctor signatures here and the overhead seems acceptable.
 */
template <class DiffusionImp, class GV, class FieldImp>
class EllipticBase<DiffusionImp, GV, FieldImp, void> : public LocalOperatorProviderBase<GV>
{
  static_assert(Stuff::is_localizable_function<DiffusionImp>::value,
                "DiffusionImp has to be derived from Stuff::LocalizableFunctionInterface!");
  typedef EllipticBase<DiffusionImp, GV, FieldImp, void> ThisType;

protected:
  typedef DiffusionImp DiffusionType;

public:
  typedef GV GridViewType;
  typedef FieldImp FieldType;
  typedef LocalOperator::Codim0Integral<LocalEvaluation::Elliptic<DiffusionType>> VolumeOperatorType;

  static const bool has_volume_operator = true;

  EllipticBase(const DiffusionType& diffusion, const size_t over_integrate = 0)
    : diffusion_(diffusion)
    , volume_operator_(over_integrate, diffusion_)
  {
  }

  EllipticBase(const ThisType& other) = default;

private:
  const DiffusionType& diffusion_;

public:
  const VolumeOperatorType volume_operator_;
}; // class EllipticBase


} // namespace internal
} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_ELLIPTIC_INTERNAL_HH
