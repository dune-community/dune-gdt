// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_L2_INTERNAL_HH
#define DUNE_GDT_PRODUCTS_L2_INTERNAL_HH

#include <type_traits>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/interfaces.hh>

#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localevaluation/product.hh>

#include "base-internal.hh"

namespace Dune {
namespace GDT {
namespace Products {
namespace internal {


/**
 * \brief Base class for all weighted L2 products.
 * \note  Most likely you do not want to use this class directly, but Products::WeightedL2Localizable,
 *        Products::WeightedL2Assemblable or Products::WeightedL2 instead!
 */
template <class GV, class FunctionImp, class FieldImp>
class WeightedL2Base : public LocalOperatorProviderBase<GV>
{
  static_assert(Stuff::is_localizable_function<FunctionImp>::value,
                "FunctionImp has to be derived from Stuff::LocalizableFunctionInterface!");
  typedef WeightedL2Base<GV, FunctionImp, FieldImp> ThisType;

protected:
  typedef FunctionImp FunctionType;

public:
  typedef GV GridViewType;
  typedef FieldImp FieldType;
  typedef LocalOperator::Codim0Integral<LocalEvaluation::Product<FunctionType>> VolumeOperatorType;

  static const bool has_volume_operator = true;

  WeightedL2Base(const FunctionType& function, const size_t over_integrate = 0)
    : function_(function)
    , volume_operator_(over_integrate, function_)
  {
  }

  WeightedL2Base(const ThisType& other) = default;

private:
  const FunctionType& function_;

public:
  const VolumeOperatorType volume_operator_;
}; // WeightedL2Base


/**
 * \brief Base class for all L2 products.
 *
 *        This class is implemented using WeightedL2Base with a Stuff::Functions::Constant of value 1 as the weight.
 * \note  Most likely you do not want to use this class directly, but Products::L2Localizable, Products::L2Assemblable
 *        or Products::L2 instead!
 */
template <class GV, class FieldImp>
class L2Base : DSC::ConstStorageProvider<Stuff::Functions::Constant<typename GV::template Codim<0>::Entity,
                                                                    typename GV::ctype, GV::dimension, FieldImp, 1>>,
               public WeightedL2Base<GV, Stuff::Functions::Constant<typename GV::template Codim<0>::Entity,
                                                                    typename GV::ctype, GV::dimension, FieldImp, 1>,
                                     FieldImp>
{
  typedef DSC::ConstStorageProvider<Stuff::Functions::Constant<
      typename GV::template Codim<0>::Entity, typename GV::ctype, GV::dimension, FieldImp, 1>> StorageBaseType;
  typedef WeightedL2Base<GV, Stuff::Functions::Constant<typename GV::template Codim<0>::Entity, typename GV::ctype,
                                                        GV::dimension, FieldImp, 1>,
                         FieldImp> L2BaseType;
  typedef L2Base<GV, FieldImp> ThisType;

public:
  L2Base(const size_t over_integrate = 0)
    : StorageBaseType(new typename L2BaseType::FunctionType(1))
    , L2BaseType(this->storage_access(), over_integrate)
    , over_integrate_(over_integrate)
  {
  }

  /**
   * \note We need the manual copy ctor bc of the Stuff::Common::ConstStorageProvider
   */
  L2Base(const ThisType& other)
    : StorageBaseType(new typename L2BaseType::FunctionType(1))
    , L2BaseType(this->storage_access(), other.over_integrate_)
    , over_integrate_(other.over_integrate_)
  {
  }

private:
  const size_t over_integrate_; //!< needed to provide manual copy ctor
}; // L2Base


} // namespace internal
} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_L2_INTERNAL_HH
