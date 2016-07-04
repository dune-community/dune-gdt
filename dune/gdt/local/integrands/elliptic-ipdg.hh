// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2013 - 2016)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_ELLIPTIC_IPDG_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_ELLIPTIC_IPDG_HH

#include <tuple>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/densematrix.hh>

#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/functions/interfaces.hh>

#include "elliptic.hh"
#include "interfaces.hh"
#include "sipdg.hh"

namespace Dune {
namespace GDT {

/**
 *  \brief      Contains local integrands for the symmetric weighted interior penalty discontinuous Galerkin (SWIPDG)
 *              discretization.
 *
 *              For the choice of penalization and the role of the user input see Epshteyn, Riviere (2007):
 *              "Estimation of penalty parameters for symmetric interior penalty Galerkin methods"
 *              For the coice of the weighting see Ern, Stephansen, Zunino (2007): "A discontinuous Galerkin method with
 *              weighted averages for advection-diffusion equations with locally small and anisotropic diffusivity"
 */
namespace LocalEllipticIpdgIntegrands {


enum class Method
{
  ipdg,
  nipdg,
  sipdg,
  swipdg,
  swipdg_affine_factor
};


static constexpr Method default_method = Method::swipdg;


// forwards
template <class DiffusionFactorImp, class DiffusionTensorImp = void, Method method = default_method>
class Inner;


template <class DiffusionFactorImp, class DiffusionTensorImp = void, Method method = default_method>
class BoundaryLHS;


template <class DirichletImp, class DiffusionFactorImp, class DiffusionTensorImp = void, Method method = default_method>
class BoundaryRHS;


namespace internal {


template <class DiffusionFactorImp, class DiffusionTensorImp, Method method>
class InnerTraits : public GDT::internal::LocalEllipticIntegrandTraits<DiffusionFactorImp, DiffusionTensorImp>
{
public:
  typedef Inner<DiffusionFactorImp, DiffusionTensorImp, method> derived_type;
};


template <class DiffusionFactorImp, class DiffusionTensorImp, Method method>
class BoundaryLHSTraits : public GDT::internal::LocalEllipticIntegrandTraits<DiffusionFactorImp, DiffusionTensorImp>
{
public:
  typedef BoundaryLHS<DiffusionFactorImp, DiffusionTensorImp, method> derived_type;
};


template <class DirichletImp, class DiffusionFactorImp, class DiffusionTensorImp, Method method>
class BoundaryRHSTraits
{
  static_assert(Stuff::is_localizable_function<DirichletImp>::value, "DirichletImp has to be a localizable function!");
  typedef GDT::internal::LocalEllipticIntegrandTraits<DiffusionFactorImp, DiffusionTensorImp> EllipticType;

public:
  typedef BoundaryRHS<DirichletImp, DiffusionFactorImp, DiffusionTensorImp, method> derived_type;
  typedef DirichletImp DirichletType;
  typedef typename EllipticType::DiffusionFactorType DiffusionFactorType;
  typedef typename EllipticType::DiffusionTensorType DiffusionTensorType;
  typedef std::tuple<std::shared_ptr<typename DirichletType::LocalfunctionType>,
                     std::shared_ptr<typename DiffusionFactorType::LocalfunctionType>,
                     std::shared_ptr<typename DiffusionTensorType::LocalfunctionType>>
      LocalfunctionTupleType;
  typedef typename EllipticType::EntityType EntityType;
  typedef typename EllipticType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = EllipticType::dimDomain;
}; // class BoundaryRHSTraits


} // namespace internal


/**
 *  see Epshteyn, Riviere, 2007 for the meaning of beta
 * \note The FieldVector< R, dimDomain - 1 > type for the intersection will probably fail for dimDomain 1
 */
template <class DiffusionFactorImp, class DiffusionTensorImp, Method method>
class Inner
    : public LocalFaceIntegrandInterface<internal::InnerTraits<DiffusionFactorImp, DiffusionTensorImp, method>, 4>
{
  typedef LocalEllipticIntegrand<DiffusionFactorImp, DiffusionTensorImp> EllipticType;
  typedef Inner<DiffusionFactorImp, DiffusionTensorImp, method> ThisType;

public:
  typedef internal::InnerTraits<DiffusionFactorImp, DiffusionTensorImp, method> Traits;
  typedef typename Traits::DiffusionFactorType DiffusionFactorType;
  typedef typename Traits::DiffusionTensorType DiffusionTensorType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const size_t dimDomain = Traits::dimDomain;

  Inner(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
        const double beta = LocalSipdgIntegrands::internal::default_beta(dimDomain))
    : elliptic_(diffusion_factor, diffusion_tensor)
    , beta_(beta)
  {
  }

  Inner(const DiffusionFactorType& diffusion_factor,
        const double beta = LocalSipdgIntegrands::internal::default_beta(dimDomain))
    : elliptic_(diffusion_factor)
    , beta_(beta)
  {
  }

  template <
      typename DiffusionType // This disables the ctor if dimDomain == 1, since factor and tensor are then identical
      ,
      typename = typename std::enable_if<(std::is_same<DiffusionType, DiffusionTensorType>::value) // and the ctors
                                         && (dimDomain > 1) && sizeof(DiffusionType)>::type> // ambiguous.
  Inner(const DiffusionType& diffusion, const double beta = LocalSipdgIntegrands::internal::default_beta(dimDomain))
    : elliptic_(diffusion)
    , beta_(beta)
  {
  }

  Inner(const ThisType& other) = default;
  Inner(ThisType&& source)     = default;

  /// \name Required by LocalFaceIntegrandInterface< ..., 4 >.
  /// \{

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(elliptic_.diffusion_factor().local_function(entity),
                           elliptic_.diffusion_tensor().local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t order(
      const LocalfunctionTupleType& local_functions_en, const LocalfunctionTupleType& local_functions_ne,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base_en,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base_en,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base_ne,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base_ne) const
  {
    const auto local_diffusion_factor_en = std::get<0>(local_functions_en);
    const auto local_diffusion_tensor_en = std::get<1>(local_functions_en);
    const auto local_diffusion_factor_ne = std::get<0>(local_functions_ne);
    const auto local_diffusion_tensor_ne = std::get<1>(local_functions_ne);
    return order(*local_diffusion_factor_en,
                 *local_diffusion_tensor_en,
                 *local_diffusion_factor_ne,
                 *local_diffusion_tensor_ne,
                 test_base_en,
                 ansatz_base_en,
                 test_base_ne,
                 ansatz_base_ne);
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class IntersectionType, class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  void
  evaluate(const LocalfunctionTupleType& local_functions_en, const LocalfunctionTupleType& local_functions_ne,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base_en,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base_en,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base_ne,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base_ne,
           const IntersectionType& intersection, const Dune::FieldVector<DomainFieldType, dimDomain - 1>& local_point,
           Dune::DynamicMatrix<R>& ret_en_en, Dune::DynamicMatrix<R>& ret_ne_ne, Dune::DynamicMatrix<R>& ret_en_ne,
           Dune::DynamicMatrix<R>& ret_ne_en) const
  {
    const auto local_diffusion_factor_en = std::get<0>(local_functions_en);
    const auto local_diffusion_tensor_en = std::get<1>(local_functions_en);
    const auto local_diffusion_factor_ne = std::get<0>(local_functions_ne);
    const auto local_diffusion_tensor_ne = std::get<1>(local_functions_ne);
    evaluate(*local_diffusion_factor_en,
             *local_diffusion_tensor_en,
             *local_diffusion_factor_ne,
             *local_diffusion_tensor_ne,
             test_base_en,
             ansatz_base_en,
             test_base_ne,
             ansatz_base_ne,
             intersection,
             local_point,
             ret_en_en,
             ret_ne_ne,
             ret_en_ne,
             ret_ne_en);
  }

  /// \}
private:
  // The Helper struct and private order/evaluate methods are required to provide varaints of order and evaluate for the
  // single diffusion case.

  template <bool single_diffusion, bool is_factor, class Anyone = void>
  struct Helper
  {
    static_assert(AlwaysFalse<Anyone>::value, "These variants of order and evaluate are only available for the single "
                                              "diffusion case (i.e., if DiffusionTensorImp is void)!");
  };

  template <class Anyone>
  struct Helper<true, true, Anyone>
  {
    template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
    static size_t
    order(const ThisType& ths,
          const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1>& local_diffusion_factor_en,
          const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1>& local_diffusion_factor_ne,
          const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base_en,
          const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base_en,
          const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base_ne,
          const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base_ne)
    {
      const auto local_functions_en        = ths.localFunctions(local_diffusion_factor_en.entity());
      const auto local_functions_ne        = ths.localFunctions(local_diffusion_factor_ne.entity());
      const auto local_diffusion_tensor_en = std::get<1>(local_functions_en);
      const auto local_diffusion_tensor_ne = std::get<1>(local_functions_ne);
      return ths.order(local_diffusion_factor_en,
                       *local_diffusion_tensor_en,
                       local_diffusion_factor_ne,
                       *local_diffusion_tensor_ne,
                       test_base_en,
                       ansatz_base_en,
                       test_base_ne,
                       ansatz_base_ne);
    } // size_t order(...)

    template <class R, class IntersectionType>
    static void evaluate(
        const ThisType& ths,
        const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& local_diffusion_factor_en,
        const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& local_diffusion_factor_ne,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& test_base_en,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& ansatz_base_en,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& test_base_ne,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& ansatz_base_ne,
        const IntersectionType& intersection, const Dune::FieldVector<DomainFieldType, dimDomain - 1>& local_point,
        Dune::DynamicMatrix<R>& ret_en_en, Dune::DynamicMatrix<R>& ret_ne_ne, Dune::DynamicMatrix<R>& ret_en_ne,
        Dune::DynamicMatrix<R>& ret_ne_en)
    {
      const auto local_functions_en        = ths.localFunctions(local_diffusion_factor_en.entity());
      const auto local_functions_ne        = ths.localFunctions(local_diffusion_factor_ne.entity());
      const auto local_diffusion_tensor_en = std::get<1>(local_functions_en);
      const auto local_diffusion_tensor_ne = std::get<1>(local_functions_ne);
      ths.evaluate(local_diffusion_factor_en,
                   *local_diffusion_tensor_en,
                   local_diffusion_factor_ne,
                   *local_diffusion_tensor_ne,
                   test_base_en,
                   ansatz_base_en,
                   test_base_ne,
                   ansatz_base_ne,
                   intersection,
                   local_point,
                   ret_en_en,
                   ret_ne_ne,
                   ret_en_ne,
                   ret_ne_en);
    }
  }; // struct Helper< true, true, ... >

  template <class Anyone>
  struct Helper<true, false, Anyone>
  {
    template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
    static size_t
    order(const ThisType& ths, const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain,
                                                                   dimDomain>& local_diffusion_tensor_en,
          const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain, dimDomain>&
              local_diffusion_tensor_ne,
          const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base_en,
          const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base_en,
          const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base_ne,
          const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base_ne)
    {
      const auto local_functions_en        = ths.localFunctions(local_diffusion_tensor_en.entity());
      const auto local_functions_ne        = ths.localFunctions(local_diffusion_tensor_ne.entity());
      const auto local_diffusion_factor_en = std::get<0>(local_functions_en);
      const auto local_diffusion_factor_ne = std::get<0>(local_functions_ne);
      return ths.order(*local_diffusion_factor_en,
                       local_diffusion_tensor_en,
                       *local_diffusion_factor_ne,
                       local_diffusion_tensor_ne,
                       test_base_en,
                       ansatz_base_en,
                       test_base_ne,
                       ansatz_base_ne);
    } // size_t order(...)

    template <class R, class IntersectionType>
    static void
    evaluate(const ThisType& ths, const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R,
                                                                      dimDomain, dimDomain>& local_diffusion_tensor_en,
             const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain, dimDomain>&
                 local_diffusion_tensor_ne,
             const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& test_base_en,
             const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& ansatz_base_en,
             const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& test_base_ne,
             const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& ansatz_base_ne,
             const IntersectionType& intersection, const Dune::FieldVector<DomainFieldType, dimDomain - 1>& local_point,
             Dune::DynamicMatrix<R>& ret_en_en, Dune::DynamicMatrix<R>& ret_ne_ne, Dune::DynamicMatrix<R>& ret_en_ne,
             Dune::DynamicMatrix<R>& ret_ne_en)
    {
      const auto local_functions_en        = ths.localFunctions(local_diffusion_tensor_en.entity());
      const auto local_functions_ne        = ths.localFunctions(local_diffusion_tensor_ne.entity());
      const auto local_diffusion_factor_en = std::get<0>(local_functions_en);
      const auto local_diffusion_factor_ne = std::get<0>(local_functions_ne);
      ths.evaluate(*local_diffusion_factor_en,
                   local_diffusion_tensor_en,
                   *local_diffusion_factor_ne,
                   local_diffusion_tensor_ne,
                   test_base_en,
                   ansatz_base_en,
                   test_base_ne,
                   ansatz_base_ne,
                   intersection,
                   local_point,
                   ret_en_en,
                   ret_ne_ne,
                   ret_en_ne,
                   ret_ne_en);
    }
  }; // struct Helper< true, false, ... >

public:
  /// \name Redirects for single diffusion (either factor or tensor, but not both).
  /// \{

  template <class R, size_t rD, size_t rCD, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t order(
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rD, rCD>& local_diffusion_en,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rD, rCD>& local_diffusion_ne,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base_en,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base_en,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base_ne,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base_ne) const
  {
    return Helper < std::is_same<DiffusionTensorImp, void>::value, rD == 1
                                                                       && rCD
                                                                              == 1 > ::order(*this,
                                                                                             local_diffusion_en,
                                                                                             local_diffusion_ne,
                                                                                             test_base_en,
                                                                                             ansatz_base_en,
                                                                                             test_base_ne,
                                                                                             ansatz_base_ne);
  } // ... order(...)

  template <class R, size_t rD, size_t rCD, size_t rT, size_t rCT, size_t rA, size_t rCA, class IntersectionType>
  void
  evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rD, rCD>& local_diffusion_en,
           const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rD, rCD>& local_diffusion_ne,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base_en,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base_en,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base_ne,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base_ne,
           const IntersectionType& intersection, const Dune::FieldVector<DomainFieldType, dimDomain - 1>& local_point,
           Dune::DynamicMatrix<R>& ret_en_en, Dune::DynamicMatrix<R>& ret_ne_ne, Dune::DynamicMatrix<R>& ret_en_ne,
           Dune::DynamicMatrix<R>& ret_ne_en) const
  {
    Helper<std::is_same<DiffusionTensorImp, void>::value, rD == 1 && rCD == 1>::evaluate(*this,
                                                                                         local_diffusion_en,
                                                                                         local_diffusion_ne,
                                                                                         test_base_en,
                                                                                         ansatz_base_en,
                                                                                         test_base_ne,
                                                                                         ansatz_base_ne,
                                                                                         intersection,
                                                                                         local_point,
                                                                                         ret_en_en,
                                                                                         ret_ne_ne,
                                                                                         ret_en_ne,
                                                                                         ret_ne_en);
  } // ... evaluate(...)

  /// \}
  /// \name Actual Implementation of order.
  /// \{

  template <class R, size_t rDF, size_t rCDF, size_t rDT, size_t rCDT, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t order(
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDF, rCDF>&
          local_diffusion_factor_en,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDT, rCDT>&
          local_diffusion_tensor_en,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDF, rCDF>&
          local_diffusion_factor_ne,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDT, rCDT>&
          local_diffusion_tensor_ne,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base_en,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base_en,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base_ne,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base_ne) const
  {
    return std::max(local_diffusion_factor_en.order(), local_diffusion_factor_ne.order())
           + std::max(local_diffusion_tensor_en.order(), local_diffusion_tensor_ne.order())
           + std::max(test_base_en.order(), test_base_ne.order())
           + std::max(ansatz_base_en.order(), ansatz_base_ne.order());
  } // size_t order(...)

  /// \}
  /// \name Actual Implementation of evaluate.
  /// \{

  template <class R, class IntersectionType>
  void evaluate(
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& local_diffusion_factor_en,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain, dimDomain>&
          local_diffusion_tensor_en,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& local_diffusion_factor_ne,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain, dimDomain>&
          local_diffusion_tensor_ne,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& test_base_en,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& ansatz_base_en,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& test_base_ne,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& ansatz_base_ne,
      const IntersectionType& intersection, const Dune::FieldVector<DomainFieldType, dimDomain - 1>& local_point,
      Dune::DynamicMatrix<R>& ret_en_en, Dune::DynamicMatrix<R>& ret_ne_ne, Dune::DynamicMatrix<R>& ret_en_ne,
      Dune::DynamicMatrix<R>& ret_ne_en) const
  {
    typedef Stuff::Common::FieldMatrix<R, dimDomain, dimDomain> TensorType;
    // clear ret
    ret_en_en *= 0.0;
    ret_ne_ne *= 0.0;
    ret_en_ne *= 0.0;
    ret_ne_en *= 0.0;
    // convert local point (which is in intersection coordinates) to entity/neighbor coordinates
    const auto local_point_en = intersection.geometryInInside().global(local_point);
    const auto local_point_ne = intersection.geometryInOutside().global(local_point);
    const auto normal         = intersection.unitOuterNormal(local_point);
    // evaluate local function
    const auto local_diffusion_factor_value_en       = local_diffusion_factor_en.evaluate(local_point_en);
    const TensorType local_diffusion_tensor_value_en = local_diffusion_tensor_en.evaluate(local_point_en);
    const auto local_diffusion_factor_value_ne       = local_diffusion_factor_ne.evaluate(local_point_ne);
    const TensorType local_diffusion_tensor_value_ne = local_diffusion_tensor_ne.evaluate(local_point_ne);
    const auto diffusion_value_en = local_diffusion_tensor_value_en * local_diffusion_factor_value_en;
    const auto diffusion_value_ne = local_diffusion_tensor_value_ne * local_diffusion_factor_value_ne;
    //    // this evaluation has to be linear wrt the diffusion factor, so no other averaging method is allowed here!
    //    const auto local_diffusion_factor = (local_diffusion_factor_en + local_diffusion_factor_ne) * 0.5;
    // compute penalty factor (see Epshteyn, Riviere, 2007)
    const size_t max_polorder = std::max(
        test_base_en.order(), std::max(ansatz_base_en.order(), std::max(test_base_ne.order(), ansatz_base_ne.order())));
    const R sigma = LocalSipdgIntegrands::internal::inner_sigma(max_polorder);
    // compute weighting (see Ern, Stephansen, Zunino 2007)
    const R delta_plus  = normal * (/*local_diffusion_tensor_ne*/ diffusion_value_ne * normal);
    const R delta_minus = normal * (/*local_diffusion_tensor_en*/ diffusion_value_en * normal);
    const R gamma       = (delta_plus * delta_minus) / (delta_plus + delta_minus);
    const R penalty = (/*local_diffusion_factor **/ sigma * gamma) / std::pow(intersection.geometry().volume(), beta_);
    const R weight_plus  = delta_minus / (delta_plus + delta_minus);
    const R weight_minus = delta_plus / (delta_plus + delta_minus);
    // evaluate bases
    // * entity
    //   * test
    const size_t rows_en         = test_base_en.size();
    const auto test_values_en    = test_base_en.evaluate(local_point_en);
    const auto test_gradients_en = test_base_en.jacobian(local_point_en);
    //   * ansatz
    const size_t cols_en           = ansatz_base_en.size();
    const auto ansatz_values_en    = ansatz_base_en.evaluate(local_point_en);
    const auto ansatz_gradients_en = ansatz_base_en.jacobian(local_point_en);
    // * neighbor
    //   * test
    const size_t rows_ne         = test_base_ne.size();
    const auto test_values_ne    = test_base_ne.evaluate(local_point_ne);
    const auto test_gradients_ne = test_base_ne.jacobian(local_point_ne);
    //   * ansatz
    const size_t cols_ne           = ansatz_base_ne.size();
    const auto ansatz_values_ne    = ansatz_base_ne.evaluate(local_point_ne);
    const auto ansatz_gradients_ne = ansatz_base_ne.jacobian(local_point_ne);
    // compute the evaluations
    assert(ret_en_en.rows() >= rows_en && ret_en_en.cols() >= cols_en);
    assert(ret_en_ne.rows() >= rows_en && ret_en_ne.cols() >= cols_ne);
    assert(ret_ne_en.rows() >= rows_ne && ret_ne_en.cols() >= cols_en);
    assert(ret_ne_ne.rows() >= rows_ne && ret_ne_ne.cols() >= cols_ne);
    // loop over all entity test basis functions
    for (size_t ii = 0; ii < rows_en; ++ii) {
      auto& ret_en_en_row = ret_en_en[ii];
      auto& ret_en_ne_row = ret_en_ne[ii];
      // loop over all entity ansatz basis functions
      for (size_t jj = 0; jj < cols_en; ++jj) {
        // consistency term
        ret_en_en_row[jj] +=
            -weight_minus * ((diffusion_value_en * ansatz_gradients_en[jj][0]) * normal) * test_values_en[ii];
        // symmetry term
        ret_en_en_row[jj] +=
            -weight_minus * ansatz_values_en[jj] * ((diffusion_value_en * test_gradients_en[ii][0]) * normal);
        // penalty term
        ret_en_en_row[jj] += penalty * ansatz_values_en[jj] * test_values_en[ii];
      } // loop over all entity ansatz basis functions
      // loop over all neighbor ansatz basis functions
      for (size_t jj = 0; jj < cols_ne; ++jj) {
        // consistency term
        ret_en_ne_row[jj] +=
            -weight_plus * ((diffusion_value_ne * ansatz_gradients_ne[jj][0]) * normal) * test_values_en[ii];
        // symmetry term
        ret_en_ne_row[jj] +=
            weight_minus * ansatz_values_ne[jj] * ((diffusion_value_en * test_gradients_en[ii][0]) * normal);
        // penalty term
        ret_en_ne_row[jj] += -1.0 * penalty * ansatz_values_ne[jj] * test_values_en[ii];
      } // loop over all neighbor ansatz basis functions
    } // loop over all entity test basis functions
    // loop over all neighbor test basis functions
    for (size_t ii = 0; ii < rows_ne; ++ii) {
      auto& ret_ne_en_row = ret_ne_en[ii];
      auto& ret_ne_ne_row = ret_ne_ne[ii];
      // loop over all entity ansatz basis functions
      for (size_t jj = 0; jj < cols_en; ++jj) {
        // consistency term
        ret_ne_en_row[jj] +=
            weight_minus * ((diffusion_value_en * ansatz_gradients_en[jj][0]) * normal) * test_values_ne[ii];
        // symmetry term
        ret_ne_en_row[jj] +=
            -weight_plus * ansatz_values_en[jj] * ((diffusion_value_ne * test_gradients_ne[ii][0]) * normal);
        // penalty term
        ret_ne_en_row[jj] += -1.0 * penalty * ansatz_values_en[jj] * test_values_ne[ii];
      } // loop over all entity ansatz basis functions
      // loop over all neighbor ansatz basis functions
      for (size_t jj = 0; jj < cols_ne; ++jj) {
        // consistency term
        ret_ne_ne_row[jj] +=
            weight_plus * ((diffusion_value_ne * ansatz_gradients_ne[jj][0]) * normal) * test_values_ne[ii];
        // symmetry term
        ret_ne_ne_row[jj] +=
            weight_plus * ansatz_values_ne[jj] * ((diffusion_value_ne * test_gradients_ne[ii][0]) * normal);
        // penalty term
        ret_ne_ne_row[jj] += penalty * ansatz_values_ne[jj] * test_values_ne[ii];
      } // loop over all neighbor ansatz basis functions
    } // loop over all neighbor test basis functions
  } // ... evaluate(...)

  /// \}

private:
  const EllipticType elliptic_;
  const double beta_;
}; // Inner


template <class DiffusionFactorImp, class DiffusionTensorImp, Method method>
class BoundaryLHS
    : public LocalFaceIntegrandInterface<internal::BoundaryLHSTraits<DiffusionFactorImp, DiffusionTensorImp, method>, 2>
{
public:
  typedef LocalEllipticIntegrand<DiffusionFactorImp, DiffusionTensorImp> EllipticType;
  typedef BoundaryLHS<DiffusionFactorImp, DiffusionTensorImp, method> ThisType;

public:
  typedef internal::BoundaryLHSTraits<DiffusionFactorImp, DiffusionTensorImp, method> Traits;
  typedef typename Traits::DiffusionFactorType DiffusionFactorType;
  typedef typename Traits::DiffusionTensorType DiffusionTensorType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const size_t dimDomain = Traits::dimDomain;

  BoundaryLHS(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
              const double beta = LocalSipdgIntegrands::internal::default_beta(dimDomain))
    : elliptic_(diffusion_factor, diffusion_tensor)
    , beta_(beta)
  {
  }

  BoundaryLHS(const DiffusionFactorType& diffusion_factor,
              const double beta = LocalSipdgIntegrands::internal::default_beta(dimDomain))
    : elliptic_(diffusion_factor)
    , beta_(beta)
  {
  }

  template <
      typename DiffusionType // This disables the ctor if dimDomain == 1, since factor and tensor are then identical
      ,
      typename = typename std::enable_if<(std::is_same<DiffusionType, DiffusionTensorType>::value) // and the ctors
                                         && (dimDomain > 1) && sizeof(DiffusionType)>::type> // ambiguous.
  BoundaryLHS(const DiffusionType& diffusion,
              const double beta = LocalSipdgIntegrands::internal::default_beta(dimDomain))
    : elliptic_(diffusion)
    , beta_(beta)
  {
  }

  BoundaryLHS(const ThisType& other) = default;
  BoundaryLHS(ThisType&& source)     = default;

  /// \name Required by LocalFaceIntegrandInterface< ..., 2 >

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(elliptic_.diffusion_factor().local_function(entity),
                           elliptic_.diffusion_tensor().local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t
  order(const LocalfunctionTupleType& local_functions,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base) const
  {
    const auto local_diffusion_factor = std::get<0>(local_functions);
    const auto local_diffusion_tensor = std::get<1>(local_functions);
    return order(*local_diffusion_factor, *local_diffusion_tensor, test_base, ansatz_base);
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class IntersectionType, class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  void evaluate(const LocalfunctionTupleType& local_functions,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base,
                const IntersectionType& intersection,
                const Dune::FieldVector<DomainFieldType, dimDomain - 1>& local_point, Dune::DynamicMatrix<R>& ret) const
  {
    const auto local_diffusion_factor = std::get<0>(local_functions);
    const auto local_diffusion_tensor = std::get<1>(local_functions);
    evaluate(*local_diffusion_factor, *local_diffusion_tensor, test_base, ansatz_base, intersection, local_point, ret);
  }

  /// \}
private:
  // The Helper struct and private order/evaluate methods are required to provide varaints of order and evaluate for the
  // single diffusion case.

  template <bool single_diffusion, bool is_factor, class Anyone = void>
  struct Helper
  {
    static_assert(AlwaysFalse<Anyone>::value, "These variants of order and evaluate are only available for the single "
                                              "diffusion case (i.e., if DiffusionTensorImp is void)!");
  };

  template <class Anyone>
  struct Helper<true, true, Anyone>
  {
    template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
    static size_t
    order(const ThisType ths,
          const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1>& local_diffusion_factor,
          const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base,
          const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base)
    {
      const auto local_functions        = ths.localFunctions(local_diffusion_factor.entity());
      const auto local_diffusion_tensor = std::get<1>(local_functions);
      return ths.order(local_diffusion_factor, *local_diffusion_tensor, test_base, ansatz_base);
    }

    template <class R, class IntersectionType>
    static void evaluate(
        const ThisType& ths,
        const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& local_diffusion_factor,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& test_base,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& ansatz_base,
        const IntersectionType& intersection, const Dune::FieldVector<DomainFieldType, dimDomain - 1>& local_point,
        Dune::DynamicMatrix<R>& ret)
    {
      const auto local_functions        = ths.localFunctions(local_diffusion_factor.entity());
      const auto local_diffusion_tensor = std::get<1>(local_functions);
      ths.evaluate(
          local_diffusion_factor, *local_diffusion_tensor, test_base, ansatz_base, intersection, local_point, ret);
    }
  }; // struct Helper< true, true, ... >

  template <class Anyone>
  struct Helper<true, false, Anyone>
  {
    template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
    static size_t
    order(const ThisType& ths, const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain,
                                                                   dimDomain>& local_diffusion_tensor,
          const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base,
          const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base)
    {
      const auto local_functions        = ths.localFunctions(local_diffusion_tensor.entity());
      const auto local_diffusion_factor = std::get<0>(local_functions);
      return ths.order(*local_diffusion_factor, local_diffusion_tensor, test_base, ansatz_base);
    }

    template <class R, class IntersectionType>
    static void
    evaluate(const ThisType& ths, const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R,
                                                                      dimDomain, dimDomain>& local_diffusion_tensor,
             const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& test_base,
             const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& ansatz_base,
             const IntersectionType& intersection, const Dune::FieldVector<DomainFieldType, dimDomain - 1>& local_point,
             Dune::DynamicMatrix<R>& ret)
    {
      const auto local_functions        = ths.localFunctions(local_diffusion_tensor.entity());
      const auto local_diffusion_factor = std::get<0>(local_functions);
      ths.evaluate(
          *local_diffusion_factor, local_diffusion_tensor, test_base, ansatz_base, intersection, local_point, ret);
    }
  }; // struct Helper< true, false, ... >

public:
  /// \name Redirects for single diffusion (either factor or tensor, but not both).
  /// \{

  template <class R, size_t rD, size_t rCD, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t
  order(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rD, rCD>& local_diffusion,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base) const
  {
    return Helper < std::is_same<DiffusionTensorImp, void>::value,
           rD == 1 && rCD == 1 > ::order(*this, local_diffusion, test_base, ansatz_base);
  }

  template <class R, size_t rD, size_t rCD, size_t rT, size_t rCT, size_t rA, size_t rCA, class IntersectionType>
  void
  evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rD, rCD>& local_diffusion,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base,
           const IntersectionType& intersection, const Dune::FieldVector<DomainFieldType, dimDomain - 1>& local_point,
           Dune::DynamicMatrix<R>& ret) const
  {
    Helper<std::is_same<DiffusionTensorImp, void>::value, rD == 1 && rCD == 1>::evaluate(
        *this, local_diffusion, test_base, ansatz_base, intersection, local_point, ret);
  }

  /// \}
  /// \name Actual implementation of order.
  /// \{

  template <class R, size_t rDF, size_t rCDF, size_t rDT, size_t rCDT, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t order(
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDF, rCDF>& local_diffusion_factor,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDT, rCDT>& local_diffusion_tensor,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base) const
  {
    return local_diffusion_factor.order() + local_diffusion_tensor.order() + test_base.order() + ansatz_base.order();
  }

  /// \}
  /// \name Actual implementation of evaluate.
  /// \{

  template <class R, class IntersectionType>
  void
  evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& local_diffusion_factor,
           const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain, dimDomain>&
               local_diffusion_tensor,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& test_base,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& ansatz_base,
           const IntersectionType& intersection, const Dune::FieldVector<DomainFieldType, dimDomain - 1>& local_point,
           Dune::DynamicMatrix<R>& ret) const
  {
    typedef Stuff::Common::FieldMatrix<R, dimDomain, dimDomain> TensorType;
    // clear ret
    ret *= 0.0;
    // get local point (which is in intersection coordinates) in entity coordinates
    const auto local_point_entity = intersection.geometryInInside().global(local_point);
    const auto normal             = intersection.unitOuterNormal(local_point);
    // evaluate local function
    const auto diffusion_factor_value       = local_diffusion_factor.evaluate(local_point_entity);
    const TensorType diffusion_tensor_value = local_diffusion_tensor.evaluate(local_point_entity);
    const auto diffusion_value              = diffusion_tensor_value * diffusion_factor_value;
    // compute penalty (see Epshteyn, Riviere, 2007)
    const size_t max_polorder = std::max(test_base.order(), ansatz_base.order());
    const R sigma             = LocalSipdgIntegrands::internal::boundary_sigma(max_polorder);
    // compute weighting (see Ern, Stephansen, Zunino 2007)
    const R gamma   = normal * (/*diffusion_tensor_value*/ diffusion_value * normal);
    const R penalty = (/*diffusion_factor_value **/ sigma * gamma) / std::pow(intersection.geometry().volume(), beta_);
    // evaluate bases
    // * test
    const size_t rows         = test_base.size();
    const auto test_values    = test_base.evaluate(local_point_entity);
    const auto test_gradients = test_base.jacobian(local_point_entity);
    // * ansatz
    const size_t cols           = ansatz_base.size();
    const auto ansatz_values    = ansatz_base.evaluate(local_point_entity);
    const auto ansatz_gradients = ansatz_base.jacobian(local_point_entity);
    // compute products
    assert(ret.rows() >= rows && ret.cols() >= cols);
    // loop over all test basis functions
    for (size_t ii = 0; ii < rows; ++ii) {
      auto& retRow = ret[ii];
      // loop over all ansatz basis functions
      for (size_t jj = 0; jj < cols; ++jj) {
        // consistency term
        retRow[jj] += -1.0 * ((diffusion_value * ansatz_gradients[jj][0]) * normal) * test_values[ii];
        // symmetry term
        retRow[jj] += -1.0 * ansatz_values[jj] * ((diffusion_value * test_gradients[ii][0]) * normal);
        // penalty term
        retRow[jj] += penalty * ansatz_values[jj] * test_values[ii];
      } // loop over all ansatz basis functions
    } // loop over all test basis functions
  } // void evaluate(...)

  /// \}

private:
  const EllipticType elliptic_;
  const double beta_;
}; // class BoundaryLHS


template <class DirichletImp, class DiffusionFactorImp, class DiffusionTensorImp, Method method>
class BoundaryRHS : public LocalFaceIntegrandInterface<internal::BoundaryRHSTraits<DirichletImp, DiffusionFactorImp,
                                                                                   DiffusionTensorImp, method>,
                                                       1>
{
  typedef LocalEllipticIntegrand<DiffusionFactorImp, DiffusionTensorImp> EllipticType;
  typedef BoundaryRHS<DirichletImp, DiffusionFactorImp, DiffusionTensorImp, method> ThisType;

public:
  typedef internal::BoundaryRHSTraits<DirichletImp, DiffusionFactorImp, DiffusionTensorImp, method> Traits;
  typedef typename Traits::DirichletType DirichletType;
  typedef typename Traits::DiffusionFactorType DiffusionFactorType;
  typedef typename Traits::DiffusionTensorType DiffusionTensorType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const size_t dimDomain = Traits::dimDomain;

  BoundaryRHS(const DirichletType& dirichlet, const DiffusionFactorType& diffusion_factor,
              const DiffusionTensorType& diffusion_tensor,
              const double beta = LocalSipdgIntegrands::internal::default_beta(dimDomain))
    : dirichlet_(dirichlet)
    , elliptic_(diffusion_factor, diffusion_tensor)
    , beta_(beta)
  {
  }

  BoundaryRHS(const DirichletType& dirichlet, const DiffusionFactorType& diffusion_factor,
              const double beta = LocalSipdgIntegrands::internal::default_beta(dimDomain))
    : dirichlet_(dirichlet)
    , elliptic_(diffusion_factor)
    , beta_(beta)
  {
  }

  template <
      typename DiffusionType // This disables the ctor if dimDomain == 1, since factor and tensor are then identical
      ,
      typename = typename std::enable_if<(std::is_same<DiffusionType, DiffusionTensorType>::value) // and the ctors
                                         && (dimDomain > 1) && sizeof(DiffusionType)>::type> // ambiguous.
  BoundaryRHS(const DirichletType& dirichlet, const DiffusionType& diffusion,
              const double beta = LocalSipdgIntegrands::internal::default_beta(dimDomain))
    : dirichlet_(dirichlet)
    , elliptic_(diffusion)
    , beta_(beta)
  {
  }

  BoundaryRHS(const ThisType& other) = default;
  BoundaryRHS(ThisType&& source)     = default;

  /// \name Required by LocalFaceIntegrandInterface< ..., 1 >.
  /// \{

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(dirichlet_.local_function(entity),
                           elliptic_.diffusion_factor().local_function(entity),
                           elliptic_.diffusion_tensor().local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class R, size_t r, size_t rC>
  size_t
  order(const LocalfunctionTupleType& local_functions,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, rC>& test_base) const
  {
    const auto local_dirichlet        = std::get<0>(local_functions);
    const auto local_diffusion_factor = std::get<1>(local_functions);
    const auto local_diffusion_tensor = std::get<2>(local_functions);
    return order(*local_dirichlet, *local_diffusion_factor, *local_diffusion_tensor, test_base);
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class IntersectionType, class R, size_t r, size_t rC>
  void evaluate(const LocalfunctionTupleType& local_functions,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, r, rC>& test_base,
                const IntersectionType& intersection,
                const Dune::FieldVector<DomainFieldType, dimDomain - 1>& local_point, Dune::DynamicVector<R>& ret) const
  {
    const auto local_dirichlet        = std::get<0>(local_functions);
    const auto local_diffusion_factor = std::get<1>(local_functions);
    const auto local_diffusion_tensor = std::get<2>(local_functions);
    evaluate(
        *local_dirichlet, *local_diffusion_factor, *local_diffusion_tensor, test_base, intersection, local_point, ret);
  }

  /// \}
  /// \name Actual implementation of order.
  /// \{

  template <class R, size_t rDF, size_t rCDF, size_t rDT, size_t rCDT, size_t rLR, size_t rCLR, size_t rT, size_t rCT>
  size_t order(
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rLR, rCLR>& local_dirichlet,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDF, rCDF>& local_diffusion_factor,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rDT, rCDT>& local_diffusion_tensor,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base) const
  {
    const size_t test_order          = test_base.order();
    const size_t test_gradient_order = std::max(ssize_t(test_order) - 1, ssize_t(0));
    const size_t diffusionOrder      = local_diffusion_factor.order() + local_diffusion_tensor.order();
    const size_t dirichletOrder      = local_dirichlet.order();
    return std::max(test_order + dirichletOrder, diffusionOrder + test_gradient_order + dirichletOrder);
  } // ... order(...)

  /// \}
  /// \name Actual implementation of evaluate.
  /// \{

  template <class R, class IntersectionType>
  void
  evaluate(const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& local_dirichlet,
           const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& local_diffusion_factor,
           const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain, dimDomain>&
               local_diffusion_tensor,
           const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& test_base,
           const IntersectionType& intersection, const Dune::FieldVector<DomainFieldType, dimDomain - 1>& local_point,
           Dune::DynamicVector<R>& ret) const
  {
    typedef Stuff::Common::FieldMatrix<R, dimDomain, dimDomain> TensorType;
    // clear ret
    ret *= 0.0;
    // get local point (which is in intersection coordinates) in entity coordinates
    const auto local_point_entity = intersection.geometryInInside().global(local_point);
    const auto normal             = intersection.unitOuterNormal(local_point);
    // evaluate local functions
    const auto dirichlet_value              = local_dirichlet.evaluate(local_point_entity);
    const auto diffusion_factor_value       = local_diffusion_factor.evaluate(local_point_entity);
    const TensorType diffusion_tensor_value = local_diffusion_tensor.evaluate(local_point_entity);
    // compute penalty (see Epshteyn, Riviere, 2007)
    const size_t polorder = test_base.order();
    const R sigma         = LocalSipdgIntegrands::internal::boundary_sigma(polorder);
    // compute weighting (see Ern, Stephansen, Zunino 2007)
    const R gamma   = normal * (diffusion_tensor_value * normal);
    const R penalty = (diffusion_factor_value * sigma * gamma) / std::pow(intersection.geometry().volume(), beta_);
    // compute diffusion value
    const auto diffusion_value = diffusion_tensor_value * diffusion_factor_value;
    // evaluate basis
    const size_t size         = test_base.size();
    const auto test_values    = test_base.evaluate(local_point_entity);
    const auto test_gradients = test_base.jacobian(local_point_entity);
    // compute
    assert(ret.size() >= size);
    // loop over all test basis functions
    for (size_t ii = 0; ii < size; ++ii) {
      // symmetry term
      ret[ii] += -1.0 * dirichlet_value * ((diffusion_value * test_gradients[ii][0]) * normal);
      // penalty term
      ret[ii] += penalty * dirichlet_value * test_values[ii];
    } // loop over all test basis functions
  } // ... evaluate(...)

  /// \}

  const DirichletType& dirichlet_;
  const EllipticType elliptic_;
  const double beta_;
}; // class BoundaryRHS


} // namespace LocalEllipticIpdgIntegrands
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_ELLIPTIC_IPDG_HH
