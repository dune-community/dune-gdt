// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_OPERATORS_WEIGHTED_L2_HH
#define DUNE_GDT_OPERATORS_WEIGHTED_L2_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/operators/integrals.hh>

#include "matrix-based.hh"

namespace Dune {
namespace GDT {


#if 0
// //////////////////////////// //
// WeightedL2LocalizableProduct //
// //////////////////////////// //

template <class WeightFunctionType,
          class GridLayer,
          class Range,
          class Source = Range,
          class Field = typename Range::RangeFieldType>
class WeightedL2LocalizableProduct : public LocalizableProductBase<GridLayer, Range, Source, Field>
{
  typedef LocalizableProductBase<GridLayer, Range, Source, Field> BaseType;

public:
  WeightedL2LocalizableProduct(const WeightFunctionType& weight,
                               GridLayer grd_layr,
                               const Range& rng,
                               const Source& src,
                               const XT::Common::Parameter& param = {})
    : BaseType(grd_layr, rng, src)
    , local_weighted_l2_operator_(weight, param)
  {
    this->append(local_weighted_l2_operator_);
  }

  WeightedL2LocalizableProduct(const size_t over_integrate,
                               const WeightFunctionType& weight,
                               GridLayer grd_layr,
                               const Range& rng,
                               const Source& src,
                               const XT::Common::Parameter& param = {})
    : BaseType(grd_layr, rng, src)
    , local_weighted_l2_operator_(over_integrate, weight, param)
  {
    this->append(local_weighted_l2_operator_);
  }

private:
  const LocalVolumeIntegralOperator<LocalProductIntegrand<WeightFunctionType>,
                                    typename Range::LocalfunctionType,
                                    typename Source::LocalfunctionType,
                                    Field>
      local_weighted_l2_operator_;
}; // class WeightedL2LocalizableProduct


// //////////////////////////////////// //
// make_weighted_l2_localizable_product //
// //////////////////////////////////// //

/**
 * \sa WeightedL2LocalizableProduct
 */
template <class WeightFunctionType, class GridLayerType, class RangeType, class SourceType>
typename std::
    enable_if<XT::Functions::is_localizable_function<WeightFunctionType>::value
                  && XT::Grid::is_layer<GridLayerType>::value
                  && XT::Functions::is_localizable_function<RangeType>::value
                  && XT::Functions::is_localizable_function<SourceType>::value,
              std::unique_ptr<WeightedL2LocalizableProduct<WeightFunctionType, GridLayerType, RangeType, SourceType>>>::
        type
        make_weighted_l2_localizable_product(const WeightFunctionType& weight,
                                             const GridLayerType& grid_layer,
                                             const RangeType& range,
                                             const SourceType& source,
                                             const size_t over_integrate = 0,
                                             const XT::Common::Parameter& param = {})
{
  return Dune::XT::Common::
      make_unique<WeightedL2LocalizableProduct<WeightFunctionType, GridLayerType, RangeType, SourceType>>(
          over_integrate, weight, grid_layer, range, source, param);
}
#endif // 0


// ////////////////////////////// //
// WeightedL2VolumeMatrixOperator //
// ////////////////////////////// //


/**
 * \note See MatrixBasedOperator and OperatorInterface for a description of the template arguments.
 *
 * \note We only provide the most general ctors here and provide make_weighted_l2_matrix_operator for convenience.
 *
 * \sa OperatorInterface
 * \sa MatrixBasedOperator
 * \sa make_weighted_l2_matrix_operator
 */
template <class M,
          class AGV,
          size_t r = 1,
          size_t rC = 1,
          class SF = double,
          class SGV = AGV,
          class F = double,
          class RF = double,
          class RGV = SGV,
          class SV = typename XT::LA::Container<typename M::ScalarType, M::vector_type>::VectorType,
          class RV = SV>
class WeightedL2VolumeMatrixOperator : public MatrixBasedOperator<M, AGV, r, rC, SF, SGV, F, r, rC, RF, RGV, SV, RV>
{
  using ThisType = WeightedL2VolumeMatrixOperator<M, AGV, r, rC, SF, SGV, F, RF, RGV, SV, RV>;
  using BaseType = MatrixBasedOperator<M, AGV, r, rC, SF, SGV, F, r, rC, RF, RGV, SV, RV>;

public:
  using typename BaseType::E;

  using typename BaseType::DofFieldType;
  using typename BaseType::MatrixType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::ElementFilterType;
  using typename BaseType::ApplyOnAllElements;

  using WeightFunctionType = XT::Functions::GridFunctionInterface<E, 1, 1, F>;

private:
  using LocalTwoFormType = LocalElementIntegralOperator<E, r, rC, RF, DofFieldType, r, rC, SF>;
  using LocalIntegrandType = LocalElementProductIntegrand<E, r, rC, RF, DofFieldType, r, rC, SF>;

public:
  /// \name Ctors which accept an existing matrix into which to assemble.
  /// \{

  WeightedL2VolumeMatrixOperator(AGV assembly_grid_view,
                                 const SourceSpaceType& source_spc,
                                 const RangeSpaceType& range_spc,
                                 MatrixType& mat,
                                 const WeightFunctionType& weight_function,
                                 const size_t over_integrate = 0,
                                 const XT::Common::Parameter& param = {},
                                 const ElementFilterType& filter = ApplyOnAllElements())
    : BaseType(assembly_grid_view, source_spc, range_spc, mat)
  {
    this->append(LocalTwoFormType(LocalIntegrandType(weight_function), over_integrate), param, filter);
  }

  /// \}
  /// \name Ctors which create an appropriate matrix into which to assemble.
  /// \{

  WeightedL2VolumeMatrixOperator(AGV assembly_grid_view,
                                 const SourceSpaceType& source_spc,
                                 const RangeSpaceType& range_spc,
                                 const WeightFunctionType& weight_function,
                                 const size_t over_integrate = 0,
                                 const XT::Common::Parameter& param = {},
                                 const ElementFilterType& filter = ApplyOnAllElements())
    : BaseType(assembly_grid_view, source_spc, range_spc, Stencil::element)
  {
    this->append(LocalTwoFormType(LocalIntegrandType(weight_function), over_integrate), param, filter);
  }

  /// \}
}; // class WeightedL2VolumeMatrixOperator


// /////////////////////////////////////// //
// make_weighted_l2_volume_matrix_operator //
// /////////////////////////////////////// //

/// \name Variants of make_weighted_l2_volume_matrix_operator for a given matrix.
/// \{

template <class AGV, class SGV, size_t r, size_t rC, class SF, class RGV, class RF, class M, class F>
WeightedL2VolumeMatrixOperator<typename XT::LA::MatrixInterface<M>::derived_type,
                               GridView<AGV>,
                               r,
                               rC,
                               SF,
                               SGV,
                               F,
                               RF,
                               RGV>
make_weighted_l2_volume_matrix_operator(
    GridView<AGV> assembly_grid_view,
    const SpaceInterface<SGV, r, rC, SF>& source_space,
    const SpaceInterface<RGV, r, rC, RF>& range_space,
    XT::LA::MatrixInterface<M>& matrix,
    const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GridView<AGV>>, 1, 1, F>& weight_function,
    const size_t over_integrate = 0,
    const XT::Common::Parameter& param = {},
    const XT::Grid::ElementFilter<GridView<AGV>>& filter = XT::Grid::ApplyOn::AllElements<GridView<AGV>>())
{
  return WeightedL2VolumeMatrixOperator<typename XT::LA::MatrixInterface<M>::derived_type,
                                        GridView<AGV>,
                                        r,
                                        rC,
                                        SF,
                                        SGV,
                                        F,
                                        RF,
                                        RGV>(
      assembly_grid_view, source_space, range_space, matrix.as_imp(), weight_function, over_integrate, param, filter);
} // ... make_weighted_l2_volume_matrix_operator(...)

template <class AGV, class GV, size_t r, size_t rC, class SF, class M, class F>
WeightedL2VolumeMatrixOperator<typename XT::LA::MatrixInterface<M>::derived_type, GridView<AGV>, r, rC, SF, GV, F>
make_weighted_l2_volume_matrix_operator(
    GridView<AGV> assembly_grid_view,
    const SpaceInterface<GV, r, rC, SF>& space,
    XT::LA::MatrixInterface<M>& matrix,
    const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GridView<AGV>>, 1, 1, F>& weight_function,
    const size_t over_integrate = 0,
    const XT::Common::Parameter& param = {},
    const XT::Grid::ElementFilter<GridView<AGV>>& filter = XT::Grid::ApplyOn::AllElements<GridView<AGV>>())
{
  return WeightedL2VolumeMatrixOperator<typename XT::LA::MatrixInterface<M>::derived_type,
                                        GridView<AGV>,
                                        r,
                                        rC,
                                        SF,
                                        GV,
                                        F>(
      assembly_grid_view, space, space, matrix.as_imp(), weight_function, over_integrate, param, filter);
}

/// \todo Properly implement all the methods below by adding weight_function, over_integrate, param, filter. Do not
/// forget to add these variants in test/operators/operators__weighted_l2_volume_operator.cc.


// template <class GV, size_t r, size_t rC, class F, class M>
// WeightedL2VolumeMatrixOperator<typename XT::LA::MatrixInterface<M>::derived_type, GV, r, rC, F>
// make_weighted_l2_volume_matrix_operator(const SpaceInterface<GV, r, rC, F>& space, XT::LA::MatrixInterface<M>&
// matrix)
//{
//  return WeightedL2VolumeMatrixOperator<typename XT::LA::MatrixInterface<M>::derived_type, GV, r, rC, F>(
//      space, matrix.as_imp());
//}

///// \}
///// \name Variants of make_weighted_l2_volume_matrix_operator, where an appropriate matrix is created .
///// \{

///**
// * \note Use as in
//\code
// auto op = make_weighted_l2_volume_matrix_operator<MatrixType>(source_space, range_space);
//\endcode
// */
// template <class MatrixType,
//          class AGV,
//          class SGV,
//          size_t s_r,
//          size_t s_rC,
//          class SF,
//          class RGV,
//          size_t r_r,
//          size_t r_rC,
//          class RF>
// typename std::enable_if<XT::LA::is_matrix<MatrixType>::value,
//                        WeightedL2VolumeMatrixOperator<MatrixType,
//                                                       GridView<AGV>,
//                                                       s_r,
//                                                       s_rC,
//                                                       SF,
//                                                       SGV,
//                                                       typename XT::Common::multiplication_promotion<SF, RF>::type,
//                                                       r_r,
//                                                       r_rC,
//                                                       RF,
//                                                       RGV>>::type
// make_weighted_l2_volume_matrix_operator(GridView<AGV> assembly_grid_view,
//                                        const SpaceInterface<SGV, s_r, s_rC, SF>& source_space,
//                                        const SpaceInterface<RGV, r_r, r_rC, RF>& range_space)
//{
//  return WeightedL2VolumeMatrixOperator<MatrixType,
//                                        GridView<AGV>,
//                                        s_r,
//                                        s_rC,
//                                        SF,
//                                        SGV,
//                                        typename XT::Common::multiplication_promotion<SF, RF>::type,
//                                        r_r,
//                                        r_rC,
//                                        RF,
//                                        RGV>(assembly_grid_view, source_space, range_space);
//} // ... make_weighted_l2_volume_matrix_operator(...)

///**
// * \note Use as in
//\code
// auto op = make_weighted_l2_volume_matrix_operator<MatrixType>(assembly_grid_view, space);
//\endcode
// */
// template <class MatrixType, class AGV, class GV, size_t r, size_t rC, class F>
// typename std::enable_if<XT::LA::is_matrix<MatrixType>::value,
//                        WeightedL2VolumeMatrixOperator<MatrixType, GridView<AGV>, r, rC, F, GV>>::type
// make_weighted_l2_volume_matrix_operator(GridView<AGV> assembly_grid_view,
//                                        const SpaceInterface<GV, r, rC, F>& space)
//{
//  return WeightedL2VolumeMatrixOperator<MatrixType, GridView<AGV>, r, rC, F, GV>(assembly_grid_view, space);
//}

///**
// * \note Use as in
//\code
// auto op = make_weighted_l2_volume_matrix_operator<MatrixType>(space);
//\endcode
// */
// template <class MatrixType, class GV, size_t r, size_t rC, class F>
// typename std::enable_if<XT::LA::is_matrix<MatrixType>::value,
//                        WeightedL2VolumeMatrixOperator<MatrixType, GV, r, rC, F>>::type
// make_weighted_l2_volume_matrix_operator(const SpaceInterface<GV, r, rC, F>& space)
//{
//  return WeightedL2VolumeMatrixOperator<MatrixType, GV, r, rC, F>(space);
//}

/// \}


#if 0
// ////////////////// //
// WeightedL2Operator //
// ////////////////// //

// forward, needed for the traits
template <class WeightFunctionType, class GridLayer, class Field = typename WeightFunctionType::RangeFieldType>
class WeightedL2Operator;


namespace internal {


template <class WeightFunctionType, class GridLayerType, class Field>
class WeightedL2OperatorTraits
{
public:
  typedef WeightedL2Operator<WeightFunctionType, GridLayerType, Field> derived_type;
  typedef NoJacobian JacobianType;
  typedef Field FieldType;
};


} // namespace internal


template <class WeightFunctionType, class GridLayerType, class Field>
class WeightedL2Operator
    : public OperatorInterface<internal::WeightedL2OperatorTraits<WeightFunctionType, GridLayerType, Field>>
{
  typedef OperatorInterface<internal::WeightedL2OperatorTraits<WeightFunctionType, GridLayerType, Field>> BaseType;

public:
  using typename BaseType::FieldType;

  WeightedL2Operator(const WeightFunctionType& weight, GridLayerType grid_layer, const size_t over_integrate = 0)
    : weight_(weight)
    , grid_layer_(grid_layer)
    , over_integrate_(over_integrate)
  {
  }

  template <class SourceSpaceType, class VectorType, class RangeSpaceType>
  void apply(const DiscreteFunction<SourceSpaceType, VectorType>& source,
             DiscreteFunction<RangeSpaceType, VectorType>& range,
             const XT::Common::Parameter& param = {}) const
  {
    typedef
        typename XT::LA::Container<typename VectorType::ScalarType, VectorType::Traits::sparse_matrix_type>::MatrixType
            MatrixType;
    auto op = make_weighted_l2_matrix_operator<MatrixType>(
        weight_, source.space(), range.space(), grid_layer_, over_integrate_);
    op->apply(source, range, param);
  }

  template <class E, class D, size_t d, class R, size_t r, size_t rC>
  FieldType apply2(const XT::Functions::GridFunctionInterface<E, D, d, R, r, rC>& range,
                   const XT::Functions::GridFunctionInterface<E, D, d, R, r, rC>& source,
                   const XT::Common::Parameter& param = {}) const
  {
    auto product = make_weighted_l2_localizable_product(weight_, grid_layer_, range, source, over_integrate_, param);
    return product->apply2();
  }

  using BaseType::apply_inverse;

  template <class RangeType, class SourceType>
  void
  apply_inverse(const RangeType& /*range*/, SourceType& /*source*/, const XT::Common::Configuration& /*opts*/) const
  {
    DUNE_THROW(NotImplemented, "yet");
  }

  std::vector<std::string> invert_options() const
  {
    DUNE_THROW(NotImplemented, "yet");
    return {"depends_on_the_vector_type_of_the_discrete_function"};
  }

  XT::Common::Configuration invert_options(const std::string& /*type*/) const
  {
    DUNE_THROW(NotImplemented, "yet");
  }

private:
  const WeightFunctionType& weight_;
  GridLayerType grid_layer_;
  const size_t over_integrate_;
}; // class WeightedL2Operator


// ///////////////////////// //
// make_weighted_l2_operator //
// ///////////////////////// //

template <class GridLayerType, class WeightFunctionType>
typename std::enable_if<XT::Functions::is_localizable_function<WeightFunctionType>::value
                            && XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<WeightedL2Operator<WeightFunctionType,
                                                           GridLayerType,
                                                           typename WeightFunctionType::RangeFieldType>>>::type
make_weighted_l2_operator(const GridLayerType& grid_layer,
                          const WeightFunctionType& weight,
                          const size_t over_integrate = 0)
{
  return Dune::XT::Common::
      make_unique<WeightedL2Operator<WeightFunctionType, GridLayerType, typename WeightFunctionType::RangeFieldType>>(
          weight, grid_layer, over_integrate);
}
#endif // 0


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_WEIGHTED_L2_HH
