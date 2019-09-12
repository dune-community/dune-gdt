// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_TOOLS_local_mass_matrices_HH
#define DUNE_GDT_TOOLS_local_mass_matrices_HH

#include <mutex>
#include <vector>

#include <dune/xt/common/parallel/threadstorage.hh>
#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/conversion.hh>
#include <dune/xt/la/matrix-inverter.hh>
#include <dune/xt/grid/functors/interfaces.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/mapper/finite-volume.hh>

namespace Dune {
namespace GDT {


template <class GV, size_t r = 1, size_t rC = 1, class F = double, class AGV = GV>
class LocalMassMatrixProvider
  : public XT::Grid::ElementFunctor<GV>
  , public XT::Common::ThreadResultPropagator<
        LocalMassMatrixProvider<GV, r, rC, F>,
        std::map<size_t, std::pair<XT::LA::CommonDenseMatrix<F>, XT::LA::CommonDenseMatrix<F>>>,
        XT::Common::concatenate_container<
            std::map<size_t, std::pair<XT::LA::CommonDenseMatrix<F>, XT::LA::CommonDenseMatrix<F>>>>>
{
  static_assert(XT::Grid::is_view<AGV>::value, "");

  using ThisType = LocalMassMatrixProvider<GV, r, rC, F>;
  using BaseType = XT::Grid::ElementFunctor<GV>;
  using Propagator = XT::Common::ThreadResultPropagator<
      LocalMassMatrixProvider<GV, r, rC, F>,
      std::map<size_t, std::pair<XT::LA::CommonDenseMatrix<F>, XT::LA::CommonDenseMatrix<F>>>,
      XT::Common::concatenate_container<
          std::map<size_t, std::pair<XT::LA::CommonDenseMatrix<F>, XT::LA::CommonDenseMatrix<F>>>>>;
  friend Propagator;

public:
  using AssemblyGridView = AGV;
  using SpaceType = SpaceInterface<GV, r, rC, F>;
  using typename BaseType::E;
  using typename BaseType::ElementType;

  LocalMassMatrixProvider(const AssemblyGridView& grid_view, const SpaceType& space)
    : BaseType()
    , Propagator(this)
    , grid_view_(grid_view)
    , space_(space)
    , element_mapper_(grid_view_)
    , instance_counter_(0)
    , local_basis_(space_.basis().localize())
  {}

  LocalMassMatrixProvider(const ThisType& other)
    : BaseType(other)
    , Propagator(this)
    , grid_view_(other.grid_view_)
    , space_(other.space_)
    , element_mapper_(grid_view_)
    , instance_counter_(other.instance_counter_ + 1)
    , local_basis_(space_.basis().localize())
  {}

  void apply_local(const ElementType& element) override
  {
    local_basis_->bind(element);
    const size_t id = element_mapper_.global_index(element, 0);
    const LocalElementIntegralBilinearForm<E, r, rC, F, F> local_l2_bilinear_form(
        LocalElementProductIntegrand<E, r, F, F>(1.));
    auto matrix =
        XT::LA::convert_to<XT::LA::CommonDenseMatrix<F>>(local_l2_bilinear_form.apply2(*local_basis_, *local_basis_));
    auto inverse = XT::LA::invert_matrix(matrix);
    local_mass_matrices_.insert(std::make_pair(id, std::make_pair(std::move(matrix), std::move(inverse))));
  } // ... apply_local(...)

  BaseType* copy() override final
  {
    return Propagator::copy_imp();
  }

  void finalize() override final
  {
    this->finalize_imp();
  }

  /**
   * \note This is only used to merge the thread-local results after the grid walk, you are probably interested in
   *       local_mass_matrix() and local_mass_matrix_inverse().
   */
  std::map<size_t, std::pair<XT::LA::CommonDenseMatrix<F>, XT::LA::CommonDenseMatrix<F>>> result() const
  {
    return local_mass_matrices_;
  }

  void set_result(std::map<size_t, std::pair<XT::LA::CommonDenseMatrix<F>, XT::LA::CommonDenseMatrix<F>>> res)
  {
    local_mass_matrices_ = res;
  }

  const XT::LA::CommonDenseMatrix<F>& local_mass_matrix(const ElementType& element) const
  {
    const size_t id = element_mapper_.global_index(element, 0);
    DUNE_THROW_IF(local_mass_matrices_.count(id) == 0,
                  XT::Common::Exceptions::this_should_not_happen,
                  "Missing local mass matrix for id " << id << "!");
    return local_mass_matrices_.at(id).first;
  }

  const XT::LA::CommonDenseMatrix<F>& local_mass_matrix_inverse(const ElementType& element) const
  {
    const size_t id = element_mapper_.global_index(element, 0);
    DUNE_THROW_IF(local_mass_matrices_.count(id) == 0,
                  XT::Common::Exceptions::this_should_not_happen,
                  "Missing local mass matrix for id " << id << "!");
    return local_mass_matrices_.at(id).second;
  }

private:
  const AssemblyGridView grid_view_;
  const SpaceType& space_;
  const FiniteVolumeMapper<GV> element_mapper_;
  size_t instance_counter_;
  std::unique_ptr<typename SpaceType::GlobalBasisType::LocalizedType> local_basis_;
  std::map<size_t, std::pair<XT::LA::CommonDenseMatrix<F>, XT::LA::CommonDenseMatrix<F>>> local_mass_matrices_;
}; // class LocalMassMatrixProvider


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TOOLS_local_mass_matrices_HH
