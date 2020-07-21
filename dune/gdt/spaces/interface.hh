// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2018)
//   René Fritze     (2014, 2016 - 2018)
//   René Milk       (2017)
//   Sven Kaulmann   (2014)
//   Tobias Leibner  (2014, 2016 - 2017)

#ifndef DUNE_GDT_SPACES_INTERFACE_HH
#define DUNE_GDT_SPACES_INTERFACE_HH

#include <dune/geometry/type.hh>

#include <dune/grid/utility/persistentcontainer.hh>

#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/la/container/vector-interface.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/generic/element-function.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/finite-elements/interfaces.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/integrands/conversion.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/print.hh>
#include <dune/gdt/spaces/basis/interface.hh>
#include <dune/gdt/spaces/mapper/interfaces.hh>
#include <dune/gdt/spaces/parallel/communication.hh>
#include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {


// forward
template <class V, class GV, size_t r, size_t rC, class R>
class ConstDiscreteFunction;


template <class GridView, size_t range_dim = 1, size_t range_dim_columns = 1, class RangeField = double>
class SpaceInterface : public XT::Common::WithLogger<SpaceInterface<GridView, range_dim, range_dim_columns, RangeField>>
{
  static_assert(XT::Grid::is_view<GridView>::value, "");
  using ThisType = SpaceInterface;
  using Logger = XT::Common::WithLogger<SpaceInterface<GridView, range_dim, range_dim_columns, RangeField>>;

public:
  using GridViewType = GridView;
  using GV = GridViewType;
  using ElementType = XT::Grid::extract_entity_t<GV>;
  using E = ElementType;
  using GridType = typename GridView::Grid;
  using G = GridType;
  using D = typename GridViewType::ctype;
  static const constexpr size_t d = GridViewType::dimension;
  using R = RangeField;
  static const constexpr size_t r = range_dim;
  static const constexpr size_t rC = range_dim_columns;

  using GlobalBasisType = GlobalBasisInterface<GridViewType, r, rC, R>;
  using LocalFiniteElementFamilyType = LocalFiniteElementFamilyInterface<D, d, R, r, rC>;
  using MapperType = MapperInterface<GridViewType>;

  using DofCommunicatorType = typename DofCommunicationChooser<GridViewType>::Type;

  SpaceInterface(const std::string& logging_prefix = "",
                 const std::string& logging_id_ = "",
                 const bool logging_disabled = true)
    : Logger(
        logging_prefix.empty() ? "gdt" : logging_prefix, logging_id_.empty() ? "Space" : logging_id_, logging_disabled)
    , dof_communicator_(nullptr)
    , adapted_(false)
  {}

  SpaceInterface(const ThisType& other) = default;

  SpaceInterface(ThisType&& source) = default;

  virtual ThisType* copy() const = 0;

  virtual ~SpaceInterface() = default;

  /// \name These methods provide most functionality, they have to be implemented.
  /// \{

  virtual const GridViewType& grid_view() const = 0;

  virtual const MapperType& mapper() const = 0;

  virtual const LocalFiniteElementFamilyType& finite_elements() const = 0;

  virtual const GlobalBasisType& basis() const = 0;

  /// \}
  /// \name These methods help to identify the space, they have to be implemented.
  /// \{

  virtual SpaceType type() const = 0;

  virtual int min_polorder() const = 0;

  virtual int max_polorder() const = 0;

  /**
   * To query if elements of this space (functions) are continuous (i.e., in C^0), use `continuous(0)`.
   */
  virtual bool continuous(const int diff_order) const = 0;

  virtual bool continuous_normal_components() const = 0;

  /**
   * If this returns true, every FE provided by local_finite_elements() is expected to provide lagrange_points().
   */
  virtual bool is_lagrangian() const = 0;

  /// \}
  /// \name These methods are required for grid adaptation.
  /// \{

  /**
   * \brief Computes an L^2 projection of the data defined on the finer elements, on the domain covered by the
   *        coarser element.
   *
   * \attention This implementation only makes sense if the evaluation of the global basis coincides with the evaluation
   *            of the finite elements basis.
   *
   * \note Override this method if this is not the correct choice for the space in question!
   */
  virtual void restrict_to(const ElementType& element,
                           PersistentContainer<G, std::pair<DynamicVector<R>, DynamicVector<R>>>& persistent_data) const
  {
    DUNE_THROW_IF(rC != 1,
                  Exceptions::space_error,
                  "Not implemented for matrix-valued spaces yet (due to LocalProductIntegrand)!");
    auto& element_restriction_data = persistent_data[element];
    auto& element_restriction_FE_data = element_restriction_data.first;
    auto& element_restriction_DoF_data = element_restriction_data.second;
    if (element_restriction_DoF_data.size() == 0) {
      DUNE_THROW_IF(
          element.isLeaf(), Exceptions::space_error, "Restriction to a leaf element does not make any sense!");
      // We do not have meaningful data to initialize a global FE, since element is in general not contained in the
      // spaces grid view, and the basis need only make sense on those. We thus use some default data, which need not
      // make sense as a global basis. This does not matter, however, since we only need to be able to restore that
      // same basis later on to make sense of the DoF information.
      auto element_basis = this->basis().localize();
      element_restriction_FE_data = element_basis->default_data(element.type());
      element_basis->restore(element, element_restriction_FE_data);
      auto lhs = LocalElementIntegralBilinearForm<E, r, rC, R, R>(LocalProductIntegrand<E, r, R, R>())
                     .apply2(*element_basis, *element_basis);
      DynamicVector<R> rhs(element_basis->size(), 0.);
      for (auto&& child_element : descendantElements(element, element.level() + 1)) {
        // ensure we have data on all descendant elements of the next level
        this->restrict_to(child_element, persistent_data);
        // prepare
        const auto& child_restriction_data = persistent_data[child_element];
        const auto& child_restriction_FE_data = child_restriction_data.first;
        auto child_basis = this->basis().localize();
        child_basis->restore(child_element, child_restriction_FE_data);
        XT::Functions::GenericElementFunctionSet<E, r, rC, R> element_basis_on_child_element(
            element_basis->size(),
            element_basis->order(),
            [&](const auto& x_in_child_reference_element_coordinates, auto& result, const auto&) {
              const auto x_in_physical_coordinates =
                  child_element.geometry().global(x_in_child_reference_element_coordinates);
              const auto x_in_reference_element_coordinates = element.geometry().local(x_in_physical_coordinates);
              element_basis->evaluate(x_in_reference_element_coordinates, result);
            });
        element_basis_on_child_element.bind(child_element);
        std::vector<typename GlobalBasisType::LocalizedType::RangeType> child_basis_values(child_basis->size());
        const auto& child_restriction_DoF_data = child_restriction_data.second;
        XT::Functions::GenericElementFunction<E, r, rC, R> restriction_data_as_function_on_child_element(
            child_basis->order(), [&](const auto& x, const auto&) {
              child_basis->evaluate(x, child_basis_values);
              std::remove_reference_t<decltype(child_basis_values[0])> result(0.);
              for (size_t ii = 0; ii < child_basis->size(); ++ii)
                result += child_basis_values[ii] * child_restriction_DoF_data[ii];
              return result;
            });
        // compute rhs for the L2 projection
        rhs += LocalElementIntegralFunctional<E, r, rC, R, R>(
                   /*order=*/
                   [&](const auto& basis, const auto&) {
                     return std::max(basis.order(), restriction_data_as_function_on_child_element.order());
                   },
                   /*integrand_evaluation=*/
                   [&](const auto& basis, const auto& x, auto& result, const auto&) {
                     const auto restriction_data_value = restriction_data_as_function_on_child_element.evaluate(x);
                     basis.evaluate(x, child_basis_values);
                     for (size_t ii = 0; ii < basis.size(); ++ii)
                       result = child_basis_values[ii] * restriction_data_value;
                   })
                   .apply(element_basis_on_child_element);
      }
      auto restricted_DoFs = rhs;
      XT::LA::make_solver(lhs).apply(rhs, restricted_DoFs);
      element_restriction_DoF_data.resize(restricted_DoFs.size());
      for (size_t ii = 0; ii < restricted_DoFs.size(); ++ii)
        element_restriction_DoF_data[ii] = restricted_DoFs[ii];
    }
  } // ... restrict_to(...)

  void pre_adapt() {}

protected:
  /// \note Ensure in derived classes that the basis, mapper and communicator are updated!
  virtual void update_after_adapt() {}

public:
  virtual void adapt()
  {
    if (adapted_)
      return;
    this->update_after_adapt();
    adapted_ = true;
  }

  void post_adapt()
  {
    adapted_ = false;
  }

  /**
   * \brief Computes a prolongation using the interpolation of the finite_element().
   *
   * \attention This implementation only makes sense if the evaluation of the global basis coincides with the evaluation
   *            of the finite elements basis.
   *
   * \note Override this method if this is not the correct choice for the space in question.
   */
  virtual void
  prolong_onto(const ElementType& element,
               const PersistentContainer<G, std::pair<DynamicVector<R>, DynamicVector<R>>>& persistent_data,
               DynamicVector<R>& element_dof_data) const
  {
    if (element.isNew()) {
      // This is a new element where we do not have any data, so we interpolate from a father element.
      // - first find a father element with data (does not have to be the direct father, since a conforming grid may
      //   create several levels of refinement) ...
      auto father = std::make_unique<ElementType>(element.father());
      auto father_data = persistent_data[*father];
      //   ... which we do by looking a the DoF data (must not be empty)
      while (father_data.second.size() == 0) {
        DUNE_THROW_IF(father->isLeaf(), Exceptions::space_error, "Cound not find a suitable father with data!");
        father = std::make_unique<ElementType>(father->father());
        father_data = persistent_data[*father];
      }
      const auto& father_fe_data = father_data.first;
      const auto& father_dof_data = father_data.second;
      // - restore the global basis from the fe data
      auto father_basis = this->basis().localize();
      father_basis->restore(*father, father_fe_data);
      // - get the basis for the current element (has to be available after update_after_adapt())
      const auto element_basis = this->basis().localize(element);
      // - interpolate the data from the father to the element
      std::vector<typename GlobalBasisType::LocalizedType::RangeType> father_basis_values(father_basis->size());
      DUNE_THROW_IF(father_dof_data.size() != father_basis->size(),
                    Exceptions::space_error,
                    "element: " << print(element) << "\nelement.level() = " << element.level()
                                << "\nfather: " << print(*father) << "\nfather.level() = " << father->level()
                                << "\nfather_dof_data.size() = " << father_dof_data.size()
                                << "\nfather_basis->size() = " << father_basis->size());
      element_basis->interpolate(
          [&](const auto& point_in_element_reference_element_coordinates) {
            const auto point_in_physical_coordinates =
                element.geometry().global(point_in_element_reference_element_coordinates);
            const auto point_in_father_reference_element_coordinates =
                father->geometry().local(point_in_physical_coordinates);
            father_basis->evaluate(point_in_father_reference_element_coordinates, father_basis_values);
            std::remove_reference_t<decltype(father_basis_values[0])> result(0.);
            for (size_t ii = 0; ii < father_basis->size(); ++ii)
              result += father_basis_values[ii] * father_dof_data[ii];
            return result;
          },
          father_basis->order(),
          element_dof_data);
    } else {
      // This grid element is unchanged, but the FE may have changed (if the FE depends on global information), so we
      // cannot simply copy the old data, but rather interpolate it anew.
      // - restore the original FE
      const auto& original_element_data = persistent_data[element];
      const auto& original_element_FE_data = original_element_data.first;
      const auto& original_element_DoF_data = original_element_data.second;
      auto original_basis = this->basis().localize();
      original_basis->restore(element, original_element_FE_data);
      // - get the basis for the current element (has to be available after update_after_adapt())
      const auto new_basis = this->basis().localize(element);
      std::vector<typename GlobalBasisType::LocalizedType::RangeType> original_basis_values(original_basis->size());
      DUNE_THROW_IF(original_element_DoF_data.size() != original_basis->size(),
                    Exceptions::space_error,
                    "element: " << print(element) << "\nelement.level() = " << element.level()
                                << "\noriginal_element_DoF_data.size() = " << original_element_DoF_data.size()
                                << "\noriginal_basis->size() = " << original_basis->size());
      // - interpolate the data from the father to the element (no need to map the coordinate, same geometry)
      new_basis->interpolate(
          [&](const auto& xx) {
            original_basis->evaluate(xx, original_basis_values);
            std::remove_reference_t<decltype(original_basis_values[0])> result(0.);
            for (size_t ii = 0; ii < original_basis->size(); ++ii)
              result += original_basis_values[ii] * original_element_DoF_data[ii];
            return result;
          },
          original_basis->order(),
          element_dof_data);
    }
  } // ... prolong_onto(...)

  virtual DynamicVector<R>
  prolong_onto(const ElementType& element,
               const PersistentContainer<G, std::pair<DynamicVector<R>, DynamicVector<R>>>& persistent_data) const
  {
    DynamicVector<R> element_data;
    this->prolong_onto(element, persistent_data, element_data);
    return element_data;
  }

  /// \}
  /// \name These methods are required for MPI communication, they are provided.
  /// \{

  virtual const DofCommunicatorType& dof_communicator() const
  {
    if (!dof_communicator_)
      DUNE_THROW(Exceptions::space_error,
                 "The actual space has to either implement its own dof_communicator() or call "
                 "create_communicator() in the ctor!");
    return *dof_communicator_;
  }

  /// \}
  /// \name These methods are provided for convenience.
  /// \{

  template <class V>
  bool contains(const XT::LA::VectorInterface<V>& vector) const
  {
    return vector.size() == this->mapper().size();
  }

  /**
   * \note A return value of true cannot be ultimately trusted yet.
   *
   * \sa https://github.com/dune-community/dune-gdt/issues/123
   */
  template <class V>
  bool contains(const ConstDiscreteFunction<V, GV, r, rC, R>& function) const
  {
    // this is the only case where we are sure^^
    if (&function.space() == this)
      return true;
    if (function.space().type() != this->type())
      return false;
    if (function.space().mapper().size() != this->mapper().size())
      return false;
    // the spaces might still differ (different grid views of same size), but we have no means to check this
    return true;
  }

  /// \}

protected:
  void create_communicator()
  {
    if (!dof_communicator_) {
      dof_communicator_ =
          std::shared_ptr<DofCommunicatorType>(DofCommunicationChooser<GridViewType>::create(grid_view()));
      DofCommunicationChooser<GridViewType>::prepare(*this, *dof_communicator_);
    }
  }

private:
  std::shared_ptr<DofCommunicatorType> dof_communicator_;
  bool adapted_;
}; // class SpaceInterface


template <class GV, size_t r, size_t rC, class R>
std::ostream& operator<<(std::ostream& out, const SpaceInterface<GV, r, rC, R>& space)
{
  out << "Space(" << space.type() << ", " << space.mapper().size() << " DoFs)";
  return out;
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_INTERFACE_HH
