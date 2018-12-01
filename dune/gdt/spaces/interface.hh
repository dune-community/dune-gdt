// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014, 2016 - 2018)
//   Sven Kaulmann   (2014)
//   Tobias Leibner  (2014, 2016 - 2017)

#ifndef DUNE_GDT_SPACES_INTERFACE_HH
#define DUNE_GDT_SPACES_INTERFACE_HH

#include <dune/geometry/type.hh>

#include <dune/grid/utility/persistentcontainer.hh>

#include <dune/xt/la/container/vector-interface.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/finite-elements/interfaces.hh>
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
class SpaceInterface
{
  static_assert(XT::Grid::is_view<GridView>::value, "");

public:
  using GridViewType = GridView;
  using GV = GridViewType;
  using ElementType = XT::Grid::extract_entity_t<GV>;
  using E = ElementType;
  using Grid = typename GridView::Grid;
  using G = Grid;
  using D = typename GridViewType::ctype;
  static const constexpr size_t d = GridViewType::dimension;
  using R = RangeField;
  static const constexpr size_t r = range_dim;
  static const constexpr size_t rC = range_dim_columns;

  using GlobalBasisType = GlobalBasisInterface<GridViewType, r, rC, R>;
  using MapperType = MapperInterface<GridViewType>;
  using FiniteElementType = LocalFiniteElementInterface<D, d, R, r, rC>;

  using DofCommunicatorType = typename DofCommunicationChooser<GridViewType>::Type;

  SpaceInterface()
    : dof_communicator_(nullptr)
    , pre_adapted_(false)
    , adapted_(false)
  {}

  virtual ~SpaceInterface() = default;

  /// \name These methods provide most functionality, they have to be implemented.
  /// \{

  virtual const GridViewType& grid_view() const = 0;

  virtual const MapperType& mapper() const = 0;

  virtual const GlobalBasisType& basis() const = 0;

  virtual const FiniteElementType& finite_element(const GeometryType& geometry_type) const = 0;

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
   * If this returns true, every finite_element() is expected to provide lagrange_points().
   */
  virtual bool is_lagrangian() const = 0;

  /// \}
  /// \name These methods are required for grid adaptation.
  /// \{

  virtual void restrict_to(const ElementType& /*element*/,
                           PersistentContainer<G, DynamicVector<R>>& /*persistent_data*/) const
  {
    DUNE_THROW(Exceptions::space_error, "This space does not support adaptation!");
  }

  virtual void pre_adapt()
  {
    if (pre_adapted_)
      return;
    pre_adapted_ = true;
  }

protected:
  /**
   * \todo Detect if this is called more than once per adaptation.
   */
  virtual void update_after_adapt()
  {
    DUNE_THROW(Exceptions::space_error, "This space does not support adaptation!");
  }

public:
  virtual void adapt()
  {
    DUNE_THROW_IF(!pre_adapted_, Exceptions::space_error, "You need to call pre_adapt() first!");
    if (adapted_)
      return;
    this->update_after_adapt();
  }

  virtual void post_adapt()
  {
    adapted_ = false;
    pre_adapted_ = false;
  }

  virtual void prolong_onto(const ElementType& /*element*/,
                            const PersistentContainer<G, DynamicVector<R>>& /*persistent_data*/,
                            DynamicVector<R>& /*element_data*/) const
  {
    DUNE_THROW(Exceptions::space_error, "This space does not support adaptation!");
  }

  virtual DynamicVector<R> prolong_onto(const ElementType& element,
                                        const PersistentContainer<G, DynamicVector<R>>& persistent_data) const
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
  bool pre_adapted_;
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
