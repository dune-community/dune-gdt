// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_INTERFACE_HH
#define DUNE_GDT_SPACES_INTERFACE_HH

#include <memory>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/gridview.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/common/parallel/threadstorage.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/common/ranges.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/grid/layers.hh>
#include <dune/stuff/la/container/pattern.hh>

#include "constraints.hh"

namespace Dune {
namespace GDT {


enum class ChooseSpaceBackend
{
    gdt
  , pdelab
  , fem
}; // enum class ChooseSpaceBackend


static constexpr ChooseSpaceBackend default_space_backend =
#if HAVE_DUNE_FEM
                                                            ChooseSpaceBackend::fem;
#elif HAVE_DUNE_PDELAB
                                                            ChooseSpaceBackend::pdelab;
#else
                                                            ChooseSpaceBackend::gdt;
#endif


template< ChooseSpaceBackend type >
struct ChooseGridPartView;


template<>
struct ChooseGridPartView< ChooseSpaceBackend::gdt >
{
  static const Stuff::Grid::ChoosePartView type = Stuff::Grid::ChoosePartView::view;
};


template<>
struct ChooseGridPartView< ChooseSpaceBackend::pdelab >
{
  static const Stuff::Grid::ChoosePartView type = Stuff::Grid::ChoosePartView::view;
};


template<>
struct ChooseGridPartView< ChooseSpaceBackend::fem >
{
  static const Stuff::Grid::ChoosePartView type = Stuff::Grid::ChoosePartView::part;
};


template< class Traits, size_t domainDim, size_t rangeDim, size_t rangeDimCols = 1 >
class SpaceInterface
  : public Stuff::CRTPInterface< SpaceInterface< Traits, domainDim, rangeDim, rangeDimCols >, Traits >
{
public:
  typedef typename Traits::derived_type         derived_type;
  static const int                              polOrder = Traits::polOrder;
  typedef typename Traits::BackendType          BackendType;
  typedef typename Traits::MapperType           MapperType;
  typedef typename Traits::BaseFunctionSetType  BaseFunctionSetType;
  typedef typename Traits::CommunicatorType     CommunicatorType;
  typedef typename Traits::GridViewType         GridViewType;
  typedef typename Traits::RangeFieldType       RangeFieldType;
  static const size_t                           dimDomain = domainDim;
  static const size_t                           dimRange = rangeDim;
  static const size_t                           dimRangeCols = rangeDimCols;

private:
  static_assert(dimDomain > 0, "dimDomain has to be positive");
  static_assert(dimRange > 0, "dimRange has to be positive");
  static_assert(dimRangeCols > 0, "dimRangeCols has to be positive");
  static_assert(std::is_base_of< GridView< typename GridViewType::Traits >, GridViewType >::value,
                "GridViewType has to be derived from GridView!");
  static_assert(GridViewType::dimension == dimDomain, "Dimension of GridView has to match dimDomain");

public:
  typedef typename GridViewType::ctype              DomainFieldType;
  typedef FieldVector< DomainFieldType, dimDomain > DomainType;

  typedef typename GridViewType::template Codim< 0 >::Entity     EntityType;
  typedef typename GridViewType::Intersection                    IntersectionType;
  typedef Stuff::Grid::BoundaryInfoInterface< IntersectionType > BoundaryInfoType;
  typedef Dune::Stuff::LA::SparsityPatternDefault                PatternType;

  static const Stuff::Grid::ChoosePartView part_view_type = Traits::part_view_type;

  static const bool needs_grid_view = Traits::needs_grid_view;

public:
  /**
   * \defgroup interface ´´These methods have to be implemented!''
   * @{
   **/

  const GridViewType& grid_view() const
  {
    CHECK_CRTP(this->as_imp().grid_view());
    return this->as_imp().grid_view();
  }

  const BackendType& backend() const
  {
    CHECK_CRTP(this->as_imp().backend());
    return this->as_imp().backend();
  }

  const MapperType& mapper() const
  {
    CHECK_CRTP(this->as_imp().mapper());
    return this->as_imp().mapper();
  }

  BaseFunctionSetType base_function_set(const EntityType& entity) const
  {
    CHECK_CRTP(this->as_imp().base_function_set(entity));
    return this->as_imp().base_function_set(entity);
  }

  CommunicatorType& communicator() const
  {
    CHECK_CRTP(this->as_imp().communicator());
    return this->as_imp().communicator();
  }

  /**
   *  \brief Computes local constraints.
   *
   *  \note  Any derived class has to implement this method, even if it does not support any kind of constraints!
   *         In that case just provide exactly the following method:\code
template< class S, size_t d, size_t r, size_t rC, class ConstraintsType >
void local_constraints(const SpaceInterface< S, d, r, rC > >&, const EntityType&, ConstraintsType&) const
{
  static_assert(AlwaysFalse< S >::value, "Not implemented for these constraints!");
}
\endcode
   */
  template< class S, size_t d, size_t r, size_t rC, class C >
  void local_constraints(const SpaceInterface< S, d, r, rC >& ansatz_space,
                         const EntityType& entity,
                         Spaces::ConstraintsInterface< C >& ret) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().local_constraints(ansatz_space.as_imp(), entity, ret.as_imp()));
  }

  template< class C >
  void local_constraints(const EntityType& entity, Spaces::ConstraintsInterface< C >& ret) const
  {
    local_constraints(*this, entity, ret);
  }

  /**
   *  \brief  Computes the appropriate sparsity pattern.
   *  \note   This method can be implemented in a derived class by a forward to one of the methods provided by this class, namely compute_volume_pattern(), compute_face_pattern() or compute_face_and_volume_pattern().
   */
  template< class G, class S, size_t d, size_t r, size_t rC >
  PatternType compute_pattern(const GridView< G >& local_grid_view, const SpaceInterface< S, d, r, rC >& ansatz_space) const
  {
    CHECK_CRTP(this->as_imp().compute_pattern(local_grid_view, ansatz_space.as_imp()));
    return this->as_imp().compute_pattern(local_grid_view, ansatz_space.as_imp());
  }
  /* @} */

  /**
   * \defgroup provided ´´These methods are provided by the interface for convenience.''
   * @{
   **/

  PatternType compute_pattern() const
  {
    return compute_pattern(*this);
  }

  template< class S, size_t d, size_t r, size_t rC >
  PatternType compute_pattern(const SpaceInterface< S, d, r, rC >& ansatz_space) const
  {
    return compute_pattern(grid_view(), ansatz_space);
  }

  template< class G >
  PatternType compute_pattern(const GridView< G >& local_grid_view) const
  {
    return compute_pattern(local_grid_view, *this);
  }

  PatternType compute_volume_pattern() const
  {
    return compute_volume_pattern(*this);
  }

  template< class S, size_t d, size_t r, size_t rC >
  PatternType compute_volume_pattern(const SpaceInterface< S, d, r, rC >& ansatz_space) const
  {
    return compute_volume_pattern(grid_view(), ansatz_space);
  }

  template< class G >
  PatternType compute_volume_pattern(const GridView< G >& local_grid_view) const
  {
    return compute_volume_pattern(local_grid_view, *this);
  }

  /**
   *  \brief  computes a sparsity pattern, where this space is the test space (rows/outer) and the other space is the
   *          ansatz space (cols/inner)
   */
  template< class G, class S, size_t d, size_t r, size_t rC >
  PatternType compute_volume_pattern(const GridView< G >& local_grid_view,
                                     const SpaceInterface< S, d, r, rC >& ansatz_space) const
  {
    PatternType pattern(mapper().size());
    Dune::DynamicVector< size_t > globalRows(mapper().maxNumDofs(), 0);
    Dune::DynamicVector< size_t > globalCols(ansatz_space.mapper().maxNumDofs(), 0);

    for (const auto& entity : DSC::entityRange(local_grid_view)) {
      const auto testBase = base_function_set(entity);
      const auto ansatzBase = ansatz_space.base_function_set(entity);
      mapper().globalIndices(entity, globalRows);
      ansatz_space.mapper().globalIndices(entity, globalCols);
      for (size_t ii = 0; ii < testBase.size(); ++ii) {
        for (size_t jj = 0; jj < ansatzBase.size(); ++jj) {
          pattern.insert(globalRows[ii], globalCols[jj]);
        }
      }
    }
    pattern.sort();
    return pattern;
  } // ... compute_volume_pattern(...)

  PatternType compute_face_and_volume_pattern() const
  {
    return compute_face_and_volume_pattern(grid_view(), *this);
  }

  template< class G >
  PatternType compute_face_and_volume_pattern(const GridView< G >& local_grid_view) const
  {
    return compute_face_and_volume_pattern(local_grid_view, *this);
  }

  template< class S, size_t d, size_t r, size_t rC >
  PatternType compute_face_and_volume_pattern(const SpaceInterface< S, d, r, rC >& ansatz_space) const
  {
    return compute_face_and_volume_pattern(grid_view(), ansatz_space);
  }

  /**
   *  \brief  computes a DG sparsity pattern, where this space is the test space (rows/outer) and the other space is the
   *          ansatz space (cols/inner)
   */
  template< class G, class S, size_t d, size_t r, size_t rC >
  PatternType compute_face_and_volume_pattern(const GridView< G >& local_grid_view,
                                              const SpaceInterface< S, d, r, rC >& ansatz_space) const
  {
    // prepare
    PatternType pattern(mapper().size());
    Dune::DynamicVector< size_t > global_rows(mapper().maxNumDofs(), 0);
    Dune::DynamicVector< size_t > global_cols(ansatz_space.mapper().maxNumDofs(), 0);
    for (const auto& entity : DSC::entityRange(local_grid_view)) {
      const auto test_base_entity = base_function_set(entity);
      const auto ansatz_base_entity = ansatz_space.base_function_set(entity);
      mapper().globalIndices(entity, global_rows);
      ansatz_space.mapper().globalIndices(entity, global_cols);
      // compute entity/entity
      for (size_t ii = 0; ii < test_base_entity.size(); ++ii) {
        for (size_t jj = 0; jj < ansatz_base_entity.size(); ++jj) {
          pattern.insert(global_rows[ii], global_cols[jj]);
        }
      }
      // walk the intersections
      const auto intersection_it_end = local_grid_view.iend(entity);
      for (auto intersection_it = local_grid_view.ibegin(entity);
           intersection_it != intersection_it_end;
           ++intersection_it) {
        const auto& intersection = *intersection_it;
        // get the neighbour
        if (intersection.neighbor() && !intersection.boundary()) {
          const auto neighbour_ptr = intersection.outside();
          const auto& neighbour = *neighbour_ptr;
          // get the basis
          const auto ansatz_base_neighbour = ansatz_space.base_function_set(neighbour);
          ansatz_space.mapper().globalIndices(neighbour, global_cols);
          // compute entity/neighbour
          for (size_t ii = 0; ii < test_base_entity.size(); ++ii) {
            for (size_t jj = 0; jj < ansatz_base_neighbour.size(); ++jj) {
              pattern.insert(global_rows[ii], global_cols[jj]);
            }
          }
        } // get the neighbour
      } // walk the intersections
    } // walk the grid view
    pattern.sort();
    return pattern;
  } // ... compute_face_and_volume_pattern(...)

  PatternType compute_face_pattern() const
  {
    return compute_face_pattern(grid_view(), *this);
  }

  template< class G >
  PatternType compute_face_pattern(const GridView< G >& local_grid_view) const
  {
    return compute_face_pattern(local_grid_view, *this);
  }

  template< class S, size_t d, size_t r, size_t rC >
  PatternType compute_face_pattern(const SpaceInterface< S, d, r, rC >& ansatz_space) const
  {
    return compute_face_pattern(grid_view(), ansatz_space);
  }

  template< class G, class S, size_t d, size_t r, size_t rC >
  PatternType compute_face_pattern(const /*GridView<*/ G /*>*/& local_grid_view,
                                   const SpaceInterface< S, d, r, rC >& ansatz_space) const
  {
    // prepare
    PatternType pattern(mapper().size());
    Dune::DynamicVector< size_t > global_rows(mapper().maxNumDofs(), 0);
    Dune::DynamicVector< size_t > global_cols(ansatz_space.mapper().maxNumDofs(), 0);
    for (const auto& entity : DSC::entityRange(local_grid_view)) {
      const auto test_base_entity = base_function_set(entity);
      mapper().globalIndices(entity, global_rows);
      // walk the intersections
      const auto intersection_it_end = local_grid_view.iend(entity);
      for (auto intersection_it = local_grid_view.ibegin(entity);
           intersection_it != intersection_it_end;
           ++intersection_it) {
        const auto& intersection = *intersection_it;
        // get the neighbour
        if (intersection.neighbor() && !intersection.boundary()) {
          const auto neighbour_ptr = intersection.outside();
          const auto& neighbour = *neighbour_ptr;
          // get the basis
          const auto ansatz_base_neighbour = ansatz_space.base_function_set(neighbour);
          ansatz_space.mapper().globalIndices(neighbour, global_cols);
          // compute entity/neighbour
          for (size_t ii = 0; ii < test_base_entity.size(); ++ii) {
            for (size_t jj = 0; jj < ansatz_base_neighbour.size(); ++jj) {
              pattern.insert(global_rows[ii], global_cols[jj]);
            }
          }
        } // get the neighbour
      } // walk the intersections
    } // walk the grid view
    pattern.sort();
    return pattern;
  } // ... compute_face_pattern(...)

private:
  class BasisVisualization
    : public Dune::VTKFunction< GridViewType >
  {
    static_assert(dimRangeCols == 1, "Not implemented for matrixvalued spaces yet!");
  public:
    typedef typename BaseFunctionSetType::DomainType  DomainType;
    typedef typename BaseFunctionSetType::RangeType   RangeType;

    BasisVisualization(const derived_type& sp, const size_t ind, const std::string nm = "basis")
      : space_(sp)
      , index_(ind)
      , values_(space_.mapper().maxNumDofs(), RangeType(0))
      , name_(nm)
    {
      if (index_ >= space_.mapper().maxNumDofs())
        DUNE_THROW(Dune::RangeError,
                   "index has to be smaller than " << space_.mapper().maxNumDofs() << "(is " << index_ << ")!");
    }

    virtual std::string name() const
    {
      return name_;
    }

    /** \defgroup vtk ´´Methods to comply with the Dune::VTKFunction interface.'' */
    /* @{ */
    virtual int ncomps() const
    {
      return dimRange;
    }

    virtual double evaluate(int component, const EntityType& entity, const DomainType& xx) const
    {
      const auto baseFunctionSet = space_.base_function_set(entity);
      if (component < 0)
        DUNE_THROW(Dune::RangeError, "component must not be negative (is " << component << ")!");
      if (component < boost::numeric_cast< int >(baseFunctionSet.size())) {
        baseFunctionSet.evaluate(xx, values_);
        assert(component < boost::numeric_cast< int >(values_.size()) && "This should not happen!");
        return values_[index_][component];
      } else if (component < boost::numeric_cast< int >(space_.mapper().maxNumDofs()))
        return 0.0;
      else
        DUNE_THROW(Dune::RangeError,
                   "component has to be smaller than " << space_.mapper().maxNumDofs() << "(is " << component << ")!");
    }
    /* @} */

  private:
    const derived_type& space_;
    const size_t index_;
    mutable std::vector< RangeType > values_;
    const std::string name_;
  }; // class BasisVisualization

public:
  void visualize(const std::string filename_prefix = "") const
  {
    std::string filename = filename_prefix;
    if (filename.empty()) {
      filename = "dune.gdt.space";
    }
    VTKWriter< GridViewType >  vtk_writer(*(grid_view()), Dune::VTK::nonconforming);
    for (size_t ii = 0; ii < mapper().maxNumDofs(); ++ii) {
      std::string number = "";
      if (ii == 1)
        number = "1st";
      else if (ii == 2)
        number = "2nd";
      else if (ii == 3)
        number = "3rd";
      else
        number = Stuff::Common::toString(ii) + "th";
      const auto iith_baseFunction = std::make_shared< BasisVisualization >(this->as_imp(*this), ii, number + " basis");
      vtk_writer.addVertexData(iith_baseFunction);
    }
    vtk_writer.write(filename);
  } // ... visualize(...)

  /* @} */
}; // class SpaceInterface


template< class Traits, size_t d, size_t r, size_t rC, size_t codim = 0 >
typename Traits::GridViewType::template Codim< codim >::Iterator begin(const Dune::GDT::SpaceInterface< Traits, d, r, rC >& space)
{
  return space.grid_view().template begin< codim >();
}

template< class Traits, size_t d, size_t r, size_t rC, size_t codim = 0 >
typename Traits::GridViewType::template Codim< codim >::Iterator end(const Dune::GDT::SpaceInterface< Traits, d, r, rC >& space)
{
  return space.grid_view().template end< codim >();
}


namespace internal {


template< class S >
struct is_space_helper
{
  DSC_has_typedef_initialize_once(Traits)
  DSC_has_static_member_initialize_once(dimDomain)
  DSC_has_static_member_initialize_once(dimRange)
  DSC_has_static_member_initialize_once(dimRangeCols)

  static const bool is_candidate = DSC_has_typedef(Traits)< S >::value
                                   && DSC_has_static_member(dimDomain)< S >::value
                                   && DSC_has_static_member(dimRange)< S >::value
                                   && DSC_has_static_member(dimRangeCols)< S >::value;
}; // class is_space_helper


} // namespace internal


template< class S, bool candidate = internal::is_space_helper< S >::is_candidate >
struct is_space
  : public std::is_base_of< SpaceInterface< typename S::Traits, S::dimDomain, S::dimRange, S::dimRangeCols >, S >
{};


template< class S >
struct is_space< S, false >
  : public std::false_type
{};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_INTERFACE_HH
