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

#include <memory>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/deprecated.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/xt/common/crtp.hh>
#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/parallel/threadstorage.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/tuple.hh>
#include <dune/xt/common/fixed_map.hh>

#include <dune/xt/la/container/pattern.hh>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/view/from-part.hh>

#include <dune/gdt/local/finite-elements/interfaces.hh>
#include <dune/gdt/spaces/mapper/interfaces.hh>
#include <dune/gdt/spaces/basis/interface.hh>
#include <dune/gdt/spaces/parallel.hh>

#include "constraints.hh"

namespace Dune {
namespace GDT {


enum class Backends
{
  gdt
};

static const XT::Common::FixedMap<Backends, std::string, 1> backend_names = {{Backends::gdt, "gdt"}};

enum class SpaceType
{
  cg,
  block_cg,
  dg,
  block_dg,
  fv,
  product_fv,
  block_fv,
  rt,
  block_rt
};


namespace internal {


template <Backends backend>
struct backend_dependent_typename
{
  typedef void type;
};


template <SpaceType tp>
struct space_type_dependent_typename
{
  typedef void type;
};


} // namespace  internal


static constexpr Backends default_space_backend = Backends::gdt;

enum class ChoosePattern
{
  volume,
  face,
  face_and_volume
};


template <Backends type>
struct layer_from_backend;

template <>
struct layer_from_backend<Backends::gdt>
{
  static const XT::Grid::Backends type = XT::Grid::Backends::view;
};


template <class Traits, size_t domainDim, size_t rangeDim, size_t rangeDimCols = 1>
class SpaceInterface : public XT::CRTPInterface<SpaceInterface<Traits, domainDim, rangeDim, rangeDimCols>, Traits>
{
public:
  typedef typename Traits::derived_type derived_type;
  static const int polOrder = Traits::polOrder;
  static const bool continuous = Traits::continuous;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::DofCommunicatorType DofCommunicatorType;
  typedef typename Traits::GridLayerType GridLayerType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const size_t dimDomain = domainDim;
  static const size_t dimRange = rangeDim;
  static const size_t dimRangeCols = rangeDimCols;
  static const constexpr Backends backend_type{Traits::backend_type};

private:
  static_assert(dimDomain > 0, "dimDomain has to be positive");
  static_assert(dimRange > 0, "dimRange has to be positive");
  static_assert(dimRangeCols > 0, "dimRangeCols has to be positive");
  static_assert(XT::Grid::is_layer<GridLayerType>::value, "");
  static_assert(GridLayerType::dimension == dimDomain, "Dimension of GridLayerType has to match dimDomain");

public:
  typedef typename GridLayerType::ctype DomainFieldType;
  typedef FieldVector<DomainFieldType, dimDomain> DomainType;

  using EntityType = XT::Grid::extract_entity_t<GridLayerType>;
  typedef Dune::XT::LA::SparsityPatternDefault PatternType;

  static const XT::Grid::Backends layer_backend = Traits::layer_backend;

  using GlobalBasisType = GlobalBasisInterface<EntityType, dimRange, dimRangeCols, RangeFieldType>;
  using MapperType = MapperInterface<GridLayerType>;
  using FiniteElementType =
      LocalFiniteElementInterface<typename GridLayerType::ctype, dimDomain, RangeFieldType, dimRange, dimRangeCols>;

  /**
   * \defgroup interface ´´These methods have to be implemented!''
   * @{
   **/

  const GridLayerType& grid_layer() const
  {
    CHECK_CRTP(this->as_imp().grid_layer());
    return this->as_imp().grid_layer();
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

  const GlobalBasisType& basis() const
  {
    CHECK_CRTP(this->as_imp().basis());
    return this->as_imp().basis();
  }

  DofCommunicatorType& dof_communicator() const
  {
    CHECK_CRTP(this->as_imp().dof_communicator());
    return this->as_imp().dof_communicator();
  }

  //! communication data handles may require to know this to setup buffers and trnasmission patterns
  static constexpr bool associates_data_with(int codim)
  {
    return derived_type::associates_data_with(codim);
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
  template <class S, size_t d, size_t r, size_t rC, class C>
  void local_constraints(const SpaceInterface<S, d, r, rC>& ansatz_space,
                         const EntityType& entity,
                         ConstraintsInterface<C>& ret) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().local_constraints(ansatz_space.as_imp(), entity, ret.as_imp()));
  }

  template <class C>
  void local_constraints(const EntityType& entity, ConstraintsInterface<C>& ret) const
  {
    local_constraints(*this, entity, ret);
  }

  /**
   *  \brief  Computes the appropriate sparsity pattern.
   *  \note   This method can be implemented in a derived class by a forward to one of the methods provided by this
   * class, namely compute_volume_pattern(), compute_face_pattern() or compute_face_and_volume_pattern().
   */
  template <class GL, class S, size_t d, size_t r, size_t rC>
  typename std::enable_if<XT::Grid::is_layer<GL>::value, PatternType>::type
  compute_pattern(const GL& grd_layr, const SpaceInterface<S, d, r, rC>& ansatz_space) const
  {
    CHECK_CRTP(this->as_imp().compute_pattern(grd_layr, ansatz_space.as_imp()));
    return this->as_imp().compute_pattern(grd_layr, ansatz_space.as_imp());
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

  template <class S, size_t d, size_t r, size_t rC>
  PatternType compute_pattern(const SpaceInterface<S, d, r, rC>& ansatz_space) const
  {
    return compute_pattern(grid_layer(), ansatz_space);
  }

  template <class GL>
  typename std::enable_if<XT::Grid::is_layer<GL>::value, PatternType>::type compute_pattern(const GL& grd_layr) const
  {
    return compute_pattern(grd_layr, *this);
  }

  PatternType compute_volume_pattern() const
  {
    return compute_volume_pattern(*this);
  }

  template <class S, size_t d, size_t r, size_t rC>
  PatternType compute_volume_pattern(const SpaceInterface<S, d, r, rC>& ansatz_space) const
  {
    return compute_volume_pattern(grid_layer(), ansatz_space);
  }

  template <class GL>
  typename std::enable_if<XT::Grid::is_layer<GL>::value, PatternType>::type
  compute_volume_pattern(const GL& grd_layr) const
  {
    return compute_volume_pattern(grd_layr, *this);
  }

  /**
   *  \brief  computes a sparsity pattern, where this space is the test space (rows/outer) and the other space is the
   *          ansatz space (cols/inner)
   */
  template <class GL, class S, size_t d, size_t r, size_t rC>
  typename std::enable_if<XT::Grid::is_layer<GL>::value, PatternType>::type
  compute_volume_pattern(const GL& grd_layr, const SpaceInterface<S, d, r, rC>& ansatz_space) const
  {
    PatternType pattern(mapper().size());
    Dune::DynamicVector<size_t> globalRows(mapper().maxNumDofs(), 0);
    Dune::DynamicVector<size_t> globalCols(ansatz_space.mapper().maxNumDofs(), 0);

    for (const auto& entity : elements(grd_layr)) {
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
    return compute_face_and_volume_pattern(grid_layer(), *this);
  }

  template <class GL>
  typename std::enable_if<XT::Grid::is_layer<GL>::value, PatternType>::type
  compute_face_and_volume_pattern(const GL& grd_layr) const
  {
    return compute_face_and_volume_pattern(grd_layr, *this);
  }

  template <class S, size_t d, size_t r, size_t rC>
  PatternType compute_face_and_volume_pattern(const SpaceInterface<S, d, r, rC>& ansatz_space) const
  {
    return compute_face_and_volume_pattern(grid_layer(), ansatz_space);
  }

  /**
   *  \brief  computes a DG sparsity pattern, where this space is the test space (rows/outer) and the other space is the
   *          ansatz space (cols/inner)
   */
  template <class GL, class S, size_t d, size_t r, size_t rC>
  typename std::enable_if<XT::Grid::is_layer<GL>::value, PatternType>::type
  compute_face_and_volume_pattern(const GL& grd_layr, const SpaceInterface<S, d, r, rC>& ansatz_space) const
  {
    // prepare
    PatternType pattern(mapper().size());
    Dune::DynamicVector<size_t> global_rows(mapper().maxNumDofs(), 0);
    Dune::DynamicVector<size_t> global_cols(ansatz_space.mapper().maxNumDofs(), 0);
    for (const auto& entity : elements(grd_layr)) {
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
      const auto intersection_it_end = grd_layr.iend(entity);
      for (auto intersection_it = grd_layr.ibegin(entity); intersection_it != intersection_it_end; ++intersection_it) {
        const auto& intersection = *intersection_it;
        // get the neighbour
        if (intersection.neighbor() && !intersection.boundary()) {
          const auto neighbour = intersection.outside();
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
    } // walk the grid layer
    pattern.sort();
    return pattern;
  } // ... compute_face_and_volume_pattern(...)

  PatternType compute_face_pattern() const
  {
    return compute_face_pattern(grid_layer(), *this);
  }

  template <class GL>
  typename std::enable_if<XT::Grid::is_layer<GL>::value, PatternType>::type
  compute_face_pattern(const GL& grd_layr) const
  {
    return compute_face_pattern(grd_layr, *this);
  }

  template <class S, size_t d, size_t r, size_t rC>
  PatternType compute_face_pattern(const SpaceInterface<S, d, r, rC>& space) const
  {
    return compute_face_pattern(grid_layer(), space);
  }

  template <class GL, class S, size_t d, size_t r, size_t rC>
  typename std::enable_if<XT::Grid::is_layer<GL>::value, PatternType>::type
  compute_face_pattern(const GL& grd_layr, const SpaceInterface<S, d, r, rC>& ansatz_space) const
  {
    // prepare
    PatternType pattern(mapper().size());
    Dune::DynamicVector<size_t> global_rows(mapper().maxNumDofs(), 0);
    Dune::DynamicVector<size_t> global_cols(ansatz_space.mapper().maxNumDofs(), 0);
    for (const auto& entity : elements(grd_layr)) {
      const auto test_base_entity = base_function_set(entity);
      mapper().globalIndices(entity, global_rows);
      // walk the intersections
      const auto intersection_it_end = grd_layr.iend(entity);
      for (auto intersection_it = grd_layr.ibegin(entity); intersection_it != intersection_it_end; ++intersection_it) {
        const auto& intersection = *intersection_it;
        // get the neighbour
        if (intersection.neighbor() && !intersection.boundary()) {
          const auto neighbour = intersection.outside();
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
    } // walk the grid layer
    pattern.sort();
    return pattern;
  } // ... compute_face_pattern(...)

private:
  template <class GV>
  class BasisVisualization : public Dune::VTKFunction<GV>
  {
    static_assert(dimRangeCols == 1, "Not implemented for matrixvalued spaces yet!");

  public:
    typedef typename GlobalBasisType::DomainType DomainType;
    typedef typename GlobalBasisType::RangeType RangeType;

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
      if (component < boost::numeric_cast<int>(baseFunctionSet.size())) {
        baseFunctionSet.evaluate(xx, values_);
        assert(component < boost::numeric_cast<int>(values_.size()) && "This should not happen!");
        return values_[index_][component];
      } else if (component < boost::numeric_cast<int>(space_.mapper().maxNumDofs()))
        return 0.0;
      else
        DUNE_THROW(Dune::RangeError,
                   "component has to be smaller than " << space_.mapper().maxNumDofs() << "(is " << component << ")!");
    }
    /* @} */

  private:
    const derived_type& space_;
    const size_t index_;
    mutable std::vector<RangeType> values_;
    const std::string name_;
  }; // class BasisVisualization

public:
  void visualize(const std::string filename_prefix = "") const
  {
    const std::string filename = filename_prefix.empty() ? "dune.gdt.space" : filename_prefix;
    const auto tmp_storage = XT::Grid::make_tmp_view(grid_layer());
    const auto& grd_vw = tmp_storage.access();
    using GridViewType = std::decay_t<decltype(grd_vw)>;
    VTKWriter<GridViewType> vtk_writer(grd_vw, Dune::VTK::nonconforming);
    for (size_t ii = 0; ii < mapper().maxNumDofs(); ++ii) {
      std::string number = "";
      if (ii == 1)
        number = "1st";
      else if (ii == 2)
        number = "2nd";
      else if (ii == 3)
        number = "3rd";
      else
        number = XT::Common::to_string(ii) + "th";
      const auto iith_baseFunction =
          std::make_shared<BasisVisualization<GridViewType>>(this->as_imp(*this), ii, number + " basis");
      vtk_writer.addVertexData(iith_baseFunction);
    }
    vtk_writer.write(filename);
  } // ... visualize(...)

  /* @} */
}; // class SpaceInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_INTERFACE_HH
