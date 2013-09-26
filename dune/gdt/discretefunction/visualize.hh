// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_DISCRETEFUNCTION_VISUALIZE_HH
#define DUNE_GDT_DISCRETEFUNCTION_VISUALIZE_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#elif defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <memory>
#include <type_traits>

#include <dune/common/fvector.hh>
#include <dune/common/dynvector.hh>

#include <dune/grid/io/file/vtk/function.hh>

namespace Dune {
namespace GDT {
namespace DiscreteFunction {


template <class SpaceImp>
class BasisVisualization : public Dune::VTKFunction<typename SpaceImp::GridPartType::GridViewType>
{
public:
  typedef SpaceImp SpaceType;
  typedef typename SpaceType::EntityType EntityType;
  typedef typename SpaceType::BaseFunctionSetType::DomainType DomainType;

private:
  typedef typename SpaceType::BaseFunctionSetType::RangeType RangeType;

public:
  BasisVisualization(const SpaceType& sp, const size_t ind, const std::string nm = "basis")
    : space_(sp)
    , index_(ind)
    , values_(space_.mapper().maxNumDofs(), RangeType(0))
    , name_(nm)
  {
    if (index_ >= space_.mapper().maxNumDofs())
      DUNE_THROW(Dune::RangeError,
                 "index has to be smaller than " << space_.mapper().maxNumDofs() << "(is " << index_ << ")!");
  }

  const SpaceType& space() const
  {
    return space_;
  }

  virtual std::string name() const
  {
    return name_;
  }

  /** \defgroup vtk ´´Methods to comply with the Dune::VTKFunction interface.'' */
  /* @{ */
  virtual int ncomps() const
  {
    return SpaceType::dimRange;
  }

  virtual double evaluate(int component, const EntityType& entity, const DomainType& x) const
  {
    const auto baseFunctionSet = space_.baseFunctionSet(entity);
    if (component < 0)
      DUNE_THROW(Dune::RangeError, "component must not be negative (is " << component << ")!");
    if (component < baseFunctionSet.size()) {
      baseFunctionSet.evaluate(x, values_);
      assert(component < values_.size() && "This should not happen!");
      return values_[index_][component];
    } else if (component < space_.mapper().maxNumDofs())
      return 0.0;
    else
      DUNE_THROW(Dune::RangeError,
                 "component has to be smaller than " << space_.mapper().maxNumDofs() << "(is " << component << ")!");
  }
  /* @} */

private:
  const SpaceType& space_;
  const size_t index_;
  mutable std::vector<RangeType> values_;
  const std::string name_;
}; // class BasisVisualization


template <class SpaceImp, class VectorImp>
class VisualizationAdapter : public Dune::VTKFunction<typename SpaceImp::GridPartType::GridViewType>
{
public:
  typedef SpaceImp SpaceType;
  typedef VectorImp VectorType;

  typedef typename SpaceType::EntityType EntityType;

  typedef typename SpaceType::BaseFunctionSetType::DomainType DomainType;
  typedef typename SpaceType::BaseFunctionSetType::RangeType RangeType;

  VisualizationAdapter(const SpaceType& sp, const VectorType& vec, const std::string nm = "discrete_function")
    : space_(sp)
    , name_(nm)
    , vector_(vec)
    , indices_(space_.mapper().maxNumDofs(), 0)
    , basis_values_(space_.mapper().maxNumDofs(), RangeType(0))
  {
    if (vector_.size() != space_.mapper().size())
      DUNE_THROW(Dune::RangeError,
                 "The sice of vec (" << vec.size() << ") does not match the size of the space ("
                                     << space_.mapper().size()
                                     << ")!");
  }

  ~VisualizationAdapter()
  {
  }

  const SpaceType& space() const
  {
    return space_;
  }

  virtual std::string name() const
  {
    return name_;
  }

  /** \defgroup vtk ´´Methods to comply with the Dune::VTKFunction interface.'' */
  /* @{ */
  virtual int ncomps() const
  {
    return SpaceType::dimRange;
  }

  virtual double evaluate(int component, const EntityType& entity, const DomainType& xx) const
  {
    if (component < 0)
      DUNE_THROW(Dune::RangeError, "component must not be negative (is " << component << ")!");
    RangeType value(0);
    if (component >= ncomps())
      DUNE_THROW(Dune::RangeError,
                 "component has to be smaller than ncomps() = " << ncomps() << " (is " << component << ")!");
    const auto baseFunctionSet = space_.baseFunctionSet(entity);
    space_.mapper().globalIndices(entity, indices_);
    baseFunctionSet.evaluate(xx, basis_values_);
    for (size_t ii = 0; ii < baseFunctionSet.size(); ++ii) {
      basis_values_[ii] *= vector_.get(indices_[ii]);
      value += basis_values_[ii];
    }
    return value[component];
  }
  /* @} */

private:
  const SpaceType& space_;
  const std::string name_;
  const VectorType& vector_;
  mutable Dune::DynamicVector<size_t> indices_;
  mutable std::vector<RangeType> basis_values_;
}; // class VisualizationAdapter


} // namespace DiscreteFunction
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETEFUNCTION_VISUALIZE_HH
