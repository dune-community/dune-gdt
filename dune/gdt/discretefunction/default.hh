// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Kirsten Weber

#ifndef DUNE_GDT_DISCRETEFUNCTION_DEFAULT_HH
#define DUNE_GDT_DISCRETEFUNCTION_DEFAULT_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#elif defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/io/file/vtk/function.hh>

#include <dune/stuff/la/container/interface.hh>
#include <dune/gdt/space/interface.hh>

#include "local.hh"

namespace Dune {
namespace GDT {


template <class SpaceImp, class VectorImp>
class DiscreteFunctionDefaultConst : public Dune::VTKFunction<typename SpaceImp::GridPartType::GridViewType>,
                                     public Dune::Stuff::LocalizableFunction
{
public:
  typedef typename SpaceInterface<typename SpaceImp::Traits>::derived_type SpaceType;
  typedef typename Dune::Stuff::LA::VectorInterface<typename VectorImp::Traits>::derived_type VectorType;

  typedef typename SpaceType::EntityType EntityType;

  typedef DiscreteFunctionLocalConst<SpaceType, VectorType> ConstLocalFunctionType;
  typedef typename ConstLocalFunctionType::DomainType DomainType;
  typedef typename ConstLocalFunctionType::RangeType RangeType;

  DiscreteFunctionDefaultConst(const SpaceType& sp, const std::shared_ptr<const VectorType> vec,
                               const std::string nm = "discrete_function")
    : space_(sp)
    , name_(nm)
    , vector_(vec)
  {
    assert(vector_->size() == space_.mapper().size() && "Given vector has wrong size!");
  }

  const SpaceType& space() const
  {
    return space_;
  }

  virtual std::string name() const
  {
    return name_;
  }

  ConstLocalFunctionType localFunction(const EntityType& entity) const
  {
    return ConstLocalFunctionType(*this, entity);
  }

  std::shared_ptr<const VectorType> vector() const
  {
    return vector_;
  }

  /** \defgroup vtk ´´Methods to comply with the Dune::VTKFunction interface.'' */
  /* @{ */
  virtual int ncomps() const
  {
    return SpaceType::dimRange;
  }

  virtual double evaluate(int component, const EntityType& entity, const DomainType& x) const
  {
    RangeType ret(0);
    localFunction(entity).evaluate(x, ret);
    return ret[component];
  }
  /* @} */

private:
  const SpaceType& space_;
  const std::string name_;
  const std::shared_ptr<const VectorType> vector_;
}; // class DiscreteFunctionDefaultConst


template <class SpaceImp, class VectorImp>
class DiscreteFunctionDefault : public DiscreteFunctionDefaultConst<SpaceImp, VectorImp>
{
  typedef DiscreteFunctionDefaultConst<SpaceImp, VectorImp> BaseType;

public:
  typedef typename SpaceInterface<typename SpaceImp::Traits>::derived_type SpaceType;
  typedef typename Dune::Stuff::LA::VectorInterface<typename VectorImp::Traits>::derived_type VectorType;

  typedef typename SpaceType::EntityType EntityType;

  typedef DiscreteFunctionLocal<SpaceType, VectorType> LocalFunctionType;
  typedef typename LocalFunctionType::DomainType DomainType;

  DiscreteFunctionDefault(const SpaceType& sp, std::shared_ptr<VectorType> vec,
                          const std::string nm = "discrete_function")
    : BaseType(sp, vec, nm)
    , nonConstVector_(vec)
  {
  }

  DiscreteFunctionDefault(const SpaceType& sp, const std::string nm = "discrete_function")
    : BaseType(sp, std::make_shared<VectorType>(sp.mapper().size()), nm)
    , nonConstVector_(std::make_shared<VectorType>(const_cast<VectorType&>(
          *(BaseType::vector())))) // <- I am not sure if this is safe to do. I would expect the BaseType to live long
  // enough, but maybe someone else has a better Idea!
  {
  }

  LocalFunctionType localFunction(const EntityType& entity)
  {
    return LocalFunctionType(*this, entity);
  }

  std::shared_ptr<VectorType> vector()
  {
    return nonConstVector_;
  }

  //  virtual std::string name() const
  //  {
  //    return BaseType::name();
  //  }

  //  /** \defgroup vtk ´´Methods to comply with the Dune::VTKFunction interface.'' */
  //  /* @{ */
  //  virtual int ncomps() const
  //  {
  //    return BaseType::ncomps();
  //  }

  //  virtual double evaluate(int component, const EntityType& entity, const DomainType& x) const
  //  {
  //    BaseType::evaluate(component, entity, x);
  //  }
  //  /* @} */

private:
  std::shared_ptr<VectorType> nonConstVector_;
}; // class DiscreteFunctionDefault


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETEFUNCTION_DEFAULT_HH
