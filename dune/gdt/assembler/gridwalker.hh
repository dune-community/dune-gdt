// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_ASSEMBLER_GRIDWALKER_HH
#define DUNE_GDT_ASSEMBLER_GRIDWALKER_HH

#include <vector>
#include <memory>
#include <type_traits>

#include <dune/stuff/grid/boundaryinfo.hh>

namespace Dune {
namespace GDT {
namespace Functor {


template <class GridPartImp>
class Codim0
{
public:
  typedef GridPartImp GridPartType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;

  virtual ~Codim0()
  {
  }

  virtual void prepare()
  {
  }

  virtual void apply_local(const EntityType& entity) = 0;

  virtual void finalize()
  {
  }
}; // class Codim0


template <class GridPartImp>
class Codim1
{
public:
  typedef GridPartImp GridPartType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef typename GridPartType::IntersectionType IntersectionType;

  virtual ~Codim1()
  {
  }

  virtual void prepare()
  {
  }

  virtual void apply_local(const IntersectionType& intersection, const EntityType& inside_entity,
                           const EntityType& outside_entity) = 0;

  virtual void finalize()
  {
  }
}; // class Codim1


/**
 *  \brief Interface for functors to tell on which entity to apply.
 *
 *  The derived class has to provide a method with the following signature:
 *  \code
virtual bool apply_on(const GridPartType& grid_part, const EntityType& entity) const
{
  ...
}
\endcode
 */
template <class GridPartImp>
class ApplyOnWhichEntity
{
public:
  typedef GridPartImp GridPartType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;

  virtual ~ApplyOnWhichEntity()
  {
  }

  virtual bool apply_on(const GridPartType& /*grid_part*/, const EntityType& /*entity*/) const = 0;
}; // class ApplyOnWhichEntity


/**
 *  \brief Selects all entities.
 */
template <class GridPartImp>
class ApplyOnAllEntities : public ApplyOnWhichEntity<GridPartImp>
{
public:
  typedef GridPartImp GridPartType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;

  virtual bool apply_on(const GridPartType& /*grid_part*/, const EntityType& /*entity*/) const /*DS_OVERRIDE*/
  {
    return true;
  }
}; // class ApplyOnAllEntities


/**
 *  \brief Interface for functors to tell on which intersection to apply.
 *
 *  The derived class has to provide a method with the following signature:
 *  \code
virtual bool apply_on(const GridPartType& grid_part, const IntersectionType& intersection) const
{
  ...
}
\endcode
 */
template <class GridPartImp>
class ApplyOnWhichIntersection
{
public:
  typedef GridPartImp GridPartType;
  typedef typename GridPartType::IntersectionType IntersectionType;

  virtual ~ApplyOnWhichIntersection<GridPartImp>()
  {
  }

  virtual bool apply_on(const GridPartType& /*grid_part*/, const IntersectionType& /*intersection*/) const = 0;
}; // class ApplyOnWhichIntersection< GridPartImp >


/**
 *  \brief Selects each inner intersection.
 *
 *  To decide if this in an inner intersection,
\code
intersection.neighbor() && !intersection.boundary()
\endcode
 *  is used.
 */
template <class GridPartImp>
class ApplyOnInnerIntersections : public ApplyOnWhichIntersection<GridPartImp>
{
public:
  typedef GridPartImp GridPartType;
  typedef typename GridPartType::IntersectionType IntersectionType;

  virtual bool apply_on(const GridPartType& /*grid_part*/, const IntersectionType& intersection) const /*DS_OVERRIDE*/
  {
    return intersection.neighbor() && !intersection.boundary();
  }
}; // class ApplyOnInnerIntersections


/**
 *  \brief Selects each inner intersection only once.
 *
 *  To decide if this in an inner intersection,
\code
intersection.neighbor() && !intersection.boundary()
\endcode
 *  is used, and true is returned, if the index of the inside() entity is smaller than the index of the outside()
 *  entity.
 */
template <class GridPartImp>
class ApplyOnInnerIntersectionsPrimally : public ApplyOnWhichIntersection<GridPartImp>
{
public:
  typedef GridPartImp GridPartType;
  typedef typename GridPartType::IntersectionType IntersectionType;

  virtual bool apply_on(const GridPartType& grid_part, const IntersectionType& intersection) const /*DS_OVERRIDE*/
  {
    if (intersection.neighbor() && !intersection.boundary()) {
      const auto insideEntityPtr    = intersection.inside();
      const auto& insideEntity      = *insideEntityPtr;
      const auto outsideNeighborPtr = intersection.outside();
      const auto& outsideNeighbor = *outsideNeighborPtr;
      return grid_part.indexSet().index(insideEntity) < grid_part.indexSet().index(outsideNeighbor);
    } else
      return false;
  }
}; // class ApplyOnInnerIntersections


template <class GridPartImp>
class ApplyOnBoundaryIntersections : public ApplyOnWhichIntersection<GridPartImp>
{
public:
  typedef GridPartImp GridPartType;
  typedef typename GridPartType::IntersectionType IntersectionType;

  virtual bool apply_on(const GridPartType& /*grid_part*/, const IntersectionType& intersection) const /*DS_OVERRIDE*/
  {
    return intersection.boundary();
  }
}; // class ApplyOnBoundaryIntersections


template <class GridPartImp>
class ApplyOnDirichletIntersections : public ApplyOnWhichIntersection<GridPartImp>
{
public:
  typedef GridPartImp GridPartType;
  typedef typename GridPartType::IntersectionType IntersectionType;
  typedef Stuff::GridboundaryInterface<IntersectionType> BoundaryInfoType;

  ApplyOnDirichletIntersections(const BoundaryInfoType& boundary_info)
    : boundary_info_(boundary_info)
  {
  }

  virtual bool apply_on(const GridPartType& /*grid_part*/, const IntersectionType& intersection) const /*DS_OVERRIDE*/
  {
    return boundary_info_.dirichlet(intersection);
  }

private:
  const BoundaryInfoType& boundary_info_;
}; // class ApplyOnDirichletIntersections


template <class GridPartImp>
class ApplyOnNeumannIntersections : public ApplyOnWhichIntersection<GridPartImp>
{
public:
  typedef GridPartImp GridPartType;
  typedef typename GridPartType::IntersectionType IntersectionType;
  typedef Stuff::GridboundaryInterface<IntersectionType> BoundaryInfoType;

  ApplyOnNeumannIntersections(const BoundaryInfoType& boundary_info)
    : boundary_info_(boundary_info)
  {
  }

  virtual bool apply_on(const GridPartType& /*grid_part*/, const IntersectionType& intersection) const /*DS_OVERRIDE*/
  {
    return boundary_info_.neumann(intersection);
  }

private:
  const BoundaryInfoType& boundary_info_;
}; // class ApplyOnNeumannIntersections


} // namespace Functor


template <class GridPartImp>
class GridWalker
{
  typedef GridWalker<GridPartImp> ThisType;

public:
  typedef GridPartImp GridPartType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef typename GridPartType::IntersectionType IntersectionType;

private:
  typedef Stuff::GridboundaryInterface<IntersectionType> BoundaryInfoType;
  typedef Functor::Codim0<GridPartType> Codim0FunctorType;
  typedef Functor::Codim1<GridPartType> Codim1FunctorType;

  class Codim0FunctorWrapper : public Functor::Codim0<GridPartType>
  {
  public:
    Codim0FunctorWrapper(Codim0FunctorType& wrapped_functor)
      : wrapped_functor_(wrapped_functor)
    {
    }

    virtual void prepare() DS_FINAL
    {
      wrapped_functor_.prepare();
    }

    virtual void apply_local(const EntityType& entity) DS_FINAL
    {
      wrapped_functor_.apply_local(entity);
    }

    virtual void finalize() DS_FINAL
    {
      wrapped_functor_.finalize();
    }

  private:
    Codim0FunctorType& wrapped_functor_;
  }; // class Codim0FunctorWrapper

  class Codim1FunctorWrapper : public Functor::Codim1<GridPartType>
  {
  public:
    Codim1FunctorWrapper(Codim1FunctorType& wrapped_functor)
      : wrapped_functor_(wrapped_functor)
    {
    }

    virtual void prepare() DS_FINAL
    {
      wrapped_functor_.prepare();
    }

    virtual void apply_local(const IntersectionType& intersection, const EntityType& inside_entity,
                             const EntityType& outside_entity) DS_FINAL
    {
      wrapped_functor_.apply_local(intersection, inside_entity, outside_entity);
    }

    virtual void finalize() DS_FINAL
    {
      wrapped_functor_.finalize();
    }

  private:
    Codim1FunctorType& wrapped_functor_;
  }; // class Codim1FunctorWrapper

public:
  GridWalker(const GridPartType& grid_part)
    : grid_part_(grid_part)
  {
  }

  void add(Functor::Codim0<GridPartType>& functor)
  {
    codim0_functors_.emplace_back(new Codim0FunctorWrapper(functor));
  }

  void add(Functor::Codim1<GridPartType>& functor)
  {
    codim1_functors_.emplace_back(new Codim1FunctorWrapper(functor));
  }

  void walk() const
  {
    // prepare functors
    for (auto& functor : codim0_functors_)
      functor->prepare();

    // only do something, if we have to
    if ((codim0_functors_.size() + codim1_functors_.size()) > 0) {
      // walk the grid
      const auto entity_it_end = grid_part_.template end<0>();
      for (auto entity_it = grid_part_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
        const EntityType& entity = *entity_it;

        // apply codim0 functors
        for (auto& functor : codim0_functors_)
          functor->apply_local(entity);

        // only walk the intersections, if there are codim1 functors present
        if (codim1_functors_.size()) {
          // walk the intersections
          const auto intersection_it_end = grid_part_.iend(entity);
          for (auto intersection_it = grid_part_.ibegin(entity); intersection_it != intersection_it_end;
               ++intersection_it) {
            const auto& intersection   = *intersection_it;
            const auto neighbor_ptr    = intersection.outside();
            const EntityType& neighbor = *neighbor_ptr;

            // apply codim1 functors
            for (auto& functor : codim1_functors_)
              functor->apply_local(intersection, entity, neighbor);

          } // walk the intersections
        } // only walk the intersections, if there are codim1 functors present
      } // walk the grid
    } // only do something, if we have to

    // finalize functors
    for (auto& functor : codim0_functors_)
      functor->finalize();
  } // ... walk()

private:
  const GridPartType& grid_part_;
  std::vector<std::unique_ptr<Codim0FunctorType>> codim0_functors_;
  std::vector<std::unique_ptr<Codim1FunctorType>> codim1_functors_;
}; // class GridWalker


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_GRIDWALKER_HH
