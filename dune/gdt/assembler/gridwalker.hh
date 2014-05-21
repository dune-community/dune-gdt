// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
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


template <class GridViewImp>
class Codim0
{
public:
  typedef GridViewImp GridViewType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;

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


template <class GridViewImp>
class Codim1
{
public:
  typedef GridViewImp GridViewType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::Intersection IntersectionType;

  virtual ~Codim1()
  {
  }

  virtual void prepare()
  {
  }

  virtual void apply_local(const IntersectionType& /*intersection*/, const EntityType& /*inside_entity*/,
                           const EntityType& /*outside_entity*/) = 0;

  virtual void finalize()
  {
  }
}; // class Codim1


template <class GridViewImp>
class Codim0And1
{
public:
  typedef GridViewImp GridViewType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::Intersection IntersectionType;

  virtual ~Codim0And1()
  {
  }

  virtual void prepare()
  {
  }

  virtual void apply_local(const EntityType& entity) = 0;

  virtual void apply_local(const IntersectionType& /*intersection*/, const EntityType& /*inside_entity*/,
                           const EntityType& /*outside_entity*/) = 0;

  virtual void finalize()
  {
  }
}; // class Codim0And1


} // namespace Functor
namespace ApplyOn {


/**
 *  \brief Interface for functors to tell on which entity to apply.
 *
 *  The derived class has to provide a method with the following signature:
 *  \code
virtual bool apply_on(const GridViewType& grid_view, const EntityType& entity) const
{
  ...
}
\endcode
 */
template <class GridViewImp>
class WhichEntity
{
public:
  typedef GridViewImp GridViewType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;

  virtual ~WhichEntity()
  {
  }

  virtual bool apply_on(const GridViewType& /*grid_view*/, const EntityType& /*entity*/) const = 0;
}; // class WhichEntity


/**
 *  \brief Selects all entities.
 */
template <class GridViewImp>
class AllEntities : public WhichEntity<GridViewImp>
{
public:
  typedef GridViewImp GridViewType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;

  virtual bool apply_on(const GridViewType& /*grid_view*/, const EntityType& /*entity*/) const DS_OVERRIDE DS_FINAL
  {
    return true;
  }
}; // class AllEntities


/**
 *  \brief Selects entities which have a boundary intersection.
 */
template <class GridViewImp>
class BoundaryEntities : public WhichEntity<GridViewImp>
{
public:
  typedef GridViewImp GridViewType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;

  virtual bool apply_on(const GridViewType& /*grid_view*/, const EntityType& entity) const DS_OVERRIDE DS_FINAL
  {
    return entity.hasBoundaryIntersections();
  }
}; // class BoundaryEntities


/**
 *  \brief Interface for functors to tell on which intersection to apply.
 *
 *  The derived class has to provide a method with the following signature:
 *  \code
virtual bool apply_on(const GridViewType& grid_view, const IntersectionType& intersection) const
{
  ...
}
\endcode
 */
template <class GridViewImp>
class WhichIntersection
{
public:
  typedef GridViewImp GridViewType;
  typedef typename GridViewType::Intersection IntersectionType;

  virtual ~WhichIntersection<GridViewImp>()
  {
  }

  virtual bool apply_on(const GridViewType& /*grid_view*/, const IntersectionType& /*intersection*/) const = 0;
}; // class WhichIntersection< GridViewImp >


/**
 *  \brief Selects all intersections.
 */
template <class GridViewImp>
class AllIntersections : public WhichIntersection<GridViewImp>
{
public:
  typedef GridViewImp GridViewType;
  typedef typename GridViewType::Intersection IntersectionType;

  virtual bool apply_on(const GridViewType& /*grid_view*/,
                        const IntersectionType& /*intersection*/) const DS_OVERRIDE DS_FINAL
  {
    return true;
  }
}; // class AllIntersections


/**
 *  \brief Selects each inner intersection.
 *
 *  To decide if this in an inner intersection,
\code
intersection.neighbor() && !intersection.boundary()
\endcode
 *  is used.
 */
template <class GridViewImp>
class InnerIntersections : public WhichIntersection<GridViewImp>
{
public:
  typedef GridViewImp GridViewType;
  typedef typename GridViewType::Intersection IntersectionType;

  virtual bool apply_on(const GridViewType& /*grid_view*/,
                        const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
  {
    return intersection.neighbor() && !intersection.boundary();
  }
}; // class InnerIntersections


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
template <class GridViewImp>
class InnerIntersectionsPrimally : public WhichIntersection<GridViewImp>
{
public:
  typedef GridViewImp GridViewType;
  typedef typename GridViewType::Intersection IntersectionType;

  virtual bool apply_on(const GridViewType& grid_view, const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
  {
    if (intersection.neighbor() && !intersection.boundary()) {
      const auto insideEntityPtr    = intersection.inside();
      const auto& insideEntity      = *insideEntityPtr;
      const auto outsideNeighborPtr = intersection.outside();
      const auto& outsideNeighbor = *outsideNeighborPtr;
      return grid_view.indexSet().index(insideEntity) < grid_view.indexSet().index(outsideNeighbor);
    } else
      return false;
  }
}; // class InnerIntersections


template <class GridViewImp>
class BoundaryIntersections : public WhichIntersection<GridViewImp>
{
public:
  typedef GridViewImp GridViewType;
  typedef typename GridViewType::Intersection IntersectionType;

  virtual bool apply_on(const GridViewType& /*grid_view*/,
                        const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
  {
    return intersection.boundary();
  }
}; // class BoundaryIntersections


template <class GridViewImp>
class DirichletIntersections : public WhichIntersection<GridViewImp>
{
public:
  typedef GridViewImp GridViewType;
  typedef typename GridViewType::Intersection IntersectionType;
  typedef Stuff::Grid::BoundaryInfoInterface<IntersectionType> BoundaryInfoType;

  DirichletIntersections(const BoundaryInfoType& boundary_info)
    : boundary_info_(boundary_info)
  {
  }

  virtual bool apply_on(const GridViewType& /*grid_view*/,
                        const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
  {
    return boundary_info_.dirichlet(intersection);
  }

private:
  const BoundaryInfoType& boundary_info_;
}; // class DirichletIntersections


template <class GridViewImp>
class NeumannIntersections : public WhichIntersection<GridViewImp>
{
public:
  typedef GridViewImp GridViewType;
  typedef typename GridViewType::Intersection IntersectionType;
  typedef Stuff::Grid::BoundaryInfoInterface<IntersectionType> BoundaryInfoType;

  NeumannIntersections(const BoundaryInfoType& boundary_info)
    : boundary_info_(boundary_info)
  {
  }

  virtual bool apply_on(const GridViewType& /*grid_view*/,
                        const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
  {
    return boundary_info_.neumann(intersection);
  }

private:
  const BoundaryInfoType& boundary_info_;
}; // class NeumannIntersections


} // namespace ApplyOn


template <class GridViewImp>
class GridWalker : public Functor::Codim0And1<GridViewImp>
{
  typedef GridWalker<GridViewImp> ThisType;

public:
  typedef GridViewImp GridViewType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::Intersection IntersectionType;

protected:
  typedef Stuff::Grid::BoundaryInfoInterface<IntersectionType> BoundaryInfoType;

  class Codim0Object : public Functor::Codim0<GridViewType>
  {
  public:
    ~Codim0Object()
    {
    }
    virtual bool apply_on(const GridViewType& grid_view, const EntityType& entity) const = 0;
  };


  template <class Codim0FunctorType>
  class Codim0FunctorWrapper : public Codim0Object
  {
  public:
    Codim0FunctorWrapper(Codim0FunctorType& wrapped_functor, const ApplyOn::WhichEntity<GridViewType>* where)
      : wrapped_functor_(wrapped_functor)
      , where_(where)
    {
    }

    virtual ~Codim0FunctorWrapper()
    {
    }

    virtual void prepare() DS_OVERRIDE DS_FINAL
    {
      wrapped_functor_.prepare();
    }

    virtual bool apply_on(const GridViewType& grid_view, const EntityType& entity) const DS_OVERRIDE DS_FINAL
    {
      return where_->apply_on(grid_view, entity);
    }

    virtual void apply_local(const EntityType& entity) DS_OVERRIDE DS_FINAL
    {
      wrapped_functor_.apply_local(entity);
    }

    virtual void finalize() DS_OVERRIDE DS_FINAL
    {
      wrapped_functor_.finalize();
    }

  private:
    Codim0FunctorType& wrapped_functor_;
    std::unique_ptr<const ApplyOn::WhichEntity<GridViewType>> where_;
  }; // class Codim0FunctorWrapper


  class Codim1Object : public Functor::Codim1<GridViewType>
  {
  public:
    ~Codim1Object()
    {
    }
    virtual bool apply_on(const GridViewType& grid_view, const IntersectionType& intersection) const = 0;
  };


  template <class Codim1FunctorType>
  class Codim1FunctorWrapper : public Codim1Object
  {
  public:
    Codim1FunctorWrapper(Codim1FunctorType& wrapped_functor, const ApplyOn::WhichIntersection<GridViewType>* where)
      : wrapped_functor_(wrapped_functor)
      , where_(where)
    {
    }

    virtual void prepare() DS_OVERRIDE DS_FINAL
    {
      wrapped_functor_.prepare();
    }

    virtual bool apply_on(const GridViewType& grid_view,
                          const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
    {
      return where_->apply_on(grid_view, intersection);
    }

    virtual void apply_local(const IntersectionType& intersection, const EntityType& inside_entity,
                             const EntityType& outside_entity) DS_OVERRIDE DS_FINAL
    {
      wrapped_functor_.apply_local(intersection, inside_entity, outside_entity);
    }

    virtual void finalize() DS_OVERRIDE DS_FINAL
    {
      wrapped_functor_.finalize();
    }

  private:
    Codim1FunctorType& wrapped_functor_;
    std::unique_ptr<const ApplyOn::WhichIntersection<GridViewType>> where_;
  }; // class Codim1FunctorWrapper


  class GridWalkerWrapper : public Codim0Object, public Codim1Object
  {
  public:
    GridWalkerWrapper(ThisType& grid_walker, const ApplyOn::WhichEntity<GridViewType>* which_entities)
      : grid_walker_(grid_walker)
      , which_entities_(which_entities)
      , which_intersections_(new ApplyOn::AllIntersections<GridViewType>())
    {
    }

    GridWalkerWrapper(ThisType& grid_walker, const ApplyOn::WhichIntersection<GridViewType>* which_intersections)
      : grid_walker_(grid_walker)
      , which_entities_(new ApplyOn::AllEntities<GridViewType>())
      , which_intersections_(which_intersections)
    {
    }

    virtual ~GridWalkerWrapper()
    {
    }

    virtual void prepare() DS_OVERRIDE DS_FINAL
    {
      grid_walker_.prepare();
    }

    virtual bool apply_on(const GridViewType& grid_view, const EntityType& entity) const DS_OVERRIDE DS_FINAL
    {
      return which_entities_->apply_on(grid_view, entity) && grid_walker_.apply_on(entity);
    }

    virtual bool apply_on(const GridViewType& grid_view,
                          const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
    {
      return which_intersections_->apply_on(grid_view, intersection) && grid_walker_.apply_on(intersection);
    }

    virtual void apply_local(const EntityType& entity) DS_OVERRIDE DS_FINAL
    {
      grid_walker_.apply_local(entity);
    }

    virtual void apply_local(const IntersectionType& intersection, const EntityType& inside_entity,
                             const EntityType& outside_entity) DS_OVERRIDE DS_FINAL
    {
      grid_walker_.apply_local(intersection, inside_entity, outside_entity);
    }

    virtual void finalize() DS_OVERRIDE DS_FINAL
    {
      grid_walker_.finalize();
    }

  private:
    ThisType& grid_walker_;
    std::unique_ptr<const ApplyOn::WhichEntity<GridViewType>> which_entities_;
    std::unique_ptr<const ApplyOn::WhichIntersection<GridViewType>> which_intersections_;
  }; // class GridWalkerWrapper


public:
  GridWalker(const GridViewType& grid_view)
    : grid_view_(grid_view)
  {
  }

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  void add(Functor::Codim0<GridViewType>& functor,
           const ApplyOn::WhichEntity<GridViewType>* where = new ApplyOn::AllEntities<GridViewType>())
  {
    codim0_functors_.emplace_back(new Codim0FunctorWrapper<Functor::Codim0<GridViewType>>(functor, where));
  }

  void add(Functor::Codim1<GridViewType>& functor,
           const ApplyOn::WhichIntersection<GridViewType>* where = new ApplyOn::AllIntersections<GridViewType>())
  {
    codim1_functors_.emplace_back(new Codim1FunctorWrapper<Functor::Codim1<GridViewType>>(functor, where));
  }

  void add(Functor::Codim0And1<GridViewType>& functor,
           const ApplyOn::WhichEntity<GridViewType>* which_entities = new ApplyOn::AllEntities<GridViewType>(),
           const ApplyOn::WhichIntersection<GridViewType>* which_intersections =
               new ApplyOn::AllIntersections<GridViewType>())
  {
    codim0_functors_.emplace_back(new Codim0FunctorWrapper<Functor::Codim0And1<GridViewType>>(functor, which_entities));
    codim1_functors_.emplace_back(
        new Codim1FunctorWrapper<Functor::Codim0And1<GridViewType>>(functor, which_intersections));
  }

  void add(Functor::Codim0And1<GridViewType>& functor,
           const ApplyOn::WhichIntersection<GridViewType>* which_intersections,
           const ApplyOn::WhichEntity<GridViewType>* which_entities = new ApplyOn::AllEntities<GridViewType>())
  {
    codim0_functors_.emplace_back(new Codim0FunctorWrapper<Functor::Codim0And1<GridViewType>>(functor, which_entities));
    codim1_functors_.emplace_back(
        new Codim1FunctorWrapper<Functor::Codim0And1<GridViewType>>(functor, which_intersections));
  }

  void add(ThisType& other,
           const ApplyOn::WhichEntity<GridViewType>* which_entities = new ApplyOn::AllEntities<GridViewType>(),
           const ApplyOn::WhichIntersection<GridViewType>* which_intersections =
               new ApplyOn::AllIntersections<GridViewType>())
  {
    if (&other == this)
      DUNE_THROW_COLORFULLY(Stuff::Exceptions::internal_error, "Do not add a GridWalker to itself!");
    codim0_functors_.emplace_back(new GridWalkerWrapper(other, which_entities));
    codim1_functors_.emplace_back(new GridWalkerWrapper(other, which_intersections));
  } // ... add(...)

  void add(ThisType& other, const ApplyOn::WhichIntersection<GridViewType>* which_intersections,
           const ApplyOn::WhichEntity<GridViewType>* which_entities = new ApplyOn::AllEntities<GridViewType>())
  {
    if (&other == this)
      DUNE_THROW_COLORFULLY(Stuff::Exceptions::internal_error, "Do not add a GridWalker to itself!");
    codim0_functors_.emplace_back(new GridWalkerWrapper(other, which_entities));
    codim1_functors_.emplace_back(new GridWalkerWrapper(other, which_intersections));
  } // ... add(...)

  void clear()
  {
    codim0_functors_ = std::vector<std::unique_ptr<Codim0Object>>();
    codim1_functors_ = std::vector<std::unique_ptr<Codim1Object>>();
  } // ... clear()

  virtual void prepare()
  {
    for (auto& functor : codim0_functors_)
      functor->prepare();
    for (auto& functor : codim1_functors_)
      functor->prepare();
  } // ... prepare()

  bool apply_on(const EntityType& entity) const
  {
    for (const auto& functor : codim0_functors_)
      if (functor->apply_on(grid_view_, entity))
        return true;
    return false;
  } // ... apply_on(...)

  bool apply_on(const IntersectionType& intersection) const
  {
    for (const auto& functor : codim1_functors_)
      if (functor->apply_on(grid_view_, intersection))
        return true;
    return false;
  } // ... apply_on(...)

  virtual void apply_local(const EntityType& entity)
  {
    for (auto& functor : codim0_functors_)
      if (functor->apply_on(grid_view_, entity))
        functor->apply_local(entity);
  } // ... apply_local(...)

  virtual void apply_local(const IntersectionType& intersection, const EntityType& inside_entity,
                           const EntityType& outside_entity)
  {
    for (auto& functor : codim1_functors_)
      if (functor->apply_on(grid_view_, intersection))
        functor->apply_local(intersection, inside_entity, outside_entity);
  } // ... apply_local(...)

  virtual void finalize()
  {
    for (auto& functor : codim0_functors_)
      functor->finalize();
    for (auto& functor : codim1_functors_)
      functor->finalize();
  } // ... finalize()

  void walk(const bool clear_stack = true)
  {
    // prepare functors
    prepare();

    // only do something, if we have to
    if ((codim0_functors_.size() + codim1_functors_.size()) > 0) {
      // walk the grid
      const auto entity_it_end = grid_view_.template end<0>();
      for (auto entity_it = grid_view_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
        const EntityType& entity = *entity_it;

        // apply codim0 functors
        apply_local(entity);

        // only walk the intersections, if there are codim1 functors present
        if (codim1_functors_.size() > 0) {
          // walk the intersections
          const auto intersection_it_end = grid_view_.iend(entity);
          for (auto intersection_it = grid_view_.ibegin(entity); intersection_it != intersection_it_end;
               ++intersection_it) {
            const auto& intersection = *intersection_it;

            // apply codim1 functors
            if (intersection.neighbor()) {
              const auto neighbor_ptr = intersection.outside();
              const auto& neighbor = *neighbor_ptr;
              apply_local(intersection, entity, neighbor);
            } else
              apply_local(intersection, entity, entity);

          } // walk the intersections
        } // only walk the intersections, if there are codim1 functors present
      } // walk the grid
    } // only do something, if we have to

    // finalize functors
    finalize();

    // clear the stack of functors
    if (clear_stack)
      clear();
  } // ... walk(...)

protected:
  const GridViewType& grid_view_;
  std::vector<std::unique_ptr<Codim0Object>> codim0_functors_;
  std::vector<std::unique_ptr<Codim1Object>> codim1_functors_;
}; // class GridWalker


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_GRIDWALKER_HH
