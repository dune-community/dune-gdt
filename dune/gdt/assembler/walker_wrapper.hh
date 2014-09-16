#ifndef DUNE_GDT_ASSEMBLER_WALKER_WRAPPER_HH
#define DUNE_GDT_ASSEMBLER_WALKER_WRAPPER_HH

#include <dune/gdt/assembler/apply_on.hh>
#include <dune/gdt/assembler/functors.hh>

namespace Dune {
namespace GDT {

template <class GridViewType>
class Codim0Object
  : public Functor::Codim0< GridViewType >
{
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
public:
  ~ Codim0Object() {}
  virtual bool apply_on(const GridViewType& grid_view, const EntityType& entity) const = 0;
};

template<class GridViewType, class Codim0FunctorType>
class Codim0FunctorWrapper
  : public Codim0Object<GridViewType>
{
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
public:
  Codim0FunctorWrapper(Codim0FunctorType& wrapped_functor,
                       const ApplyOn::WhichEntity< GridViewType >* where)
    : wrapped_functor_(wrapped_functor)
    , where_(where)
  {}

  virtual ~Codim0FunctorWrapper() {}

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
  std::unique_ptr< const ApplyOn::WhichEntity< GridViewType > > where_;
}; // class Codim0FunctorWrapper

template <class GridViewType>
class Codim1Object
  : public Functor::Codim1< GridViewType >
{
  typedef typename GridViewType::Intersection IntersectionType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
public:
  ~ Codim1Object() {}
  virtual bool apply_on(const GridViewType& grid_view, const IntersectionType& intersection) const = 0;
};


template<class GridViewType, class Codim1FunctorType>
class Codim1FunctorWrapper
  : public Codim1Object<GridViewType>
{
  typedef typename GridViewType::Intersection IntersectionType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
public:
  Codim1FunctorWrapper(Codim1FunctorType& wrapped_functor,
                       const ApplyOn::WhichIntersection< GridViewType >* where)
    : wrapped_functor_(wrapped_functor)
    , where_(where)
  {}

  virtual void prepare() DS_OVERRIDE DS_FINAL
  {
    wrapped_functor_.prepare();
  }

  virtual bool apply_on(const GridViewType& grid_view, const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
  {
    return where_->apply_on(grid_view, intersection);
  }

  virtual void apply_local(const IntersectionType& intersection,
                           const EntityType& inside_entity,
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
  std::unique_ptr< const ApplyOn::WhichIntersection< GridViewType > > where_;
}; // class Codim1FunctorWrapper

template<class GridViewType, class WalkerType>
class GridWalkerWrapper
  : public Codim0Object<GridViewType>
  , public Codim1Object<GridViewType>
{
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::Intersection IntersectionType;
public:
  GridWalkerWrapper(WalkerType& grid_walker,
                    const ApplyOn::WhichEntity< GridViewType >* which_entities)
    : grid_walker_(grid_walker)
    , which_entities_(which_entities)
    , which_intersections_(new ApplyOn::AllIntersections< GridViewType >())
  {}

  GridWalkerWrapper(WalkerType& grid_walker,
                    const ApplyOn::WhichIntersection< GridViewType >* which_intersections)
    : grid_walker_(grid_walker)
    , which_entities_(new ApplyOn::AllEntities< GridViewType >())
    , which_intersections_(which_intersections)
  {}

  virtual ~GridWalkerWrapper() {}

  virtual void prepare() DS_OVERRIDE DS_FINAL
  {
    grid_walker_.prepare();
  }

  virtual bool apply_on(const GridViewType& grid_view, const EntityType& entity) const DS_OVERRIDE DS_FINAL
  {
    return which_entities_->apply_on(grid_view, entity) && grid_walker_.apply_on(entity);
  }

  virtual bool apply_on(const GridViewType& grid_view, const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
  {
    return which_intersections_->apply_on(grid_view, intersection) && grid_walker_.apply_on(intersection);
  }

  virtual void apply_local(const EntityType& entity) DS_OVERRIDE DS_FINAL
  {
    grid_walker_.apply_local(entity);
  }

  virtual void apply_local(const IntersectionType& intersection,
                           const EntityType& inside_entity,
                           const EntityType& outside_entity) DS_OVERRIDE DS_FINAL
  {
    grid_walker_.apply_local(intersection, inside_entity, outside_entity);
  }

  virtual void finalize() DS_OVERRIDE DS_FINAL
  {
    grid_walker_.finalize();
  }

private:
  WalkerType& grid_walker_;
  std::unique_ptr< const ApplyOn::WhichEntity< GridViewType > > which_entities_;
  std::unique_ptr< const ApplyOn::WhichIntersection< GridViewType > > which_intersections_;
}; // class GridWalkerWrapper

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_WALKER_WRAPPER_HH
