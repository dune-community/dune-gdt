#ifndef DUNE_GDT_ASSEMBLER_FUNCTORS_HH
#define DUNE_GDT_ASSEMBLER_FUNCTORS_HH

namespace Dune {
namespace GDT {
namespace Functor {

template <class GridViewImp>
class Codim0
{
public:
  typedef GridViewImp GridViewType;
  typedef typename Stuff::Grid::Entity<GridViewType>::Type EntityType;

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
  typedef typename Stuff::Grid::Entity<GridViewType>::Type EntityType;
  typedef typename Stuff::Grid::Intersection<GridViewType>::Type IntersectionType;

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
  typedef typename Stuff::Grid::Entity<GridViewType>::Type EntityType;
  typedef typename Stuff::Grid::Intersection<GridViewType>::Type IntersectionType;

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
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_FUNCTORS_HH
