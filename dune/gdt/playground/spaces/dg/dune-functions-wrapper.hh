#ifndef DUNE_GDT_PLAYGROUND_SPACES_DG_DUNE_FUNCTIONS_WRAPPER_HH
#define DUNE_GDT_PLAYGROUND_SPACES_DG_DUNE_FUNCTIONS_WRAPPER_HH

#include <dune/common/typetraits.hh>

#include <dune/functions/functionspacebases/lagrangedgbasis.hh>

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/playground/spaces/mapper/dune-functions-wrapper.hh>
#include <dune/gdt/playground/spaces/basefunctionset/dune-functions-wrapper.hh>

namespace Dune {
namespace GDT {

#if HAVE_DUNE_FUNCTIONS


// forward, to be used in the traits and to allow for specialization
template <class GL, int p, class R, size_t r, size_t rC = 1>
class DuneFunctionsDgSpaceWrapper
{
  static_assert(Dune::AlwaysFalse<GL>::value, "Untested for these dimensions!");
};


namespace internal {


template <class GL, int p, class R, size_t r, size_t rC>
class DuneFunctionsDgSpaceWrapperTraits
{
  static_assert(XT::Grid::is_view<GL>::value, "We Probably need to use TemporaryGridView from dune-xt-grid!");

public:
  typedef DuneFunctionsDgSpaceWrapper<GL, p, R, r, rC> derived_type;
  static const int polOrder = p;
  static const bool continuous = false;
  static const XT::Grid::Backends layer_backend = XT::Grid::extract_layer_backend<GL>::value;
  typedef Functions::LagrangeDGBasis<GL, p> BackendType;
  typedef DuneFunctionsMapperWrapper<GL, p, R, r, rC> MapperType;
  typedef DuneFunctionsBaseFunctionSetWrapper<GL, p, R, r, rC> BaseFunctionSetType;
  typedef double CommunicatorType;
  typedef GL GridLayerType;
  typedef R RangeFieldType;
}; // class DuneFunctionsDgSpaceWrapperTraits


} // namespace internal


template <class GL, int p, class R>
class DuneFunctionsDgSpaceWrapper<GL, p, R, 1, 1>
    : public SpaceInterface<internal::DuneFunctionsDgSpaceWrapperTraits<GL, p, R, 1, 1>, GL::dimension, 1, 1>
{
  typedef DuneFunctionsDgSpaceWrapper<GL, p, R, 1, 1> ThisType;
  typedef SpaceInterface<internal::DuneFunctionsDgSpaceWrapperTraits<GL, p, R, 1, 1>, GL::dimension, 1, 1> BaseType;

public:
  typedef internal::DuneFunctionsDgSpaceWrapperTraits<GL, p, R, 1, 1> Traits;

  using typename BaseType::BackendType;
  using typename BaseType::BaseFunctionSetType;
  using typename BaseType::CommunicatorType;
  using typename BaseType::EntityType;
  using typename BaseType::GridLayerType;
  using typename BaseType::MapperType;

  DuneFunctionsDgSpaceWrapper(GridLayerType grd_layr)
    : grid_layer_(new GridLayerType(grd_layr))
    , backend_(new BackendType(*grid_layer_))
    , mapper_(new MapperType(backend_))
    , communicator_(new CommunicatorType(0.))
  {
  }

  DuneFunctionsDgSpaceWrapper(const ThisType& other) = default;
  DuneFunctionsDgSpaceWrapper(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  const GridLayerType& grid_layer() const
  {
    return *grid_layer_;
  }

  GridLayerType& grid_layer()
  {
    return *grid_layer_;
  }

  const BackendType& backend() const
  {
    return *backend_;
  }

  const MapperType& mapper() const
  {
    return *mapper_;
  }

  BaseFunctionSetType base_function_set(const EntityType& entity) const
  {
    return BaseFunctionSetType(backend_, entity);
  }

  CommunicatorType& communicator() const
  {
    return *communicator_;
  }

private:
  std::shared_ptr<GridLayerType> grid_layer_;
  const std::shared_ptr<const BackendType> backend_;
  const std::shared_ptr<const MapperType> mapper_;
  mutable std::shared_ptr<CommunicatorType> communicator_;
}; // class DuneFunctionsDgSpaceWrapper


#else // HAVE_DUNE_FUNCTIONS


template <class GL, int p, class R, size_t r, size_t rC = 1>
class DuneFunctionsDgSpaceWrapper
{
  static_assert(Dune::AlwaysFalse<GL>::value, "You are missing dune-functions!");
};


#endif // HAVE_DUNE_FUNCTIONS

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_SPACES_DG_DUNE_FUNCTIONS_WRAPPER_HH
