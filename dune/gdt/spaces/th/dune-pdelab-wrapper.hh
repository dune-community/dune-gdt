template <class GridViewImp, int polynomialOrder, size_t domainDim, class RangeFieldImp>
class TaylorHoodSpace
    : public DefaultProductSpace<DunePdelabCgSpaceWrapper<GridViewImp, polynomialOrder, RangeFieldImp, domainDim, 1>,
                                 DunePdelabCgSpaceWrapper<GridViewImp, polynomialOrder - 1, RangeFieldImp, 1, 1>>
{
public:
  typedef DunePdelabCgSpaceWrapper<GridViewImp, polynomialOrder, RangeFieldImp, domainDim, 1> VelocitySpaceType;
  typedef DunePdelabCgSpaceWrapper<GridViewImp, polynomialOrder - 1, RangeFieldImp, 1, 1> PressureSpaceType;
  typedef DefaultProductSpace<VelocitySpaceType, PressureSpaceType> BaseType;
  TaylorHoodSpace(const GridViewImp& grid_view)
    : BaseType(VelocitySpaceType(grid_view), PressureSpaceType(grid_view))
  {
  }
};
