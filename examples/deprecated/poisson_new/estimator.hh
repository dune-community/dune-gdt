#ifndef ESTIMATOR_HH
#define ESTIMATOR_HH

//- Dune-fem includes
#include <dune/fem/quadrature/caching/twistutility.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>
#include <dune/fem/operator/matrix/blockmatrix.hh>
#include <dune/fem/space/dgspace.hh>


// Estimator
// ---------
template <class DiscreteFunction>
class Estimator
{
  typedef Estimator<DiscreteFunction> ThisType;

public:
  typedef DiscreteFunction DiscreteFunctionType;

  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
  typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;

  typedef typename GridPartType::GridType GridType;
  typedef typename GridPartType::IndexSetType IndexSetType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;

  typedef typename IntersectionIteratorType::Intersection IntersectionType;

  typedef typename GridType::template Codim<0>::Entity ElementType;
  typedef typename GridType::template Codim<0>::EntityPointer ElementPointerType;
  typedef typename ElementType::Geometry GeometryType;
  static const int dimension = GridType::dimension;

  typedef Dune::CachingQuadrature<GridPartType, 0> ElementQuadratureType;
  typedef Dune::CachingQuadrature<GridPartType, 1> FaceQuadratureType;
  typedef Dune::FieldMatrix<double, dimension, dimension> JacobianInverseType;
  typedef std::vector<double> ErrorIndicatorType;

private:
  const DiscreteFunctionType& uh_;
  const DiscreteFunctionSpaceType& dfSpace_;
  GridPartType& gridPart_;
  const IndexSetType& indexSet_;
  GridType& grid_;
  ErrorIndicatorType indicator_;

public:
  explicit Estimator(const DiscreteFunctionType& uh)
    : uh_(uh)
    , dfSpace_(uh.space())
    , gridPart_(dfSpace_.gridPart())
    , indexSet_(gridPart_.indexSet())
    , grid_(gridPart_.grid())
    , indicator_(indexSet_.size(0))
  {
  }

  void clear()
  {
    typedef typename ErrorIndicatorType::iterator IteratorType;
    const IteratorType end = indicator_.end();
    for (IteratorType it = indicator_.begin(); it != end; ++it)
      *it = 0.0;
  }

  template <class RHSFunctionType>
  double estimate(const RHSFunctionType& rhs)
  {
    clear();
    const IteratorType end = dfSpace_.end();
    for (IteratorType it = dfSpace_.begin(); it != end; ++it)
      estimateLocal(rhs, *it);

    double error                                                   = 0.0;
    const typename ErrorIndicatorType::const_iterator endEstimator = indicator_.end();
    for (typename ErrorIndicatorType::const_iterator it = indicator_.begin(); it != endEstimator; ++it)
      error += *it;
    std::cout << "Estimated H1-Error: " << sqrt(error) << std::endl;
    return sqrt(error);
  }

  template <class RHSFunctionType>
  bool estimateAndMark(const RHSFunctionType& rhs, const double tolerance)
  {
    const double error = estimate(rhs);
    return (error < tolerance ? false : mark(0.9 * tolerance));
  }

  //! calculates the error estimator and also marks the
  //! grid for adaptation
  template <class RHSFunctionType>
  static bool estimateAndMark(const DiscreteFunctionType& uh, const RHSFunctionType& rhs, const double tolerance)
  {
    ThisType estimator(uh);
    return estimator.estimateAndMark(rhs, tolerance);
  }

  //! mark all elements due to given tolerance
  bool mark(const double tolerance) const
  {
    // get local tolerance
    const double localTol2 = tolerance * tolerance / (double)indexSet_.size(0);

    bool marked = false;
    // loop over all elements
    const IteratorType end = dfSpace_.end();
    for (IteratorType it = dfSpace_.begin(); it != end; ++it) {
      const ElementType& entity = *it;
      // check local error indicator
      if (indicator_[indexSet_.index(entity)] > localTol2) {
// mark entity for refinement
#ifdef OLD_DUNE_GRID_VERSION
        grid_.mark(1, it);
#else
        grid_.mark(1, entity);
#endif
        marked = true;
      }
    }
    return marked;
  }

private:
  //! caclulate error on element
  template <class RHSFunctionType>
  void estimateLocal(const RHSFunctionType& rhs, const ElementType& entity)
  {
    const typename ElementType::Geometry& geometry = entity.geometry();

    const double volume            = geometry.volume();
    const double h2                = (dimension == 2 ? volume : std::pow(volume, 2.0 / (double)dimension));
    const int index                = indexSet_.index(entity);
    const LocalFunctionType uLocal = uh_.localFunction(entity);

    ElementQuadratureType quad(entity, 2 * (dfSpace_.order() + 1));
    const int numQuadraturePoints = quad.nop();
    for (int qp = 0; qp < numQuadraturePoints; ++qp) {
      const typename ElementQuadratureType::CoordinateType& x = quad.point(qp);
      const double weight                                     = quad.weight(qp) * geometry.integrationElement(x);

      RangeType y;
      rhs.evaluate(geometry.global(x), y);

      RangeType tmp;
      laplacianLocal(uLocal, quad[qp], tmp);
      y += tmp;

      indicator_[index] += h2 * weight * (y * y);
    }

    IntersectionIteratorType end = gridPart_.iend(entity);
    for (IntersectionIteratorType it = gridPart_.ibegin(entity); it != end; ++it) {
      const IntersectionType& intersection = *it;
      if (intersection.neighbor())
        estimateIntersection(intersection, entity, uLocal);
    }
    // std::cout<<"ELT:"<<index<<"Error="<<indicator_[index]<<"\n";
  }

  //! caclulate error on boundary intersections
  void estimateIntersection(const IntersectionType& intersection, const ElementType& inside,
                            const LocalFunctionType& uInside)
  {
    const ElementPointerType pOutside = intersection.outside();
    const ElementType& outside        = *pOutside;

    const int quadOrder = 2 * (dfSpace_.order() - 1);

    const int insideIndex  = indexSet_.index(inside);
    const int outsideIndex = indexSet_.index(outside);

    const bool isOutsideInterior = (outside.partitionType() == Dune::InteriorEntity);
    if (!isOutsideInterior || (insideIndex < outsideIndex)) {
      const LocalFunctionType uOutside = uh_.localFunction(outside);

      double error;

#ifdef OLD_DUNE_GRID_VERSION

      typedef TwistUtility<GridType> TwistUtilityType;
      if (TwistUtilityType::conforming(gridPart_.grid(), intersection))
#else
      if (intersection.conforming())
#endif
      {
        FaceQuadratureType quadInside(gridPart_, intersection, quadOrder, FaceQuadratureType::INSIDE);
        FaceQuadratureType quadOutside(gridPart_, intersection, quadOrder, FaceQuadratureType::OUTSIDE);
        error = estimateIntersection(intersection, quadInside, quadOutside, uInside, uOutside);
      } else {
        typedef typename FaceQuadratureType::NonConformingQuadratureType NonConformingQuadratureType;
        NonConformingQuadratureType quadInside(gridPart_, intersection, quadOrder, NonConformingQuadratureType::INSIDE);
        NonConformingQuadratureType quadOutside(
            gridPart_, intersection, quadOrder, NonConformingQuadratureType::OUTSIDE);
        error = estimateIntersection(intersection, quadInside, quadOutside, uInside, uOutside);
      }

      if (error > 0.0) {
        const double volume = 0.5 * (inside.geometry().volume() + outside.geometry().volume());
        const double h      = (dimension == 1 ? volume : std::pow(volume, 1.0 / (double)dimension));
        indicator_[insideIndex] += h * error;
        if (isOutsideInterior)
          indicator_[outsideIndex] += h * error;
      }
    }
  }

  //! caclulate error on element intersections
  template <class Quadrature>
  double estimateIntersection(const IntersectionType& intersection, const Quadrature& quadInside,
                              const Quadrature& quadOutside, const LocalFunctionType& uInside,
                              const LocalFunctionType& uOutside) const
  {
    double error                  = 0.0;
    const int numQuadraturePoints = quadInside.nop();
    for (int qp = 0; qp < numQuadraturePoints; ++qp) {
      DomainType integrationNormal    = intersection.integrationOuterNormal(quadInside.localPoint(qp));
      const double integrationElement = integrationNormal.two_norm();

      // evaluate | (d u_l * n_l) + (d u_r * n_r) | = | (d u_l - d u_r) * n_l |
      JacobianRangeType jacobianInside;
      JacobianRangeType jacobianOutside;
      uInside.jacobian(quadInside[qp], jacobianInside);
      uOutside.jacobian(quadOutside[qp], jacobianOutside);

      RangeType jump;
      jacobianInside -= jacobianOutside;
      jacobianInside.mv(integrationNormal, jump);

      error += quadInside.weight(qp) * (jump * jump) / integrationElement;
    }
    return error;
  }

  template <class PointType>
  void laplacianLocal(const LocalFunctionType& u_h, const PointType& x, RangeType& result) const
  {
    typename LocalFunctionType::HessianRangeType hessian;
    u_h.hessian(x, hessian);

    result = 0;
    for (int r = 0; r < LocalFunctionType::dimRange; ++r) {
      for (int i = 0; i < dimension; ++i)
        result[r] += hessian[r][i][i];
    }
  }
};


#endif // #ifndef ESTIMATOR_HH
