// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_LEBEDEVQUADRATURE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_LEBEDEVQUADRATURE_HH

#include <vector>

#include <boost/geometry.hpp>

#include <dune/common/fvector.hh>

#include <dune/geometry/quadraturerules.hh>

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


// Converts from (x, y, z) to (theta, phi) on the unit sphere s.t.
// (x, y, z) = (sin(theta) cos(phi), sin(theta) sin(phi), cos(theta)).
// with 0 \leq \theta \leq \pi and 0 \leq \varphi \leq 2\pi
// Note: Boost uses a different naming convention in its documentation (theta <--> phi)
template <class DomainFieldType>
struct CoordinateConverter
{
  typedef FieldVector<DomainFieldType, 3> CartesianCoordType;
  typedef FieldVector<DomainFieldType, 2> SphericalCoordType;
  typedef typename boost::geometry::model::point<DomainFieldType, 3, typename boost::geometry::cs::cartesian>
      BoostCartesianCoordType;
  typedef typename boost::geometry::model::point<DomainFieldType,
                                                 2,
                                                 typename boost::geometry::cs::spherical<boost::geometry::radian>>
      BoostSphericalCoordType;

  static SphericalCoordType to_spherical(const CartesianCoordType& x)
  {
    BoostCartesianCoordType x_boost(x[0], x[1], x[2]);
    BoostSphericalCoordType x_spherical_boost;
    boost::geometry::transform(x_boost, x_spherical_boost);
    return SphericalCoordType{boost::geometry::get<1>(x_spherical_boost), boost::geometry::get<0>(x_spherical_boost)};
  }

  static CartesianCoordType to_cartesian(const SphericalCoordType& x_spherical)
  {
    BoostSphericalCoordType x_spherical_boost(x_spherical[1], x_spherical[0]);
    BoostCartesianCoordType x_boost;
    boost::geometry::transform(x_spherical_boost, x_boost);
    return CartesianCoordType{
        boost::geometry::get<0>(x_boost), boost::geometry::get<1>(x_boost), boost::geometry::get<2>(x_boost)};
  }
};

template <class FieldType, bool cartesian = false>
struct LebedevQuadratureCreator
{
  static Dune::QuadratureRule<FieldType, 2> create(std::vector<FieldVector<FieldType, 2>>& positions,
                                                   std::vector<FieldType>& weights)
  {
    Dune::QuadratureRule<FieldType, 2> ret;
    for (size_t ii = 0; ii < weights.size(); ++ii)
      ret.emplace_back(Dune::QuadraturePoint<double, 2>(positions[ii], weights[ii]));
    return ret;
  }
};

template <class FieldType>
struct LebedevQuadratureCreator<FieldType, true>
{
  static Dune::QuadratureRule<FieldType, 3> create(std::vector<FieldVector<FieldType, 2>>& positions,
                                                   std::vector<FieldType>& weights)
  {
    Dune::QuadratureRule<FieldType, 3> ret;
    for (size_t ii = 0; ii < weights.size(); ++ii)
      ret.emplace_back(
          Dune::QuadraturePoint<double, 3>(CoordinateConverter<FieldType>::to_cartesian(positions[ii]), weights[ii]));
    return ret;
  }
};

template <class FieldType, bool cartesian = false>
class LebedevQuadrature
{
public:
  static Dune::QuadratureRule<FieldType, 2 + cartesian> get(const size_t /*requested_order*/);

private:
  static void get_positions_and_weights(const size_t /*order*/,
                                        std::vector<FieldVector<FieldType, 2>>& /*positions*/,
                                        std::vector<FieldType>& /*weights*/);

  static const std::vector<size_t> allowed_orders;
}; // class LebedevQuadrature


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_LEBEDEVQUADRATURE_HH
