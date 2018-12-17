// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_PRODUCTQUADRATURE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_PRODUCTQUADRATURE_HH

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/coordinates.hh>
#include <dune/xt/common/string.hh>

#include <dune/xt/grid/gridprovider.hh>

namespace Dune {
namespace GDT {


// use product quadrature for phi and theta to get quadrature on sphere
template <class FieldType>
Dune::QuadratureRule<FieldType, 3> product_quadrature(const size_t theta_grid_size = 100,
                                                      const size_t phi_grid_size = 100)
{
  // create grids for phi and theta
  typedef typename Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<FieldType, 1>> GridType;
  Dune::XT::Common::Configuration theta_grid_config;
  theta_grid_config["type"] = Dune::XT::Grid::cube_gridprovider_default_config()["type"];
  theta_grid_config["lower_left"] = "[0]";
  theta_grid_config["upper_right"] = "[" + Dune::XT::Common::to_string(M_PI, 15) + "]";
  theta_grid_config["num_elements"] = XT::Common::to_string(theta_grid_size);
  const auto theta_grid_ptr = Dune::XT::Grid::CubeGridProviderFactory<GridType>::create(theta_grid_config).grid_ptr();
  Dune::XT::Common::Configuration phi_grid_config;
  phi_grid_config["type"] = Dune::XT::Grid::cube_gridprovider_default_config()["type"];
  phi_grid_config["lower_left"] = "[0]";
  phi_grid_config["upper_right"] = "[" + Dune::XT::Common::to_string(2. * M_PI, 15) + "]";
  phi_grid_config["num_elements"] = phi_grid_size;
  const auto phi_grid_ptr = Dune::XT::Grid::CubeGridProviderFactory<GridType>::create(phi_grid_config).grid_ptr();
  Dune::QuadratureRule<FieldType, 2> product_quad_spherical;
  const auto phi_grid_view = phi_grid_ptr->leafGridView();
  const auto theta_grid_view = theta_grid_ptr->leafGridView();
  for (const auto& phi_entity : Dune::elements(phi_grid_view)) {
    const auto phi_quad_rule = Dune::QuadratureRules<FieldType, 1>::rule(phi_entity.geometry().type(), 60);
    for (const auto& theta_entity : Dune::elements(theta_grid_view)) {
      const auto theta_quad_rule = Dune::QuadratureRules<FieldType, 1>::rule(theta_entity.geometry().type(), 60);
      for (const auto& phi_quad_point : phi_quad_rule)
        for (const auto& theta_quad_point : theta_quad_rule)
          product_quad_spherical.push_back(Dune::QuadraturePoint<FieldType, 2>(
              {theta_entity.geometry().global(theta_quad_point.position()),
               phi_entity.geometry().global(phi_quad_point.position())},
              std::sin(theta_entity.geometry().global(theta_quad_point.position())[0]) * phi_quad_point.weight()
                  * theta_quad_point.weight() * phi_entity.geometry().volume() * theta_entity.geometry().volume()));
    }
  }
  Dune::QuadratureRule<FieldType, 3> product_quad_cartesian;
  for (const auto quad_point : product_quad_spherical)
    product_quad_cartesian.push_back(Dune::QuadraturePoint<FieldType, 3>(
        Dune::XT::Common::CoordinateConverter<FieldType>::to_cartesian(quad_point.position()), quad_point.weight()));
  return product_quad_cartesian;
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_PRODUCTQUADRATURE_HH
