// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_SAP_HH
#define DUNE_GDT_SAP_HH

#include <set>
#include <string>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/grid/boundaryinfo/interfaces.hh>
#include <dune/xt/grid/type_traits.hh>

namespace Dune {
namespace GDT {


template <size_t d, class R = double>
class SelfacoustophoreticParticleTools
{
  static_assert(d == 1 || d == 2, "Not yet implemented for 3d!");
  using ThisType = SelfacoustophoreticParticleTools<d, R>;

public:
  static const constexpr size_t m = d + 1;

  SelfacoustophoreticParticleTools(const R& p_0 = 101325,
                                   const R& rho_0 = 998.2,
                                   const R& c_f = 1482,
                                   const R& nu_s = 1.002e-3,
                                   const R& nu_b = 2.87e-3)
    : p_0_(p_0)
    , rho_0_(rho_0)
    , c_f_(c_f)
    , nu_s_(nu_s)
    , nu_b_(nu_b)
  {
  }

  XT::Common::FieldVector<R, m> to_primitive(const FieldVector<R, m>& conservative_variables) const
  {
    XT::Common::FieldVector<R, m> primitive_variables;
    const auto& rho = conservative_variables[0];
    primitive_variables[0] = rho;
    for (size_t ii = 0; ii < d; ++ii)
      primitive_variables[1 + ii] = conservative_variables[1 + ii] / rho;
    return primitive_variables;
  }

  XT::Common::FieldVector<R, m> to_conservative(const FieldVector<R, m>& primitive_variables) const
  {
    XT::Common::FieldVector<R, m> conservative_variables;
    const auto& rho = primitive_variables[0];
    conservative_variables[0] = rho;
    for (size_t ii = 0; ii < d; ++ii)
      conservative_variables[1 + ii] = primitive_variables[1 + ii] * rho;
    return conservative_variables;
  }

  XT::Common::FieldVector<R, 1> pressure_from_density(const R& density) const
  {
    return p_0_ / (rho_0_ * c_f_ * c_f_) + (density - 1);
  }

  XT::Common::FieldVector<R, d> velocity_from_conservative(const FieldVector<R, m>& conservative_variables) const
  {
    const auto& rho = conservative_variables[0];
    XT::Common::FieldVector<R, d> v;
    for (size_t ii = 0; ii < d; ++ii)
      v[ii] = conservative_variables[1 + ii] / rho;
    return v;
  }

  XT::Common::FieldMatrix<R, d, m> flux(const FieldVector<R, m>& conservative_variables) const
  {
    const auto& rho = conservative_variables[0];
    const auto p = pressure_from_density(rho);
    const auto v = velocity_from_conservative(conservative_variables);
    XT::Common::FieldMatrix<R, d, m> ret;
    for (size_t ss = 0; ss < d; ++ss) {
      auto& f_s = ret[ss];
      f_s[0] = rho * v[ss];
      for (size_t ii = 0; ii < d; ++ii)
        f_s[1 + ii] = rho * v[ii] * v[ss] + (ss == ii ? 1 : 0) * p;
    }
    return ret;
  } // ... flux(...)

  template <class D>
  XT::Common::FieldVector<R, m> flux_at_impermeable_walls(const FieldVector<R, m>& conservative_variables,
                                                          const FieldVector<D, d>& normal) const
  {
    const auto pressure = pressure_from_density(conservative_variables[0])[0];
    const auto tmp = normal * pressure;
    XT::Common::FieldVector<R, m> ret(0.);
    for (size_t ii = 0; ii < d; ++ii)
      ret[ii + 1] = tmp[ii];
    return ret;
  } // ... flux_at_impermeable_walls(...)

  XT::Common::FieldVector<XT::Common::FieldMatrix<R, m, m>, d>
  flux_jacobian(const FieldVector<R, m>& conservative_variables) const
  {
    return dim_switch<>::flux_jacobian(*this, conservative_variables);
  }

  XT::Common::FieldMatrix<R, m, m> flux_jacobi_matrix(const FieldVector<R, m>& conservative_variables,
                                                      const FieldVector<double, d>& normal) const
  {
    return dim_switch<>::flux_jacobi_matrix(*this, conservative_variables, normal);
  }

  XT::Common::FieldVector<R, m> eigenvalues_flux_jacobi_matrix(const FieldVector<R, m>& conservative_variables,
                                                               const FieldVector<double, d>& normal) const
  {
    return dim_switch<>::eigenvalues_flux_jacobi_matrix(*this, conservative_variables, normal);
  }

  XT::Common::FieldMatrix<R, m, m> eigenvaluematrix_flux_jacobi_matrix(const FieldVector<R, m>& conservative_variables,
                                                                       const FieldVector<double, d>& normal) const
  {
    const auto evs = eigenvalues_flux_jacobi_matrix(conservative_variables, normal);
    XT::Common::FieldMatrix<R, m, m> ret(0.);
    for (size_t ii = 0; ii < m; ++ii)
      ret[ii][ii] = evs[ii];
    return ret;
  }

  XT::Common::FieldMatrix<R, m, m> eigenvectors_flux_jacobi_matrix(const FieldVector<R, m>& conservative_variables,
                                                                   const FieldVector<double, d>& normal) const
  {
    return dim_switch<>::eigenvectors_flux_jacobi_matrix(*this, conservative_variables, normal);
  }

  XT::Common::FieldMatrix<R, m, m> eigenvectors_inv_flux_jacobi_matrix(const FieldVector<R, m>& conservative_variables,
                                                                       const FieldVector<double, d>& normal) const
  {
    return dim_switch<>::eigenvectors_inv_flux_jacobi_matrix(*this, conservative_variables, normal);
  }

  template <class E, class D, class GL>
  void visualize(const XT::Functions::LocalizableFunctionInterface<E, D, d, R, m>& u_conservative,
                 const GL& grid_layer,
                 const std::string filename_prefix = "",
                 const std::string filename_suffix = "",
                 const bool subsampling = false,
                 const XT::Common::Parameter& param = {})
  {
    const std::string prefix = filename_prefix.empty() ? "" : filename_prefix + "_";
    const std::string suffix = filename_suffix.empty() ? "" : "_" + filename_suffix;
    const auto density = XT::Functions::make_sliced_function<1>(u_conservative, {0}, "density");
    density.visualize(grid_layer, prefix + "density" + suffix, subsampling, VTK::appendedraw, param);
    XT::Functions::make_transformed_function<1, 1, R>(
        density, [&](const auto& rho) { return rho - 1.; }, "density_variation")
        .visualize(grid_layer, prefix + "density_variation" + suffix, subsampling, VTK::appendedraw, param);
    XT::Functions::make_sliced_function<d>(u_conservative, dim_switch<>::velocity_indices(), "density_times_velocity")
        .visualize(grid_layer, prefix + "density_times_velocity" + suffix, subsampling, VTK::appendedraw, param);
    const auto u_primitive = XT::Functions::make_transformed_function<m, 1, R>(
        u_conservative, [&](const auto& cons) { return to_primitive(cons); });
    XT::Functions::make_sliced_function<d>(u_primitive, dim_switch<>::velocity_indices(), "velocity")
        .visualize(grid_layer, prefix + "velocity" + suffix, subsampling, VTK::appendedraw, param);
    const auto pressure = XT::Functions::make_transformed_function<1, 1, R>(
        density, [&](const auto& rho) { return pressure_from_density(rho); }, "pressure");
    pressure.visualize(grid_layer, prefix + "pressure" + suffix, subsampling, VTK::appendedraw, param);
    //    XT::Functions::make_transformed_function<1, 1, R>(
    //        pressure, [&](const auto& p) { return p - p_0_; }, "pressure_variation")
    //        .visualize(grid_layer, prefix + "pressure_variation" + suffix, subsampling, VTK::appendedraw, param);
  } // ... visualize(...)

private:
  template <size_t d_ = d, typename anything = void>
  struct dim_switch;

  template <typename anything>
  struct dim_switch<1, anything>
  {
    static std::array<size_t, d> velocity_indices()
    {
      return {1};
    }

    static XT::Common::FieldVector<XT::Common::FieldMatrix<R, m, m>, d>
    flux_jacobian(const ThisType& self, const FieldVector<R, m>& conservative_variables)
    {
      const auto v = self.velocity_from_conservative(conservative_variables);
      // clang-format off
      XT::Common::FieldMatrix<R, m, m> jacobian = {{        0,     1},
                                                   {1 - v * v, 2 * v}}; // clang-format on
      return {jacobian};
    }

    static XT::Common::FieldMatrix<R, m, m> flux_jacobi_matrix(const ThisType& self,
                                                               const FieldVector<R, m>& conservative_variables,
                                                               const FieldVector<double, d>& normal)
    {
      auto A = flux_jacobian(self, conservative_variables);
      A[0] *= normal[0];
      return A;
    }

    static XT::Common::FieldVector<R, m> eigenvalues_flux_jacobi_matrix(const ThisType& self,
                                                                        const FieldVector<R, m>& conservative_variables,
                                                                        const FieldVector<double, d>& normal)
    {
      const auto v = self.velocity_from_conservative(conservative_variables);
      const auto& n = normal[0];
      return {(v - 1) * n, (v + 1) * n};
    }

    static XT::Common::FieldMatrix<R, m, m> eigenvectors_flux_jacobi_matrix(
        const ThisType& self, const FieldVector<R, m>& conservative_variables, const FieldVector<double, d>& /*normal*/)
    {
      const auto v = self.velocity_from_conservative(conservative_variables);
      // clang-format off
      return {{  1. ,   1. },
              {v - 1, v + 1}}; // clang-format on
    }

    static XT::Common::FieldMatrix<R, m, m> eigenvectors_inv_flux_jacobi_matrix(
        const ThisType& self, const FieldVector<R, m>& conservative_variables, const FieldVector<double, d>& /*normal*/)
    {
      const auto v = self.velocity_from_conservative(conservative_variables);
      // clang-format off
      return {{(1 + v) / 2, -1. / 2},
              {(1 - v) / 2,  1. / 2}}; // clang-format on
    }
  }; // struct dim_switch<1, ...>

  template <typename anything>
  struct dim_switch<2, anything>
  {
    static std::array<size_t, d> velocity_indices()
    {
      return {1, 2};
    }

    static XT::Common::FieldVector<XT::Common::FieldMatrix<R, m, m>, d>
    flux_jacobian(const ThisType& self, const FieldVector<R, m>& conservative_variables)
    {
      static_assert(AlwaysFalse<R>::value, "");
    }

    static XT::Common::FieldMatrix<R, m, m> flux_jacobi_matrix(const ThisType& self,
                                                               const FieldVector<R, m>& conservative_variables,
                                                               const FieldVector<double, d>& normal)
    {
      static_assert(AlwaysFalse<R>::value, "");
    }

    static XT::Common::FieldVector<R, m> eigenvalues_flux_jacobi_matrix(const ThisType& self,
                                                                        const FieldVector<R, m>& conservative_variables,
                                                                        const FieldVector<double, d>& normal)
    {
      static_assert(AlwaysFalse<R>::value, "");
    }

    static XT::Common::FieldMatrix<R, m, m> eigenvectors_flux_jacobi_matrix(
        const ThisType& self, const FieldVector<R, m>& conservative_variables, const FieldVector<double, d>& normal)
    {
      static_assert(AlwaysFalse<R>::value, "");
    }

    static XT::Common::FieldMatrix<R, m, m> eigenvectors_inv_flux_jacobi_matrix(
        const ThisType& self, const FieldVector<R, m>& conservative_variables, const FieldVector<double, d>& normal)
    {
      static_assert(AlwaysFalse<R>::value, "");
    }
  }; //  struct dim_switch<2, ...>

  template <typename anything>
  struct dim_switch<3, anything>
  {
    static std::array<size_t, d> velocity_indices()
    {
      return {1, 2, 3};
    }
  }; //  struct dim_switch<3, ...>

  const R p_0_;
  const R rho_0_;
  const R c_f_;
  const R nu_s_;
  const R nu_b_;
}; // class SelfacoustophoreticParticleTools


template <class GridLayerType, class IndexSetType>
void visualize_grid(const GridLayerType& grid_layer,
                    const IndexSetType& index_set,
                    XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info,
                    const std::string& filename)
{
  const auto num_elements = index_set.size(0);
  // walk the grid once to find all boundary types
  std::set<std::string> boundary_types;
  for (auto&& element : elements(grid_layer))
    for (auto&& intersection : intersections(grid_layer, element)) {
      if (intersection.boundary() && intersection.neighbor())
        boundary_types.insert("periodic_boundary");
      boundary_types.insert(boundary_info.type(intersection).id());
    }
  // prepare data structures
  std::map<std::string, std::vector<double>> name_to_data_map;
  const auto no_boundary_type = XT::Grid::NoBoundary().id();
  for (const auto& boundary_type : boundary_types)
    if (boundary_type != no_boundary_type)
      name_to_data_map[boundary_type] = std::vector<double>(num_elements, 0.);
  // walk the grid again to fill data structures
  for (auto&& element : elements(grid_layer)) {
    const auto element_index = index_set.index(element);
    for (auto&& intersection : intersections(grid_layer, element)) {
      if (intersection.boundary() && intersection.neighbor())
        name_to_data_map["periodic_boundary"][element_index] = 1.;
      const auto boundary_type = boundary_info.type(intersection).id();
      if (boundary_type != no_boundary_type)
        name_to_data_map[boundary_type][element_index] = 1.;
    }
  }
  // write data to file
  VTKWriter<GridLayerType> vtkwriter(grid_layer, VTK::nonconforming);
  for (const auto& name_and_data_pair : name_to_data_map)
    vtkwriter.addCellData(name_and_data_pair.second, name_and_data_pair.first);
  vtkwriter.write(filename, VTK::appendedraw);
} // ... visualize_grid(...)


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SAP_HH
