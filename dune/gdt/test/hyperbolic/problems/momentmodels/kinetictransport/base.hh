// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_KINETICTRANSPORTEQUATION_BASE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_KINETICTRANSPORTEQUATION_BASE_HH

#include <dune/grid/common/partitionset.hh>

#include <dune/gdt/test/hyperbolic/problems/momentmodels/basisfunctions/base.hh>

#include "../kineticequation.hh"

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


template <class BasisfunctionImp, class GridLayerImp, class U_>
class KineticTransportEquation
  : public KineticEquationImplementationInterface<BasisfunctionImp, GridLayerImp, U_>
  , public XT::Common::ParametricInterface
{
  using ThisType = KineticTransportEquation<BasisfunctionImp, GridLayerImp, U_>;
  using BaseType = KineticEquationImplementationInterface<BasisfunctionImp, GridLayerImp, U_>;

public:
  using BaseType::dimDomain;
  using BaseType::dimFlux;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  using typename BaseType::BasisfunctionType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::GridLayerType;
  using typename BaseType::IntersectionType;
  using typename BaseType::MatrixType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::StateType;

  using typename BaseType::ActualBoundaryValueType;
  using typename BaseType::ActualDirichletBoundaryValueType;
  using typename BaseType::ActualFluxType;
  using typename BaseType::ActualInitialValueType;
  using typename BaseType::ActualRhsType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::DirichletBoundaryValueType;
  using typename BaseType::FluxType;
  using typename BaseType::InitialValueType;
  using typename BaseType::RhsAffineFunctionType;
  using typename BaseType::RhsType;
  using CheckerboardType =
      XT::Functions::CheckerboardFunction<typename U_::EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1>;
  using ConstantType =
      XT::Functions::ConstantFunction<typename U_::EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1>;

  using BaseType::default_boundary_cfg;
  using BaseType::default_grid_cfg;
  KineticTransportEquation(const BasisfunctionType& basis_functions,
                           const GridLayerType& grid_layer,
                           DynamicVector<size_t> num_segments = {1, 1, 1},
                           const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                           const XT::Common::Configuration& boundary_cfg = default_boundary_cfg(),
                           const RangeFieldType psi_vac = 5e-9,
                           const XT::Common::ParameterType& parameter_type = XT::Common::ParameterType())
    : BaseType(basis_functions, grid_layer)
    , num_segments_(num_segments)
    , grid_cfg_(grid_cfg)
    , boundary_cfg_(boundary_cfg)
    , psi_vac_(psi_vac)
    , parameter_type_(parameter_type)
  {
    // by default, initially num_segments has size 3. If dimDomain < 3 we thus need to resize
    num_segments_.resize(dimDomain);
    if (parameter_type_.empty())
      parameter_type_ = XT::Common::ParameterType({std::make_pair("sigma_a", get_num_regions(num_segments_)),
                                                   std::make_pair("sigma_s", get_num_regions(num_segments_)),
                                                   std::make_pair("Q", get_num_regions(num_segments_)),
                                                   std::make_pair("CFL", 1),
                                                   std::make_pair("t_end", 1)});
  }

  virtual ~KineticTransportEquation() override {}

  virtual bool is_parametric() const override
  {
    return true;
  }

  virtual const XT::Common::ParameterType& parameter_type() const override
  {
    return parameter_type_;
  }

  virtual XT::Common::Parameter parameters() const = 0;

  template <class VectorType>
  void solve(const MatrixType& mat,
             VectorType& x,
             const VectorType& rhs,
             const BasisfunctionsInterface<DomainFieldType,
                                           BasisfunctionType::dimDomain,
                                           RangeFieldType,
                                           dimRange,
                                           dimRangeCols,
                                           dimFlux>&) const
  {
    // copy to our FieldMatrix as the DenseMatrix in dune-common has a bug in its solve method (fixed in 2.6)
    auto field_mat = std::make_unique<XT::Common::FieldMatrix<RangeFieldType, dimRange, dimRange>>(mat);
    field_mat->solve(x, rhs);
  }

  template <class VectorType, size_t refinements>
  void solve(const MatrixType& mat,
             VectorType& x,
             const VectorType& rhs,
             const PartialMomentBasis<DomainFieldType, 3, RangeFieldType, refinements, dimRangeCols, 3, 1>&) const
  {
    const size_t num_blocks = dimRange / 4;
    const size_t block_size = 4;
    XT::Common::FieldMatrix<RangeFieldType, block_size, block_size> local_mat;
    XT::Common::FieldVector<RangeFieldType, block_size> local_x, local_rhs;
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      const size_t offset = jj * block_size;
      // copy to local matrix and vector
      for (size_t rr = 0; rr < block_size; ++rr) {
        local_rhs[rr] = rhs[offset + rr];
        for (size_t cc = 0; cc < block_size; ++cc)
          local_mat[rr][cc] = mat[offset + rr][offset + cc];
      } // rr
      local_mat.solve(local_x, local_rhs);
      for (size_t rr = 0; rr < block_size; ++rr)
        x[offset + rr] = local_x[rr];
    } // jj
  }

  using XT::Common::ParametricInterface::parse_parameter;

  // flux matrix A = B M^{-1} with B_{ij} = <v h_i h_j>
  virtual FluxType* create_flux() const override
  {
    // calculate B row-wise by solving M^{T} A^T = B^T column-wise
    auto A = basis_functions_.mass_matrix_with_v();
    const auto M_T = basis_functions_.mass_matrix(); // mass matrix is symmetric
    // solve
    DynamicVector<RangeFieldType> tmp_row(M_T.N(), 0.);
    for (size_t dd = 0; dd < dimDomain; ++dd) {
      for (size_t ii = 0; ii < M_T.N(); ++ii) {
        solve(M_T, tmp_row, A[dd][ii], basis_functions_);
        A[dd][ii] = tmp_row;
      }
    }
    return new ActualFluxType(A, typename ActualFluxType::RangeType(0));
  }

  // RHS is (sigma_s/vol*G - sigma_t * I)u + Q<b>,
  // where sigma_t = sigma_s + sigma_a, G = <b><b>^T M^{-1} = <b>*c^T and
  // vol = <1> is the volume of the integration domain.
  virtual RhsType* create_rhs() const override
  {
    const auto param = parse_parameter(parameters());
    const size_t num_regions = get_num_regions(num_segments_);
    const auto sigma_a = param.get("sigma_a");
    const auto sigma_s = param.get("sigma_s");
    const auto Q = param.get("Q");
    assert(sigma_a.size() == sigma_s.size() && sigma_a.size() == Q.size() && sigma_a.size() == num_regions);
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    auto sigma_t = sigma_a;
    for (size_t ii = 0; ii < num_regions; ++ii)
      sigma_t[ii] += sigma_s[ii];
    const RangeType basis_integrated = basis_functions_.integrated(true);
    // calculate c = M^{-T} <b>
    const auto M_T = basis_functions_.mass_matrix(true); // mass matrix is symmetric
    RangeType c(0.);
    solve(M_T, c, basis_integrated, basis_functions_);
    MatrixType I(dimRange, dimRange, 0.);
    for (size_t rr = 0; rr < dimRange; ++rr)
      I[rr][rr] = 1;
    MatrixType G(dimRange, dimRange, 0.);
    for (size_t rr = 0; rr < dimRange; ++rr)
      for (size_t cc = 0; cc < dimRange; ++cc)
        G[rr][cc] = basis_integrated[rr] * c[cc];
    const auto vol = basis_functions_.unit_ball_volume(true);
    std::vector<RhsAffineFunctionType> affine_functions;
    for (size_t ii = 0; ii < num_regions; ++ii) {
      MatrixType G_scaled = G;
      G_scaled *= sigma_s[ii] / vol;
      MatrixType I_scaled = I;
      I_scaled *= sigma_t[ii];
      MatrixType A = G_scaled;
      A -= I_scaled;
      RangeType b = basis_integrated;
      b *= Q[ii];
      affine_functions.emplace_back(A, b, true, "rhs");
    } // ii
    return new ActualRhsType(lower_left, upper_right, num_segments_, affine_functions);
  } // ... create_rhs(...)

  virtual CheckerboardType get_sigma_a() const
  {
    const auto param = parse_parameter(parameters());
    const auto sigma_a = param.get("sigma_a");
    const auto num_regions = get_num_regions(num_segments_);
    std::vector<ConstantType> sigma_a_funcs;
    for (size_t ii = 0; ii < num_regions; ++ii)
      sigma_a_funcs.emplace_back(sigma_a[ii]);
    assert(sigma_a.size() == num_regions);
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    return CheckerboardType(lower_left, upper_right, num_segments_, sigma_a_funcs);
  } // ... sigma_a(...)

  virtual CheckerboardType get_sigma_s() const
  {
    const auto param = parse_parameter(parameters());
    const auto sigma_s = param.get("sigma_s");
    const auto num_regions = get_num_regions(num_segments_);
    std::vector<ConstantType> sigma_s_funcs;
    for (size_t ii = 0; ii < num_regions; ++ii)
      sigma_s_funcs.emplace_back(sigma_s[ii]);
    assert(sigma_s.size() == num_regions);
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    return CheckerboardType(lower_left, upper_right, num_segments_, sigma_s_funcs);
  } // ... sigma_s(...)

  virtual CheckerboardType get_Q() const
  {
    const auto param = parse_parameter(parameters());
    const auto Q = param.get("Q");
    const auto num_regions = get_num_regions(num_segments_);
    std::vector<ConstantType> Q_funcs;
    for (size_t ii = 0; ii < num_regions; ++ii)
      Q_funcs.emplace_back(Q[ii]);
    assert(Q.size() == num_regions);
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    return CheckerboardType(lower_left, upper_right, num_segments_, Q_funcs);
  } // ... Q(...)

  // Initial value of the kinetic equation is a constant vacuum concentration psi_vac.
  // Thus, the initial value of the n-th moment is basis_integrated * psi_vac.
  virtual InitialValueType* create_initial_values() const override
  {
    const DomainType lower_left = XT::Common::from_string<DomainType>(grid_cfg_["lower_left"]);
    const DomainType upper_right = XT::Common::from_string<DomainType>(grid_cfg_["upper_right"]);
    RangeType value = basis_functions_.integrated() * psi_vac_;
    std::vector<typename ActualInitialValueType::LocalizableFunctionType> initial_vals;
    const size_t num_regions = get_num_regions(num_segments_);
    for (size_t ii = 0; ii < num_regions; ++ii)
      initial_vals.emplace_back([=](const DomainType&, const XT::Common::Parameter&) { return value; }, 0);
    return new ActualInitialValueType(lower_left, upper_right, num_segments_, initial_vals, "initial_values");
  } // ... create_initial_values()

  // Use a constant vacuum concentration basis_integrated * psi_vac as boundary value
  virtual BoundaryValueType* create_boundary_values() const override
  {
    RangeType value = basis_functions_.integrated() * psi_vac_;
    return new ActualBoundaryValueType(XT::Grid::make_alldirichlet_boundaryinfo<IntersectionType>(),
                                       std::make_unique<ActualDirichletBoundaryValueType>(
                                           [=](const DomainType&, const XT::Common::Parameter&) { return value; }, 0));
  } // ... create_boundary_values()

  virtual RangeFieldType CFL() const override
  {
    return parameters().get("CFL")[0];
  }

  virtual RangeFieldType t_end() const override
  {
    return parameters().get("t_end")[0];
  }

  virtual XT::Common::Configuration grid_config() const override
  {
    return grid_cfg_;
  }

  virtual XT::Common::Configuration boundary_config() const override
  {
    return boundary_cfg_;
  }

  static std::string static_id()
  {
    return "kinetictransportequation";
  }

  virtual const RangeFieldType psi_vac() const
  {
    return psi_vac_;
  }

protected:
  static size_t get_num_regions(const DynamicVector<size_t>& num_segments)
  {
    return std::accumulate(num_segments.begin(), num_segments.end(), size_t(1), [](auto a, auto b) { return a * b; });
  }

  using BaseType::basis_functions_;
  DynamicVector<size_t> num_segments_;
  const XT::Common::Configuration grid_cfg_;
  const XT::Common::Configuration boundary_cfg_;
  const RangeFieldType psi_vac_;
  XT::Common::ParameterType parameter_type_;
}; // class KineticTransportEquation<...>


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_KINETICTRANSPORTEQUATION_BASE_HH
