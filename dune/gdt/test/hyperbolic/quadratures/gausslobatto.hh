// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2018)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_GAUSSLOBATTO_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_GAUSSLOBATTO_HH

#include <numeric>
#include <vector>

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/common/debug.hh>
#include <dune/xt/common/float_cmp.hh>

#include <dune/xt/data/gausslobatto/gausslobatto_data.hh>

namespace Dune {
namespace GDT {


template <class FieldType>
class GaussLobattoQuadrature
{
public:
  static Dune::QuadratureRule<FieldType, 1> get(const size_t requested_order)
  {
    size_t order = requested_order;
    if (requested_order > 197) {
      std::cerr << "Warning: Requested gauss lobatto quadrature with order " << requested_order
                << " is not available, using highest available order " << 197 << "." << std::endl;
      order = 197;
    }
    if (order < 2)
      order = 2;
    const size_t num_quad_points = (requested_order + 4) / 2;
    std::vector<std::vector<double>> data_vector;
    switch (num_quad_points) {
      case 2:
        data_vector = XT::Data::internal::GaussLobattoData<2>::get();
        break;
      case 3:
        data_vector = XT::Data::internal::GaussLobattoData<3>::get();
        break;
      case 4:
        data_vector = XT::Data::internal::GaussLobattoData<4>::get();
        break;
      case 5:
        data_vector = XT::Data::internal::GaussLobattoData<5>::get();
        break;
      case 6:
        data_vector = XT::Data::internal::GaussLobattoData<6>::get();
        break;
      case 7:
        data_vector = XT::Data::internal::GaussLobattoData<7>::get();
        break;
      case 8:
        data_vector = XT::Data::internal::GaussLobattoData<8>::get();
        break;
      case 9:
        data_vector = XT::Data::internal::GaussLobattoData<9>::get();
        break;
      case 10:
        data_vector = XT::Data::internal::GaussLobattoData<10>::get();
        break;
      case 11:
        data_vector = XT::Data::internal::GaussLobattoData<11>::get();
        break;
      case 12:
        data_vector = XT::Data::internal::GaussLobattoData<12>::get();
        break;
      case 13:
        data_vector = XT::Data::internal::GaussLobattoData<13>::get();
        break;
      case 14:
        data_vector = XT::Data::internal::GaussLobattoData<14>::get();
        break;
      case 15:
        data_vector = XT::Data::internal::GaussLobattoData<15>::get();
        break;
      case 16:
        data_vector = XT::Data::internal::GaussLobattoData<16>::get();
        break;
      case 17:
        data_vector = XT::Data::internal::GaussLobattoData<17>::get();
        break;
      case 18:
        data_vector = XT::Data::internal::GaussLobattoData<18>::get();
        break;
      case 19:
        data_vector = XT::Data::internal::GaussLobattoData<19>::get();
        break;
      case 20:
        data_vector = XT::Data::internal::GaussLobattoData<20>::get();
        break;
      case 21:
        data_vector = XT::Data::internal::GaussLobattoData<21>::get();
        break;
      case 22:
        data_vector = XT::Data::internal::GaussLobattoData<22>::get();
        break;
      case 23:
        data_vector = XT::Data::internal::GaussLobattoData<23>::get();
        break;
      case 24:
        data_vector = XT::Data::internal::GaussLobattoData<24>::get();
        break;
      case 25:
        data_vector = XT::Data::internal::GaussLobattoData<25>::get();
        break;
      case 26:
        data_vector = XT::Data::internal::GaussLobattoData<26>::get();
        break;
      case 27:
        data_vector = XT::Data::internal::GaussLobattoData<27>::get();
        break;
      case 28:
        data_vector = XT::Data::internal::GaussLobattoData<28>::get();
        break;
      case 29:
        data_vector = XT::Data::internal::GaussLobattoData<29>::get();
        break;
      case 30:
        data_vector = XT::Data::internal::GaussLobattoData<30>::get();
        break;
      case 31:
        data_vector = XT::Data::internal::GaussLobattoData<31>::get();
        break;
      case 32:
        data_vector = XT::Data::internal::GaussLobattoData<32>::get();
        break;
      case 33:
        data_vector = XT::Data::internal::GaussLobattoData<33>::get();
        break;
      case 34:
        data_vector = XT::Data::internal::GaussLobattoData<34>::get();
        break;
      case 35:
        data_vector = XT::Data::internal::GaussLobattoData<35>::get();
        break;
      case 36:
        data_vector = XT::Data::internal::GaussLobattoData<36>::get();
        break;
      case 37:
        data_vector = XT::Data::internal::GaussLobattoData<37>::get();
        break;
      case 38:
        data_vector = XT::Data::internal::GaussLobattoData<38>::get();
        break;
      case 39:
        data_vector = XT::Data::internal::GaussLobattoData<39>::get();
        break;
      case 40:
        data_vector = XT::Data::internal::GaussLobattoData<40>::get();
        break;
      case 41:
        data_vector = XT::Data::internal::GaussLobattoData<41>::get();
        break;
      case 42:
        data_vector = XT::Data::internal::GaussLobattoData<42>::get();
        break;
      case 43:
        data_vector = XT::Data::internal::GaussLobattoData<43>::get();
        break;
      case 44:
        data_vector = XT::Data::internal::GaussLobattoData<44>::get();
        break;
      case 45:
        data_vector = XT::Data::internal::GaussLobattoData<45>::get();
        break;
      case 46:
        data_vector = XT::Data::internal::GaussLobattoData<46>::get();
        break;
      case 47:
        data_vector = XT::Data::internal::GaussLobattoData<47>::get();
        break;
      case 48:
        data_vector = XT::Data::internal::GaussLobattoData<48>::get();
        break;
      case 49:
        data_vector = XT::Data::internal::GaussLobattoData<49>::get();
        break;
      case 50:
        data_vector = XT::Data::internal::GaussLobattoData<50>::get();
        break;
      case 51:
        data_vector = XT::Data::internal::GaussLobattoData<51>::get();
        break;
      case 52:
        data_vector = XT::Data::internal::GaussLobattoData<52>::get();
        break;
      case 53:
        data_vector = XT::Data::internal::GaussLobattoData<53>::get();
        break;
      case 54:
        data_vector = XT::Data::internal::GaussLobattoData<54>::get();
        break;
      case 55:
        data_vector = XT::Data::internal::GaussLobattoData<55>::get();
        break;
      case 56:
        data_vector = XT::Data::internal::GaussLobattoData<56>::get();
        break;
      case 57:
        data_vector = XT::Data::internal::GaussLobattoData<57>::get();
        break;
      case 58:
        data_vector = XT::Data::internal::GaussLobattoData<58>::get();
        break;
      case 59:
        data_vector = XT::Data::internal::GaussLobattoData<59>::get();
        break;
      case 60:
        data_vector = XT::Data::internal::GaussLobattoData<60>::get();
        break;
      case 61:
        data_vector = XT::Data::internal::GaussLobattoData<61>::get();
        break;
      case 62:
        data_vector = XT::Data::internal::GaussLobattoData<62>::get();
        break;
      case 63:
        data_vector = XT::Data::internal::GaussLobattoData<63>::get();
        break;
      case 64:
        data_vector = XT::Data::internal::GaussLobattoData<64>::get();
        break;
      case 65:
        data_vector = XT::Data::internal::GaussLobattoData<65>::get();
        break;
      case 66:
        data_vector = XT::Data::internal::GaussLobattoData<66>::get();
        break;
      case 67:
        data_vector = XT::Data::internal::GaussLobattoData<67>::get();
        break;
      case 68:
        data_vector = XT::Data::internal::GaussLobattoData<68>::get();
        break;
      case 69:
        data_vector = XT::Data::internal::GaussLobattoData<69>::get();
        break;
      case 70:
        data_vector = XT::Data::internal::GaussLobattoData<70>::get();
        break;
      case 71:
        data_vector = XT::Data::internal::GaussLobattoData<71>::get();
        break;
      case 72:
        data_vector = XT::Data::internal::GaussLobattoData<72>::get();
        break;
      case 73:
        data_vector = XT::Data::internal::GaussLobattoData<73>::get();
        break;
      case 74:
        data_vector = XT::Data::internal::GaussLobattoData<74>::get();
        break;
      case 75:
        data_vector = XT::Data::internal::GaussLobattoData<75>::get();
        break;
      case 76:
        data_vector = XT::Data::internal::GaussLobattoData<76>::get();
        break;
      case 77:
        data_vector = XT::Data::internal::GaussLobattoData<77>::get();
        break;
      case 78:
        data_vector = XT::Data::internal::GaussLobattoData<78>::get();
        break;
      case 79:
        data_vector = XT::Data::internal::GaussLobattoData<79>::get();
        break;
      case 80:
        data_vector = XT::Data::internal::GaussLobattoData<80>::get();
        break;
      case 81:
        data_vector = XT::Data::internal::GaussLobattoData<81>::get();
        break;
      case 82:
        data_vector = XT::Data::internal::GaussLobattoData<82>::get();
        break;
      case 83:
        data_vector = XT::Data::internal::GaussLobattoData<83>::get();
        break;
      case 84:
        data_vector = XT::Data::internal::GaussLobattoData<84>::get();
        break;
      case 85:
        data_vector = XT::Data::internal::GaussLobattoData<85>::get();
        break;
      case 86:
        data_vector = XT::Data::internal::GaussLobattoData<86>::get();
        break;
      case 87:
        data_vector = XT::Data::internal::GaussLobattoData<87>::get();
        break;
      case 88:
        data_vector = XT::Data::internal::GaussLobattoData<88>::get();
        break;
      case 89:
        data_vector = XT::Data::internal::GaussLobattoData<89>::get();
        break;
      case 90:
        data_vector = XT::Data::internal::GaussLobattoData<90>::get();
        break;
      case 91:
        data_vector = XT::Data::internal::GaussLobattoData<91>::get();
        break;
      case 92:
        data_vector = XT::Data::internal::GaussLobattoData<92>::get();
        break;
      case 93:
        data_vector = XT::Data::internal::GaussLobattoData<93>::get();
        break;
      case 94:
        data_vector = XT::Data::internal::GaussLobattoData<94>::get();
        break;
      case 95:
        data_vector = XT::Data::internal::GaussLobattoData<95>::get();
        break;
      case 96:
        data_vector = XT::Data::internal::GaussLobattoData<96>::get();
        break;
      case 97:
        data_vector = XT::Data::internal::GaussLobattoData<97>::get();
        break;
      case 98:
        data_vector = XT::Data::internal::GaussLobattoData<98>::get();
        break;
      case 99:
        data_vector = XT::Data::internal::GaussLobattoData<99>::get();
        break;
      case 100:
        data_vector = XT::Data::internal::GaussLobattoData<100>::get();
        break;
      default:
        DUNE_THROW(NotImplemented, "Requested order is not available!");
    }
    DXT_ASSERT(data_vector.size() == num_quad_points);
    Dune::QuadratureRule<FieldType, 1> ret;
    // stored values are for [-1, 1], so convert for [0, 1]
    for (const auto& pair : data_vector)
      ret.emplace_back((pair[0] + 1) / 2, pair[1] / 2);
#ifndef NDEBUG
    if (requested_order <= 31) {
      auto reference_quadrature =
          Dune::QuadratureRules<FieldType, 1>::rule(Dune::GeometryType(Dune::GeometryType::BasicType::simplex, 1),
                                                    static_cast<int>(requested_order),
                                                    Dune::QuadratureType::GaussLobatto);
      std::sort(reference_quadrature.begin(),
                reference_quadrature.end(),
                [](const Dune::QuadraturePoint<FieldType, 1>& a, const Dune::QuadraturePoint<FieldType, 1>& b) {
                  return a.position()[0] < b.position()[0];
                });
      for (size_t ii = 0; ii < num_quad_points; ++ii) {
        DXT_ASSERT(XT::Common::FloatCmp::eq(ret[ii].position(), reference_quadrature[ii].position()));
        DXT_ASSERT(XT::Common::FloatCmp::eq(ret[ii].weight(), reference_quadrature[ii].weight()));
      }
    }
    // check sanity of quadrature
    const FieldType summed_weights = std::accumulate(
        ret.begin(), ret.end(), 0., [](const FieldType& sum, const Dune::QuadraturePoint<FieldType, 1> quad_point) {
          return sum + quad_point.weight();
        });
    DXT_ASSERT(XT::Common::FloatCmp::eq(summed_weights, 1.));
#endif
    return ret;
  }
}; // class GaussLobattoQuadrature


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_GAUSSLOBATTO_HH
