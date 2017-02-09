// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_POINTSOURCE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_POINTSOURCE_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/xt/common/string.hh>

#include <dune/xt/functions/global.hh>

#include "planesource.hh"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Origin.h>

struct CGALWrapper
{
  // A vertex type with indices.
  template <class Refs, class Traits>
  struct My_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, typename Traits::Point_3>
  {
    typedef typename CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, typename Traits::Point_3> BaseType;
    typedef typename Traits::Point_3 Point_3;

    template <class... Args>
    My_vertex(Args&&... args)
      : BaseType(std::forward<Args>(args)...)
    {
    }

    size_t index;
  };

  // A face type with indices.
  template <class Refs>
  struct My_face : public CGAL::HalfedgeDS_face_base<Refs>
  {
    size_t index;
  };

  // An items type using the vertex and face type with indices.
  struct My_items : public CGAL::Polyhedron_items_3
  {
    template <class Refs, class Traits>
    struct Vertex_wrapper
    {
      typedef typename Traits::Point_3 Point;
      typedef My_vertex<Refs, Traits> Vertex;
    };
    template <class Refs, class Traits>
    struct Face_wrapper
    {
      typedef typename Traits::Plane_3 Plane;
      typedef My_face<Refs> Face;
    };
  };

  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Polyhedron_3<K, My_items> Polyhedron_3;
  typedef typename Polyhedron_3::Vertex_const_handle VertexHandleType;
  typedef typename Polyhedron_3::Facet_const_handle FacetHandleType;
  // define point creator
  typedef K::Point_3 Point_3;
  typedef K::Vector_3 Vector_3;

  static Polyhedron_3 create_spherical_triangulation(std::vector<Point_3> points, const size_t num_refinements = 0)

  {
    // define polyhedron to hold convex hull
    Polyhedron_3 poly;
    // generate initial octaeder
    CGAL::convex_hull_3(points.begin(), points.end(), poly);
    std::cout << "The convex hull contains " << poly.size_of_vertices() << " vertices and " << poly.size_of_facets()
              << " faces" << std::endl;
    for (size_t ii = 0; ii < num_refinements; ++ii) {
      points.clear();
      std::for_each(
          poly.points_begin(), poly.points_end(), [&points](const Point_3& point) { points.push_back(point); });
      // add edge_centers and face centers
      const auto face_it_end = poly.facets_end();
      for (auto face_it = poly.facets_begin(); face_it != face_it_end; ++face_it) {
        assert(face_it->is_triangle());
        // a circulator does not have a past-the-end concept, but starts over at the beginning
        auto halfedge_circ = face_it->facet_begin();
        auto halfedge_circ_begin = halfedge_circ;
        // insert face center
        const Point_3 sum3 = halfedge_circ->vertex()->point()
                             + Vector_3(CGAL::ORIGIN, halfedge_circ->next()->vertex()->point())
                             + Vector_3(CGAL::ORIGIN, halfedge_circ->next()->next()->vertex()->point());
        Dune::FieldVector<double, 3> center({sum3.x() / 3., sum3.y() / 3., sum3.z() / 3.});
        center /= center.two_norm(); // projection onto sphere
        //      points.push_back(Point_3(center[0], center[1], center[2]));
        // insert edge center (inserts each edge center twice, but this shouldn't be a problem due to the convex hull
        // call
        do {
          const Point_3& point = halfedge_circ->vertex()->point();
          const Point_3& next_point = halfedge_circ->next()->vertex()->point();
          const Point_3 sum2 = point + Vector_3(CGAL::ORIGIN, next_point);
          center = {sum2.x() / 2., sum2.y() / 2., sum2.z() / 2.};
          center /= center.two_norm();
          points.push_back(Point_3(center[0], center[1], center[2]));
        } while (++halfedge_circ != halfedge_circ_begin);
      } // iterate over faces
      poly.clear();
      CGAL::convex_hull_3(points.begin(), points.end(), poly);
      std::cout << "The convex hull contains " << poly.size_of_vertices() << " vertices and " << poly.size_of_facets()
                << " faces" << std::endl;
    } // refinement loop

    // add indices to vertices and facets
    size_t index = 0;
    const auto vertices_it_end = poly.vertices_end();
    for (auto vertices_it = poly.vertices_begin(); vertices_it != vertices_it_end; ++vertices_it, ++index)
      vertices_it->index = index;
    index = 0;
    const auto facets_it_end = poly.facets_end();
    for (auto facets_it = poly.facets_begin(); facets_it != facets_it_end; ++facets_it, ++index)
      facets_it->index = index;
    return poly;
  } // create_spherical_triangulation(...)

  static Polyhedron_3 create_octaeder_spherical_triangulation(const size_t num_refinements = 0)
  {
    // vertices of the octaeder
    std::vector<Point_3> points{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {-1, 0, 0}, {0, -1, 0}, {0, 0, -1}};
    return create_spherical_triangulation(points, num_refinements);
  }
}; // struct CGALWrapper


namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {

inline std::string trim(const std::string& s)
{
  auto wsfront = std::find_if_not(s.begin(), s.end(), [](int c) { return std::isspace(c); });
  auto wsback = std::find_if_not(s.rbegin(), s.rend(), [](int c) { return std::isspace(c); }).base();
  return (wsback <= wsfront ? std::string() : std::string(wsfront, wsback));
}

static const std::vector<int> allowed_degrees = {6,    14,   26,   38,   50,   74,   86,   110,  146,  170,  194,
                                                 230,  266,  302,  350,  434,  590,  770,  974,  1202, 1454, 1730,
                                                 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810};
static const std::vector<int> allowed_orders = {3,  5,  7,  9,  11, 13, 15, 17, 19, 21, 23,  25,  27,  29,  31,  35,
                                                41, 47, 53, 59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131};

Dune::QuadratureRule<double, 2> get_lebedev_quadrature(size_t requested_order)
{
  size_t index = -1;
  for (int ii = 0; ii < allowed_orders.size(); ++ii) {
    if (allowed_orders[ii] >= requested_order) {
      index = ii;
      break;
    }
  }
  int order = allowed_orders[index];
  //  int degree = allowed_degrees[index];
  char orderstring[4];
  sprintf(orderstring, "%03d", order);
  std::string filename =
      std::string("/home/tobias/Software/dune-gdt-super-2.5/dune-gdt/dune/gdt/LebedevTables/lebedev_") + orderstring
      + ".txt";

  Dune::QuadratureRule<double, 2> quadrature_rule;
  std::string current_line;
  std::ifstream quadrature_file(filename);
  while (getline(quadrature_file, current_line)) {
    current_line = trim(current_line);
    auto quadrature_values =
        XT::Common::tokenize(current_line, " ", boost::algorithm::token_compress_mode_type::token_compress_on);
    std::cout << XT::Common::to_string(quadrature_values) << std::endl;
    const double phi = XT::Common::from_string<double>(quadrature_values[0]) / 360. * 2 * M_PI;
    const double theta = XT::Common::from_string<double>(quadrature_values[1]) / 360. * 2 * M_PI;
    const double mu = std::cos(theta);
    const double weight = 4 * M_PI * XT::Common::from_string<double>(quadrature_values[2]);
    Dune::QuadraturePoint<double, 2> quadrature_point(FieldVector<double, 2>{mu, phi}, weight);
    quadrature_rule.push_back(quadrature_point);
  }
  return quadrature_rule;
}

template <class DomainType, class VertexVectorType>
bool calculate_barycentric_coordinates(const DomainType& v,
                                       const VertexVectorType& vertices,
                                       Dune::FieldVector<double, 3>& ret)
{
  Dune::FieldMatrix<double, 3, 3> gradients(0);
  for (size_t ii = 0; ii < 3; ++ii) {
    const auto& point = vertices[ii]->point();
    // copy points to gradients
    gradients[ii] = Dune::FieldVector<double, 3>({point.x(), point.y(), point.z()});
    const auto scalar_prod = v * gradients[ii];
    // if v is not on the same octant of the sphere as the vertices, return false
    // assumes the triangulation is fine enough that vertices[ii]*vertices[jj] >= 0 for all triangles
    if (XT::Common::FloatCmp::lt(scalar_prod, 0.))
      return false;
    auto v_scaled = v;
    v_scaled *= scalar_prod;
    gradients[ii] -= v_scaled;
    // scale with factor
    gradients[ii] *= std::acos(scalar_prod) / std::sqrt(1. - std::pow(scalar_prod, 2));
  } // ii
  // calculate barycentric coordinates for 0 w.r.t the points g_i
  const auto& g0 = gradients[0];
  const auto& g1 = gradients[1];
  const auto& g2 = gradients[2];
  auto g0_minus_g2 = g0;
  auto g1_minus_g2 = g1;
  g0_minus_g2 -= g2;
  g1_minus_g2 -= g2;
  Dune::FieldMatrix<double, 2, 2> A;
  Dune::FieldVector<double, 2> solution;
  Dune::FieldVector<double, 2> rhs;
  // (ii, jj) = (0, 1), (0, 2), (1, 2)
  for (size_t ii = 0; ii < 2; ++ii) {
    for (size_t jj = ii + 1; jj < 3; ++jj) {
      A[0][0] = g0_minus_g2[ii];
      A[1][0] = g0_minus_g2[jj];
      A[0][1] = g1_minus_g2[ii];
      A[1][1] = g1_minus_g2[jj];
      double det = A.determinant();
      if (XT::Common::FloatCmp::eq(det, 0.))
        break;
      rhs[0] = -g2[ii];
      rhs[1] = -g2[jj];
      A.solve(solution, rhs);
      if (XT::Common::FloatCmp::lt(solution[0], 0.) || XT::Common::FloatCmp::lt(solution[1], 0.))
        return false;
      ret[0] = solution[0];
      ret[1] = solution[1];
      ret[2] = 1. - ret[0] - ret[1];
      if (XT::Common::FloatCmp::lt(ret[2], 0.))
        return false;
      return true;
    }
  }
  return false;
}

Dune::QuadratureRule<double, 3>
get_barycentre_rule(const Dune::XT::Common::FieldVector<Dune::XT::Common::FieldVector<double, 3>, 3>& vertices)
{
  double s = 1. / 3.;
  double t = 1. / 3.;
  const auto ff = vertices[0] + (vertices[1] - vertices[0]) * t + (vertices[2] - vertices[0]) * s;
  const auto norm_ff = ff.two_norm();
  const auto bb = ff / norm_ff;
  const auto partial_s_gg =
      (vertices[2] - vertices[0]) / norm_ff - ff * (ff * (vertices[2] - vertices[0]) / std::pow(norm_ff, 3));
  const auto partial_t_gg =
      (vertices[1] - vertices[0]) / norm_ff - ff * (ff * (vertices[1] - vertices[0]) / std::pow(norm_ff, 3));
  const auto weight = Dune::PDELab::crossproduct(partial_s_gg, partial_t_gg).two_norm() / 2.;
  Dune::QuadratureRule<double, 3> ret;
  ret.push_back(Dune::QuadraturePoint<double, 3>(bb, weight));
  return ret;
}

template <class RangeFieldImp = double>
class SphericalTriangle
{
  typedef SphericalTriangle<RangeFieldImp> ThisType;

public:
  typedef typename XT::Common::FieldVector<RangeFieldImp, 3> VertexType;
  typedef typename Dune::FieldVector<VertexType, 3> VertexVectorType;
  typedef typename Dune::FieldVector<ThisType, 4> SubtrianglesVectorType;
  typedef typename Dune::QuadraturePoint<RangeFieldImp, 3> QuadraturePointType;
  typedef typename Dune::QuadratureRule<RangeFieldImp, 3> QuadratureRuleType;

  SphericalTriangle() = default;

  SphericalTriangle(const VertexVectorType& vertices)
    : vertices_(vertices)
  {
  }

  SphericalTriangle(const VertexType& vertex_1, const VertexType& vertex_2, const VertexType& vertex_3)
    : vertices_{vertex_1, vertex_2, vertex_3}
  {
  }

  SphericalTriangle(const VertexType& vertex_1,
                    const VertexType& vertex_2,
                    const VertexType& vertex_3,
                    const std::shared_ptr<QuadratureRuleType>& barycentre_rule)
    : vertices_{vertex_1, vertex_2, vertex_3}
    , barycentre_rule_(barycentre_rule)
  {
  }

  const SubtrianglesVectorType& get_subtriangles() const
  {
    initialize_subtriangles();
    return *subtriangles_;
  }

  const QuadraturePointType& get_quadrature_point() const
  {
    calculate_barycentre_rule();
    return (*barycentre_rule_)[0];
  }

private:
  void initialize_subtriangles() const
  {
    if (!subtriangles_) {
      VertexVectorType midpoints;
      for (size_t ii = 0; ii < 3; ++ii) {
        midpoints[ii] = vertices_[(1 + ii) % 3] + vertices_[(2 + ii) % 3];
        midpoints[ii] /= midpoints[ii].two_norm();
      }
      subtriangles_ = std::make_shared<SubtrianglesVectorType>();
      auto& subtriangles = *subtriangles_;
      subtriangles[0] = ThisType(vertices_[0], midpoints[1], midpoints[2]);
      subtriangles[1] = ThisType(vertices_[1], midpoints[2], midpoints[0]);
      subtriangles[2] = ThisType(vertices_[2], midpoints[0], midpoints[1]);
      subtriangles[3] = ThisType(midpoints[0], midpoints[1], midpoints[2], barycentre_rule_);
    }
  } // initialize_subtriangles()

  void calculate_barycentre_rule() const
  {
    if (!barycentre_rule_) {
      const auto ff = (vertices_[0] + vertices_[1] + vertices_[2]) / 3.;
      const auto norm_ff = ff.two_norm();
      const auto bb = ff / norm_ff;
      const auto norm_ff_3 = std::pow(norm_ff, 3);
      const auto vertices_1_minus_0 = vertices_[1] - vertices_[0];
      const auto vertices_2_minus_0 = vertices_[2] - vertices_[0];
      const auto partial_s_gg = vertices_2_minus_0 / norm_ff - ff * (ff * vertices_2_minus_0 / norm_ff_3);
      const auto partial_t_gg = vertices_1_minus_0 / norm_ff - ff * (ff * vertices_1_minus_0 / norm_ff_3);
      const auto weight = Dune::PDELab::crossproduct(partial_s_gg, partial_t_gg).two_norm() / 2.;
      barycentre_rule_ = std::make_shared<QuadratureRuleType>();
      barycentre_rule_->emplace_back(bb, weight);
    }
  }

  VertexVectorType vertices_;
  mutable std::shared_ptr<SubtrianglesVectorType> subtriangles_;
  mutable std::shared_ptr<QuadratureRuleType> barycentre_rule_;
};


template <class RangeFieldImp = double>
class SphericalTriangulation : std::vector<SphericalTriangle<RangeFieldImp>>
{
public:
  typedef SphericalTriangle<RangeFieldImp> TriangleType;
  typedef typename TriangleType::VertexVectorType VertexVectorType;
  typedef typename VertexVectorType::value_type VertexType;

  SphericalTriangulation(const typename CGALWrapper::Polyhedron_3& poly)
  {
    const auto facets_it_end = poly.facets_end();
    for (auto facets_it = poly.facets_begin(); facets_it != facets_it_end; ++facets_it) {
      VertexVectorType vertices;
      auto halfedge_circ = facets_it->facet_begin();
      auto halfedge_circ_begin = halfedge_circ;
      size_t ii = 0;
      do {
        const auto& point = halfedge_circ->vertex()->point();
        vertices[ii] = VertexType({point.x(), point.y(), point.z()});
      } while (++ii, ++halfedge_circ != halfedge_circ_begin);
      faces_.emplace_back(vertices);
    } // iterate over faces
  }

  const std::vector<TriangleType>& get_faces() const
  {
    return faces_;
  }

  std::vector<TriangleType>& get_faces()
  {
    return faces_;
  }

private:
  std::vector<TriangleType> faces_;
};

template <class DomainType, class RangeType, class DataType = double, class RangeFieldImp = double>
class AdaptiveQuadrature
{
public:
  typedef SphericalTriangulation<RangeFieldImp> TriangulationType;
  typedef SphericalTriangle<RangeFieldImp> TriangleType;
  typedef typename TriangleType::VertexVectorType VertexVectorType;
  typedef typename VertexVectorType::value_type VertexType;

  AdaptiveQuadrature(const typename CGALWrapper::Polyhedron_3& poly, const double tol = 1e-4, const double gamma = 1.)
    : triangulation_(poly)
    , tol_(tol)
    , gamma_(gamma)
  {
    for (auto& face : triangulation_.get_faces())
      current_triangles_.push_back(&face);
  }

  RangeType calculate_integral(std::function<RangeType(const DomainType&, const DataType&)> psi,
                               std::function<DataType(const DomainType&)> data_function)
  {
    // get variables in thread and resize if necessary
    //    auto& I_0_vector = *I_0_vector_;
    //    auto& I_1_vector = *I_1_vector_;
    //    auto& error_vector = *error_vector_;
    //    auto& norm_of_errors_vector = *norm_of_errors_vector_;
    //    auto& sum_of_norm_of_errors = *sum_of_norm_of_errors_;
    //    auto& error = *error_;
    //    auto& result = *result_;
    //    auto& data = *data_;
    //    auto& data_stored = *data_stored_;
    //    auto& subtriangle_data = *subtriangle_data_;
    //    auto& subtriangle_data_stored = *subtriangle_data_stored_;

    while (true) {
      sum_of_norm_of_errors_ = 0;
      size_t num_triangles = current_triangles_.size();
      std::cout << num_triangles << std::endl;
      if (num_triangles > I_0_vector_.size()) {
        I_0_vector_.resize(num_triangles);
        I_1_vector_.resize(num_triangles);
        error_vector_.resize(num_triangles);
        norm_of_errors_vector_.resize(num_triangles);
        data_.resize(num_triangles);
        data_stored_.resize(num_triangles, false);
        subtriangle_data_.resize(num_triangles);
        subtriangle_data_stored_.resize(num_triangles, false);
      }
      // loop over triangles
      for (size_t ii = 0; ii < num_triangles; ++ii) {
        auto& I_0 = I_0_vector_[ii];
        auto& I_1 = I_1_vector_[ii];
        auto& local_error = error_vector_[ii];
        auto& triangle = *(current_triangles_[ii]);
        // calculate I_0
        const auto& barycentre = triangle.get_quadrature_point();
        if (!data_stored_[ii]) {
          data_[ii] = data_function(barycentre.position());
          data_stored_[ii] = true;
        }
        I_0 = barycentre.weight() * psi(barycentre.position(), data_[ii]);
        // calculate I_1
        I_1 *= 0;
        for (const auto& subtriangle : triangle.get_subtriangles()) {
          const auto& subbarycentre = subtriangle.get_quadrature_point();
          if (!subtriangle_data_stored_[ii]) {
            subtriangle_data_[ii] = data_function(subbarycentre.position());
            subtriangle_data_stored_[ii] = true;
          }
          I_1 += subbarycentre.weight() * psi(subbarycentre.position(), data_[ii]);
        }
        // calculate E(K)
        local_error = I_1;
        local_error -= I_0;
        norm_of_errors_vector_[ii] = local_error.two_norm() * 4. / 3.;
        sum_of_norm_of_errors_ += norm_of_errors_vector_[ii];
      } // loop over triangles;

      error_ *= 0;
      result_ *= 0;
      for (size_t ii = 0; ii < num_triangles; ++ii) {
        error_ += error_vector_[ii];
        result_ += I_1_vector_[ii];
      }
      const auto error_norm = error_.two_norm();
      const auto result_norm = result_.two_norm();
      if (error_norm < tol_ * result_norm)
        break;

      for (size_t ii = 0; ii < num_triangles; ++ii) {
        if (norm_of_errors_vector_[ii] > gamma_ / num_triangles * sum_of_norm_of_errors_) {
          const auto& subtriangles = current_triangles_[ii]->get_subtriangles();
          // no need to transfer data for the midpoint subtriangle, as it is the same as for
          // the father triangle
          current_triangles_[ii] = &subtriangles[3];
          subtriangle_data_stored_[ii] = false;
          for (size_t jj = 0; jj < 3; ++jj) {
            current_triangles_.push_back(&subtriangles[jj]);
            data_.push_back(subtriangle_data_[ii][jj]);
          }
        }
      }
      // if no triangles where added, refine all triangles as error is evenly distributed
      if (current_triangles_.size() == num_triangles) {
        for (size_t ii = 0; ii < num_triangles; ++ii) {
          auto& subtriangles = current_triangles_[ii].get_subtriangles();
          current_triangles_[ii] = &subtriangles[3];
          subtriangle_data_stored_[ii] = false;
          for (size_t jj = 0; jj < 3; ++jj) {
            current_triangles_.push_back(&subtriangles[jj]);
            data_.push_back(subtriangle_data_[ii][jj]);
          }
        }
      }
    } // while(true)
    return result_;
  } // void calculate_integral

  void clear_data()
  {
    data_.clear();
    data_stored_.clear();
    subtriangle_data_.clear();
    subtriangle_data_stored_.clear();
  }

private:
  TriangulationType triangulation_;
  std::vector<RangeType> I_0_vector_, I_1_vector_, error_vector_;
  std::vector<RangeFieldImp> norm_of_errors_vector_;
  RangeType error_, result_;
  RangeFieldImp sum_of_norm_of_errors_;
  double tol_;
  double gamma_;
  std::vector<TriangleType*> current_triangles_;
  std::vector<DataType> data_;
  std::vector<bool> data_stored_;
  std::vector<Dune::FieldVector<DataType, 3>> subtriangle_data_;
  std::vector<bool> subtriangle_data_stored_;
};


template <class DomainType, class DataType, class RangeFieldImp>
class AdaptiveQuadrature<DomainType, double, DataType, RangeFieldImp>
{
public:
  typedef double RangeType;
  typedef SphericalTriangulation<RangeFieldImp> TriangulationType;
  typedef SphericalTriangle<RangeFieldImp> TriangleType;
  typedef typename TriangleType::VertexVectorType VertexVectorType;
  typedef typename VertexVectorType::value_type VertexType;

  AdaptiveQuadrature(const typename CGALWrapper::Polyhedron_3& poly, const double tol = 1e-4, const double gamma = 1.)
    : triangulation_(poly)
    , tol_(tol)
    , gamma_(gamma)
  {
    for (auto& face : triangulation_.get_faces())
      current_triangles_.push_back(&face);
  }

  RangeType calculate_integral(std::function<RangeType(const DomainType&, const DataType&)> psi,
                               std::function<DataType(const DomainType&)> data_function)
  {
    // get variables in thread and resize if necessary
    //    auto& I_0_vector = *I_0_vector_;
    //    auto& I_1_vector = *I_1_vector_;
    //    auto& error_vector = *error_vector_;
    //    auto& norm_of_errors_vector = *norm_of_errors_vector_;
    //    auto& sum_of_norm_of_errors = *sum_of_norm_of_errors_;
    //    auto& error = *error_;
    //    auto& result = *result_;
    //    auto& data = *data_;
    //    auto& data_stored = *data_stored_;
    //    auto& subtriangle_data = *subtriangle_data_;
    //    auto& subtriangle_data_stored = *subtriangle_data_stored_;

    while (true) {
      sum_of_norm_of_errors_ = 0;
      size_t num_triangles = current_triangles_.size();
      if (num_triangles > I_0_vector_.size()) {
        I_0_vector_.resize(num_triangles);
        I_1_vector_.resize(num_triangles);
        error_vector_.resize(num_triangles);
        norm_of_errors_vector_.resize(num_triangles);
        data_.resize(num_triangles);
        data_stored_.resize(num_triangles, false);
        subtriangle_data_.resize(num_triangles);
        subtriangle_data_stored_.resize(num_triangles, false);
      }
      // loop over triangles
      for (size_t ii = 0; ii < num_triangles; ++ii) {
        auto& I_0 = I_0_vector_[ii];
        auto& I_1 = I_1_vector_[ii];
        auto& local_error = error_vector_[ii];
        const auto& triangle = *(current_triangles_[ii]);
        // calculate I_0
        const auto& barycentre = triangle.get_quadrature_point();
        if (!data_stored_[ii]) {
          data_[ii] = data_function(barycentre.position());
          data_stored_[ii] = true;
        }
        I_0 = barycentre.weight() * psi(barycentre.position(), data_[ii]);
        // calculate I_1
        I_1 *= 0;
        for (const auto& subtriangle : triangle.get_subtriangles()) {
          const auto& subbarycentre = subtriangle.get_quadrature_point();
          if (!subtriangle_data_stored_[ii]) {
            subtriangle_data_[ii] = data_function(subbarycentre.position());
            subtriangle_data_stored_[ii] = true;
          }
          I_1 += subbarycentre.weight() * psi(subbarycentre.position(), data_[ii]);
        }
        // calculate E(K)
        local_error = I_1;
        local_error -= I_0;
        norm_of_errors_vector_[ii] = std::abs(local_error) * 4. / 3.;
        sum_of_norm_of_errors_ += norm_of_errors_vector_[ii];
      } // loop over triangles;

      error_ *= 0;
      result_ *= 0;
      for (size_t ii = 0; ii < num_triangles; ++ii) {
        error_ += error_vector_[ii];
        result_ += I_1_vector_[ii];
      }
      const auto error_norm = std::abs(error_);
      const auto result_norm = std::abs(result_);
      if (error_norm < tol_ * result_norm)
        break;

      for (size_t ii = 0; ii < num_triangles; ++ii) {
        if (norm_of_errors_vector_[ii] > gamma_ / num_triangles * sum_of_norm_of_errors_) {
          const auto& subtriangles = current_triangles_[ii]->get_subtriangles();
          // no need to transfer data for the midpoint subtriangle, as it is the same as for
          // the father triangle
          current_triangles_[ii] = &subtriangles[3];
          subtriangle_data_stored_[ii] = false;
          for (size_t jj = 0; jj < 3; ++jj) {
            current_triangles_.push_back(&subtriangles[jj]);
            data_.push_back(subtriangle_data_[ii][jj]);
          }
        }
      }
      // if no triangles where added, refine all triangles as error is evenly distributed
      if (current_triangles_.size() == num_triangles) {
        for (size_t ii = 0; ii < num_triangles; ++ii) {
          const auto& subtriangles = current_triangles_[ii]->get_subtriangles();
          current_triangles_[ii] = &subtriangles[3];
          subtriangle_data_stored_[ii] = false;
          for (size_t jj = 0; jj < 3; ++jj) {
            current_triangles_.push_back(&subtriangles[jj]);
            data_.push_back(subtriangle_data_[ii][jj]);
          }
        }
      }
    } // while(true)
    return result_;
  } // void calculate_integral

  void clear_data()
  {
    data_.clear();
    data_stored_.clear();
    subtriangle_data_.clear();
    subtriangle_data_stored_.clear();
  }

private:
  TriangulationType triangulation_;
  std::vector<RangeType> I_0_vector_, I_1_vector_, error_vector_;
  std::vector<RangeFieldImp> norm_of_errors_vector_;
  RangeType error_, result_;
  RangeFieldImp sum_of_norm_of_errors_;
  double tol_;
  double gamma_;
  std::vector<const TriangleType*> current_triangles_;
  std::vector<DataType> data_;
  std::vector<bool> data_stored_;
  std::vector<Dune::FieldVector<DataType, 3>> subtriangle_data_;
  std::vector<bool> subtriangle_data_stored_;
};

Dune::QuadratureRule<double, 3> get_equally_dist_quad_points_on_spherical_triangle(
    const Dune::XT::Common::FieldVector<Dune::XT::Common::FieldVector<double, 3>, 3>& vertices,
    const size_t num_refinements = 1)
{
  Dune::QuadratureRule<double, 3> ret;
  if (num_refinements == 0) {
    return get_barycentre_rule(vertices);
  } else {
    Dune::FieldVector<Dune::XT::Common::FieldVector<double, 3>, 3> midpoints;
    for (size_t ii = 0; ii < 3; ++ii) {
      midpoints[ii] = vertices[(1 + ii) % 3] + vertices[(2 + ii) % 3];
      midpoints[ii] /= midpoints[ii].two_norm();
    }
    XT::Common::FieldVector<Dune::XT::Common::FieldVector<Dune::XT::Common::FieldVector<double, 3>, 3>, 4> subtriangles;
    subtriangles[0] = Dune::XT::Common::FieldVector<Dune::XT::Common::FieldVector<double, 3>, 3>(
        {vertices[0], midpoints[1], midpoints[2]});
    subtriangles[1] = Dune::XT::Common::FieldVector<Dune::XT::Common::FieldVector<double, 3>, 3>(
        {vertices[1], midpoints[2], midpoints[0]});
    subtriangles[2] = Dune::XT::Common::FieldVector<Dune::XT::Common::FieldVector<double, 3>, 3>(
        {vertices[2], midpoints[0], midpoints[1]});
    subtriangles[3] = Dune::XT::Common::FieldVector<Dune::XT::Common::FieldVector<double, 3>, 3>(
        {midpoints[0], midpoints[1], midpoints[2]});
    for (const auto& subtriangle : subtriangles)
      for (auto&& quad_point : get_equally_dist_quad_points_on_spherical_triangle(subtriangle, num_refinements - 1))
        ret.push_back(quad_point);
    return ret;
  }
} // get_equally_dist_quad_points_on_spherical_triangle(...)

Dune::QuadratureRule<double, 3> get_equally_dist_quad_points_on_poly(const CGALWrapper::Polyhedron_3& poly,
                                                                     const size_t num_refinements = 1)
{
  Dune::QuadratureRule<double, 3> ret;
  // get vertices
  const auto facets_it_end = poly.facets_end();
  for (auto facets_it = poly.facets_begin(); facets_it != facets_it_end; ++facets_it) {
    Dune::XT::Common::FieldVector<Dune::XT::Common::FieldVector<double, 3>, 3> vertices;
    auto halfedge_circ = facets_it->facet_begin();
    auto halfedge_circ_begin = halfedge_circ;
    size_t ii = 0;
    do {
      const auto& point = halfedge_circ->vertex()->point();
      vertices[ii] = Dune::XT::Common::FieldVector<double, 3>({point.x(), point.y(), point.z()});
    } while (++ii, ++halfedge_circ != halfedge_circ_begin);
    for (auto&& quad_point :
         Dune::GDT::Hyperbolic::Problems::get_equally_dist_quad_points_on_spherical_triangle(vertices, num_refinements))
      ret.push_back(quad_point);
  } // iterate over faces
  return ret;
} // get_equally_dist_quad_points_on_poly(...)

template <class DomainType, class PolyhedronType>
std::vector<double> evaluate_spherical_barycentric_coordinates(const DomainType& v, const PolyhedronType& poly)
{
  std::vector<double> ret(poly.size_of_vertices(), 0);
  // walk over facets
  std::vector<typename PolyhedronType::Vertex_const_handle> local_vertices(3);
  const auto facets_it_end = poly.facets_end();
  for (auto facets_it = poly.facets_begin(); facets_it != facets_it_end; ++facets_it) {
    // circulate halfedges around facets
    const auto halfedge_it_begin = facets_it->facet_begin();
    auto halfedge_it = facets_it->facet_begin();
    size_t index = 0;
    do {
      local_vertices[index] = halfedge_it->vertex();
    } while (++index, ++halfedge_it != halfedge_it_begin);
    DomainType barycentric_coords(0);
    bool success = calculate_barycentric_coordinates(v, local_vertices, barycentric_coords);
    if (success) {
      for (size_t ii = 0; ii < 3; ++ii)
        ret[local_vertices[ii]->index] = barycentric_coords[ii];
      break;
    }
  } // facets
  return ret;
} // evaluate_spherical_barycentric_coordinates

template <class RangeType, class DomainType, class PolyhedronType>
RangeType evaluate_linear_partial_basis(const DomainType& v, const PolyhedronType& poly)
{
  RangeType ret(0);
  Dune::FieldMatrix<double, 3, 3> vertices_matrix;
  Dune::FieldMatrix<double, 3, 3> determinant_matrix;
  const auto facets_it_end = poly.facets_end();
  for (auto facets_it = poly.facets_begin(); facets_it != facets_it_end; ++facets_it) {
    // circulate halfedges around facets
    // halfedges are circulated counterclockwise, so if the points is inside the spherical triangle,
    // the coordinate system formed by two adjacent vertices and v is always right-handed, i.e.
    // the triple product is positive
    const auto halfedge_it_begin = facets_it->facet_begin();
    auto halfedge_it = facets_it->facet_begin();
    size_t index = 0;
    do {
      const auto& vertex = halfedge_it->vertex()->point();
      vertices_matrix[index] = Dune::FieldVector<double, 3>({vertex.x(), vertex.y(), vertex.z()});
    } while (++index, ++halfedge_it != halfedge_it_begin);
    bool v_in_this_facet = true;
    // the triple products that need to be positive are the determinants of the matrices (v1, v2, v), (v2, v3, v), (v3,
    // v1, v), where vi is the ith
    // vertex. Swapping two columns changes the sign of det, the matrices used below all have an even number of column
    // swaps
    for (size_t ii = 0; ii < 3; ++ii) {
      determinant_matrix = vertices_matrix;
      determinant_matrix[ii] = v;
      if (XT::Common::FloatCmp::lt(determinant_matrix.determinant(), 0.)) {
        v_in_this_facet = false;
        break;
      }
    }
    if (v_in_this_facet) {
      ret[4 * (facets_it->index)] = 1;
      for (size_t ii = 1; ii < 4; ++ii)
        ret[4 * (facets_it->index) + ii] = v[ii - 1];
      break;
    }
  } // facets
  return ret;
} // evaluate_linear_partial_basis


/** \see class TwoBeams in twobeams.hh */
template <class PointSourceImp,
          class E,
          class D,
          size_t d,
          class R,
          size_t r,
          size_t rC = 1,
          class PointsVectorType = FieldVector<R, r>>
class PointSourceBase : public PlaneSourceBase<PointSourceImp, E, D, d, R, r, rC, PointsVectorType>
{
  typedef PointSourceBase<PointSourceImp, E, D, d, R, r, rC, PointsVectorType> ThisType;
  typedef PlaneSourceBase<PointSourceImp, E, D, d, R, r, rC, PointsVectorType> BaseType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::precision;
  using typename BaseType::DefaultFluxType;
  using typename BaseType::DefaultInitialValueType;
  using typename BaseType::DefaultRHSType;
  using typename BaseType::DefaultBoundaryValueType;

  using typename BaseType::ConfigType;
  using typename BaseType::MatrixType;
  using typename BaseType::RangeType;
  using typename BaseType::RangeFieldType;

  static ConfigType default_grid_config()
  {
    ConfigType grid_config;
    grid_config["type"] = "provider.cube";
    grid_config["lower_left"] = "[-0.5 -0.5 -0.5]";
    grid_config["upper_right"] = "[0.5 0.5 0.5]";
    grid_config["num_elements"] = "[100 100 100]";
    grid_config["overlap_size"] = "[1 1 1 1]";
    return grid_config;
  }

  using BaseType::default_boundary_info_config;
  using BaseType::create_equidistant_points;
  using BaseType::tridiagonal_matrix_inverse;

  template <class... Args>
  PointSourceBase(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  static double CFL()
  {
    return 0.4;
  }

  static double t_end()
  {
    return 0.45;
  }

  // Initial value of the kinetic equation is psi_vac + 1/(8 pi sigma^2) * exp(-|x|^2/(2*sigma^2)).
  // Thus the initial value for the n-th moment is base_integrated_n * (psi_vac + 1/(8 pi sigma^2) *
  // exp(-|x|^2/(2*sigma^2))).
  static ConfigType create_initial_value_config(const ConfigType grid_config = default_grid_config(),
                                                const PointsVectorType& v_points = create_equidistant_points(),
                                                const RangeFieldType psi_vac = 5e-9)
  {
    static const double sigma = 0.03;
    ConfigType initial_value_config;
    initial_value_config["lower_left"] = grid_config["lower_left"];
    initial_value_config["upper_right"] = grid_config["upper_right"];
    initial_value_config["num_elements"] = "[1 1 1]";
    initial_value_config["variable"] = "x";
    initial_value_config["name"] = PointSourceImp::DefaultInitialValueType::static_id();
    const auto basis_integrated = PointSourceImp::basisfunctions_integrated(v_points);
    std::string str = "[";
    for (size_t rr = 0; rr < dimRange; ++rr) {
      if (rr > 0)
        str += " ";
      str += XT::Common::to_string(basis_integrated[rr], precision) + "*(" + XT::Common::to_string(psi_vac, precision)
             + "+" + XT::Common::to_string(1. / (8. * M_PI * sigma * sigma), precision)
             + "*exp((x[0]*x[0]+x[1]*x[1]+x[2]*x[2])/(" + XT::Common::to_string(-2. * sigma * sigma, precision) + ")))";
    } // rr
    str += "]";
    initial_value_config["values.0"] = str;
    initial_value_config["order.0"] = "20";
    return initial_value_config;
  } // ... create_initial_value_config(...)

  // RHS is (G - sigma_t * I)u + Q<b>
  // For this test case (sigma_t = sigma_s + sigma_a),
  // sigma_a = 0, sigma_s = 1, Q = 0
  static ConfigType create_rhs_config(const ConfigType grid_config = default_grid_config(),
                                      const PointsVectorType& v_points = create_equidistant_points())
  {
    ConfigType rhs_config;
    rhs_config["lower_left"] = grid_config["lower_left"];
    rhs_config["upper_right"] = grid_config["upper_right"];
    rhs_config["num_elements"] = "[1 1 1]";
    rhs_config["name"] = PointSourceImp::DefaultRHSType::static_id();
    const RangeType integrated_basis = PointSourceImp::basisfunctions_integrated(v_points);
    const MatrixType M = PointSourceImp::mass_matrix(v_points);
    const MatrixType M_inv = tridiagonal_matrix_inverse(M);
    RangeType c(0);
    M_inv.mtv(integrated_basis, c);
    Dune::FieldMatrix<RangeFieldType, dimRange, dimRange> A(0);
    Dune::FieldVector<RangeFieldType, dimRange> b(0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      for (size_t cc = 0; cc < dimRange; ++cc)
        A[rr][cc] = 0.5 * integrated_basis[rr] * c[cc] - (rr == cc);
    rhs_config["A.0"] = XT::Common::to_string(A);
    rhs_config["b.0"] = XT::Common::to_string(b);
    return rhs_config;
  } // ... create_rhs_config()
}; // class PlaneSourceBase


/** \see class TwoBeams in twobeams.hh */
template <class E, class D, size_t d, class R, size_t momentOrder>
class PointSourcePnLegendre : public PointSourceBase<PointSourcePnLegendre<E, D, d, R, momentOrder>,
                                                     E,
                                                     D,
                                                     d,
                                                     R,
                                                     (momentOrder + 1) * (momentOrder + 1),
                                                     1>
{
  typedef PointSourcePnLegendre<E, D, d, R, momentOrder> ThisType;
  typedef PointSourceBase<PointSourcePnLegendre<E, D, d, R, momentOrder>,
                          E,
                          D,
                          d,
                          R,
                          (momentOrder + 1) * (momentOrder + 1),
                          1>
      BaseType;

public:
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::precision;
  using typename BaseType::FluxRangeType;
  using typename BaseType::DefaultFluxType;
  using typename BaseType::DefaultRHSType;
  using typename BaseType::ConfigType;
  using typename BaseType::MatrixType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  typedef RangeType PointsVectorType;

  static std::string static_id()
  {
    return "PointSourcePnLegendre";
  }

  using BaseType::create_equidistant_points;
  using BaseType::default_grid_config;

  template <class... Args>
  PointSourcePnLegendre(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  // flux matrix A_i,nm = <Omega_i h_n h_m>
  static ConfigType create_flux_config(const PointsVectorType& /*v_points*/ = create_equidistant_points())
  {
    ConfigType flux_config;
    flux_config["type"] = DefaultFluxType::static_id();
    MatrixType A_0(0), A_1(0), A_2(0);
    const auto quadrature = get_lebedev_quadrature(131);
    for (int nn = 0; nn < dimRange; ++nn) {
      for (int mm = 0; mm < dimRange; ++mm) {
        auto l_and_m_nn = get_l_and_m(nn);
        auto l_and_m_mm = get_l_and_m(mm);
        for (const auto& quad_point : quadrature) {
          const auto point = quad_point.position();
          const auto mu = point[0];
          const auto phi = point[1];
          const auto weight = quad_point.weight();
          A_0[nn][mm] += std::sqrt(1. - mu * mu) * std::cos(phi)
                         * evaluate_real_spherical_harmonics(mu, phi, l_and_m_nn.first, l_and_m_nn.second)
                         * evaluate_real_spherical_harmonics(mu, phi, l_and_m_mm.first, l_and_m_mm.second) * weight;
          A_1[nn][mm] += std::sqrt(1. - mu * mu) * std::sin(phi)
                         * evaluate_real_spherical_harmonics(mu, phi, l_and_m_nn.first, l_and_m_nn.second)
                         * evaluate_real_spherical_harmonics(mu, phi, l_and_m_mm.first, l_and_m_mm.second) * weight;
          A_2[nn][mm] += mu * evaluate_real_spherical_harmonics(mu, phi, l_and_m_nn.first, l_and_m_nn.second)
                         * evaluate_real_spherical_harmonics(mu, phi, l_and_m_mm.first, l_and_m_mm.second) * weight;
        }
      }
    }
    flux_config["A.0"] = XT::Common::to_string(A_0, precision);
    flux_config["A.1"] = XT::Common::to_string(A_1, precision);
    flux_config["A.2"] = XT::Common::to_string(A_2, precision);
    flux_config["b"] = Dune::XT::Common::to_string(FluxRangeType(0));
    return flux_config;
  } // ... create_flux_config(...)

  // returns <b>, where b is the basis functions vector
  static RangeType basisfunctions_integrated(const PointsVectorType& /*v_points*/ = create_equidistant_points())
  {
    RangeType ret(0);
    ret[0] = std::sqrt(4. * M_PI);
    return ret;
  }

  // n-th component of RHS is -sigma_t u_n + (sigma_s u_0 + sqrt(4*pi)*Q) delta(n)
  // For this test case (sigma_t = sigma_s + sigma_a),
  // sigma_a = 0, sigma_s = 1, Q = 0
  // Thus, rhs[n] = (delta(n)-1)u[n]
  static ConfigType create_rhs_config(const ConfigType grid_config = default_grid_config(),
                                      const PointsVectorType& /*v_points*/ = create_equidistant_points())
  {
    ConfigType rhs_config;
    rhs_config["lower_left"] = grid_config["lower_left"];
    rhs_config["upper_right"] = grid_config["upper_right"];
    rhs_config["num_elements"] = "[1 1 1]";
    rhs_config["name"] = DefaultRHSType::static_id();
    MatrixType A(0);
    RangeType b(0);
    for (size_t rr = 0; rr < dimRange; ++rr)
      A[rr][rr] = (rr == 0) - 1;
    rhs_config["A.0"] = XT::Common::to_string(A);
    rhs_config["b.0"] = XT::Common::to_string(b);
    return rhs_config;
  } // ... create_rhs_config()
}; // class PointSourcePnLegendre


/** \see class TwoBeams in twobeams.hh */
template <class E, class D, size_t d, class R, size_t num_vertices>
class PointSourcePnHatFunctions
    : public PointSourceBase<PointSourcePnHatFunctions<E, D, d, R, num_vertices>, E, D, d, R, num_vertices, 1>
{
  typedef PointSourcePnHatFunctions<E, D, d, R, num_vertices> ThisType;
  typedef PointSourceBase<PointSourcePnHatFunctions<E, D, d, R, num_vertices>, E, D, d, R, num_vertices, 1> BaseType;

public:
  typedef typename CGALWrapper::Polyhedron_3 Polyhedron_3;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::precision;
  using typename BaseType::FluxRangeType;
  using typename BaseType::RHSType;
  using typename BaseType::FluxType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::DefaultFluxType;
  using typename BaseType::DefaultRHSType;
  typedef typename XT::Functions::GlobalLambdaFunction<E, D, d, R, dimRange> InitialValueFunctionType;
  typedef typename XT::Functions::FunctionCheckerboardFunction<InitialValueFunctionType, E, D, d, R, dimRange>
      DefaultInitialValueType;
  using typename BaseType::DefaultBoundaryValueType;
  using typename BaseType::ConfigType;
  using typename BaseType::MatrixType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::DomainType;
  typedef RangeType PointsVectorType;

  static std::string static_id()
  {
    return "PointSourcePnHatFunctions";
  }

  using BaseType::create_equidistant_points;
  using BaseType::default_grid_config;
  using BaseType::default_boundary_info_config;

  static std::unique_ptr<ThisType> create(const Dune::QuadratureRule<double, 3>& quadrature,
                                          const Polyhedron_3& poly,
                                          const ConfigType config = default_config())
  {
    const std::shared_ptr<const FluxType> flux(DefaultFluxType::create(config.sub("flux")));
    const std::shared_ptr<const RHSType> rhs(DefaultRHSType::create(config.sub("rhs")));
    const std::shared_ptr<const InitialValueType> initial_values(DefaultInitialValueType::create(
        config.sub("initial_values"), "", create_initial_value_lambda(quadrature, poly)));
    const ConfigType grid_config = config.sub("grid");
    const ConfigType boundary_info = config.sub("boundary_info");
    const std::shared_ptr<const BoundaryValueType> boundary_values(
        DefaultBoundaryValueType::create(config.sub("boundary_values")));
    return XT::Common::make_unique<ThisType>(flux, rhs, initial_values, grid_config, boundary_info, boundary_values);
  } // ... create(...)

  static ConfigType default_config(const ConfigType grid_config,
                                   const Dune::QuadratureRule<double, 3>& quadrature,
                                   const Polyhedron_3& poly,
                                   const RangeFieldType psi_vac = 5e-9)
  {
    ConfigType config;
    config.add(grid_config, "grid");
    config.add(default_boundary_info_config(), "boundary_info");
    config.add(create_flux_config(quadrature, poly), "flux");
    config.add(create_rhs_config(grid_config, quadrature, poly), "rhs");
    config.add(create_initial_value_config(grid_config, quadrature, poly, psi_vac), "initial_values");
    config.add(create_boundary_value_config(quadrature, poly, psi_vac), "boundary_values");
    return config;
  } // ... default_config(...)

  template <class... Args>
  PointSourcePnHatFunctions(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  // flux matrix A_i,nm = <Omega_i h_n h_m>
  static ConfigType create_flux_config(const Dune::QuadratureRule<double, 3>& quadrature, const Polyhedron_3& poly)
  {
    ConfigType flux_config;
    flux_config["type"] = DefaultFluxType::static_id();
    MatrixType A_0(0), A_1(0), A_2(0);
    for (const auto& quad_point : quadrature) {
      const auto point = quad_point.position();
      const auto basis_evaluated = evaluate_spherical_barycentric_coordinates(point, poly);
      const auto weight = quad_point.weight();
      for (int nn = 0; nn < dimRange; ++nn) {
        for (int mm = 0; mm < dimRange; ++mm) {
          A_0[nn][mm] += basis_evaluated[nn] * basis_evaluated[mm] * point[0] * weight;
          A_1[nn][mm] += basis_evaluated[nn] * basis_evaluated[mm] * point[1] * weight;
          A_2[nn][mm] += basis_evaluated[nn] * basis_evaluated[mm] * point[2] * weight;
        } // mm
      } // nn
    } // quadrature
    flux_config["A.0"] = XT::Common::to_string(A_0, precision);
    flux_config["A.1"] = XT::Common::to_string(A_1, precision);
    flux_config["A.2"] = XT::Common::to_string(A_2, precision);
    flux_config["b"] = Dune::XT::Common::to_string(FluxRangeType(0));
    return flux_config;
  } // ... create_flux_config(...)

  // returns <b>, where b is the basis functions vector
  static RangeType basisfunctions_integrated(const Dune::QuadratureRule<double, 3>& quadrature,
                                             const Polyhedron_3& poly)
  {
    RangeType ret(0);
    for (const auto& quad_point : quadrature) {
      const auto point = quad_point.position();
      const auto basis_evaluated = evaluate_spherical_barycentric_coordinates(point, poly);
      const auto weight = quad_point.weight();
      for (size_t nn = 0; nn < dimRange; ++nn)
        ret[nn] += basis_evaluated[nn] * weight;
    } // quadrature
    return ret;
  }

  // n-th component of RHS is -sigma_t u_n + sigma_s/(4 pi) <psi><b_n> + Q<b_n>
  // For this test case (sigma_t = sigma_s + sigma_a),
  // sigma_a = 0, sigma_s = 1, Q = 0
  // As <psi> = sum_i u_i, the rhs becomes
  // -u_n + 1/(4 pi) sum_i u_i <b_n>
  static ConfigType create_rhs_config(const ConfigType grid_config,
                                      const Dune::QuadratureRule<double, 3>& quadrature,
                                      const Polyhedron_3& poly)
  {
    const auto basis_integrated = basisfunctions_integrated(quadrature, poly);
    ConfigType rhs_config;
    rhs_config["lower_left"] = grid_config["lower_left"];
    rhs_config["upper_right"] = grid_config["upper_right"];
    rhs_config["num_elements"] = "[1 1 1]";
    rhs_config["name"] = DefaultRHSType::static_id();
    MatrixType A(1);
    A *= 1. / (4 * M_PI);
    RangeType b(0);
    for (size_t nn = 0; nn < dimRange; ++nn) {
      A[nn] *= basis_integrated[nn];
      A[nn][nn] -= 1.;
    }
    rhs_config["A.0"] = XT::Common::to_string(A);
    rhs_config["b.0"] = XT::Common::to_string(b);
    return rhs_config;
  } // ... create_rhs_config()

  // Initial value of the kinetic equation is psi_vac + 1/(8 pi sigma^2) * exp(-|x|^2/(2*sigma^2)).
  // Thus the initial value for the n-th moment is base_integrated_n * (psi_vac + 1/(8 pi sigma^2) *
  // exp(-|x|^2/(2*sigma^2))).
  static ConfigType create_initial_value_config(const ConfigType& grid_config,
                                                const Dune::QuadratureRule<double, 3>& quadrature,
                                                const Polyhedron_3& poly,
                                                const double psi_vac = 5e-9)
  {
    static const double sigma = 0.03;
    ConfigType initial_value_config;
    initial_value_config["lower_left"] = grid_config["lower_left"];
    initial_value_config["upper_right"] = grid_config["upper_right"];
    initial_value_config["num_elements"] = "[1 1 1]";
    initial_value_config["variable"] = "x";
    initial_value_config["name"] = DefaultInitialValueType::static_id();
    const auto basis_integrated = basisfunctions_integrated(quadrature, poly);
    std::string str = "[";
    for (size_t rr = 0; rr < dimRange; ++rr) {
      if (rr > 0)
        str += " ";
      str += XT::Common::to_string(basis_integrated[rr], precision) + "*(" + XT::Common::to_string(psi_vac, precision)
             + "+" + XT::Common::to_string(1. / (8. * M_PI * sigma * sigma), precision)
             + "*exp((x[0]*x[0]+x[1]*x[1]+x[2]*x[2])/(" + XT::Common::to_string(-2. * sigma * sigma, precision) + ")))";
    } // rr
    str += "]";
    initial_value_config["values.0"] = str;
    initial_value_config["order.0"] = "20";
    return initial_value_config;
  } // ... create_initial_value_config(...)

  // Initial value of the kinetic equation is psi_vac + 1/(8 pi sigma^2) * exp(-|x|^2/(2*sigma^2)).
  // Thus the initial value for the n-th moment is base_integrated_n * (psi_vac + 1/(8 pi sigma^2) *
  // exp(-|x|^2/(2*sigma^2))).
  static std::vector<std::function<RangeType(DomainType)>> create_initial_value_lambda(
      const Dune::QuadratureRule<double, 3>& quadrature, const Polyhedron_3& poly, const double psi_vac = 5e-9)
  {
    std::vector<std::function<RangeType(DomainType)>> ret;
    static const double sigma = 0.03;
    const auto basis_integrated = basisfunctions_integrated(quadrature, poly);
    ret.push_back([basis_integrated, psi_vac](const DomainType& x) {
      RangeType result = basis_integrated;
      result *= psi_vac + 1. / (8. * M_PI * sigma * sigma) * std::exp(-1. * x.two_norm() / (2. * sigma * sigma));
      return result;
    });
    return ret;
  } // ... create_initial_value_lambda(...)

  // Boundary value of kinetic equation is psi_vac at both boundaries
  // so n-th component of boundary value has to be \psi_{vac}*base_integrated_n at both boundaries.
  // Modell with constant function.
  static ConfigType create_boundary_value_config(const Dune::QuadratureRule<double, 3>& quadrature,
                                                 const Polyhedron_3& poly,
                                                 const double psi_vac = 5e-9)
  {
    ConfigType boundary_value_config;
    boundary_value_config["type"] = DefaultBoundaryValueType::static_id();
    boundary_value_config["variable"] = "x";
    boundary_value_config["order"] = "1";
    const RangeType integrated_basis = basisfunctions_integrated(quadrature, poly);
    std::string str = "[";
    for (size_t rr = 0; rr < dimRange; ++rr) {
      if (rr > 0)
        str += " ";
      str += XT::Common::to_string(integrated_basis[rr] * psi_vac, precision);
    }
    str += "]";
    boundary_value_config["expression"] = str;
    return boundary_value_config;
  } // ... create_boundary_value_config(...)

}; // class PointSourcePnHatFunctions


/** \see class TwoBeams in twobeams.hh */
template <class E, class D, size_t d, class R, size_t num_faces>
class PointSourcePnPartialMoments
    : public PointSourceBase<PointSourcePnPartialMoments<E, D, d, R, num_faces>, E, D, d, R, 4 * num_faces, 1>
{
  typedef PointSourcePnPartialMoments<E, D, d, R, num_faces> ThisType;
  typedef PointSourceBase<PointSourcePnPartialMoments<E, D, d, R, num_faces>, E, D, d, R, 4 * num_faces, 1> BaseType;

public:
  typedef typename CGALWrapper::Polyhedron_3 Polyhedron_3;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::precision;
  using typename BaseType::FluxRangeType;
  using typename BaseType::RHSType;
  using typename BaseType::FluxType;
  using typename BaseType::InitialValueType;
  using typename BaseType::BoundaryValueType;
  using typename BaseType::DefaultFluxType;
  using typename BaseType::DefaultRHSType;
  typedef typename XT::Functions::GlobalLambdaFunction<E, D, d, R, dimRange> InitialValueFunctionType;
  typedef typename XT::Functions::FunctionCheckerboardFunction<InitialValueFunctionType, E, D, d, R, dimRange>
      DefaultInitialValueType;
  using typename BaseType::DefaultBoundaryValueType;
  using typename BaseType::ConfigType;
  using typename BaseType::MatrixType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::DomainType;
  typedef RangeType PointsVectorType;

  static std::string static_id()
  {
    return "PointSourcePnPartialBasis";
  }

  using BaseType::create_equidistant_points;
  using BaseType::default_grid_config;
  using BaseType::default_boundary_info_config;

  static std::unique_ptr<ThisType> create(const Dune::QuadratureRule<double, 3>& quadrature,
                                          const Polyhedron_3& poly,
                                          const ConfigType config = default_config())
  {
    const std::shared_ptr<const FluxType> flux(DefaultFluxType::create(config.sub("flux")));
    const std::shared_ptr<const RHSType> rhs(DefaultRHSType::create(config.sub("rhs")));
    const std::shared_ptr<const InitialValueType> initial_values(DefaultInitialValueType::create(
        config.sub("initial_values"), "", create_initial_value_lambda(quadrature, poly)));
    const ConfigType grid_config = config.sub("grid");
    const ConfigType boundary_info = config.sub("boundary_info");
    const std::shared_ptr<const BoundaryValueType> boundary_values(
        DefaultBoundaryValueType::create(config.sub("boundary_values")));
    return XT::Common::make_unique<ThisType>(flux, rhs, initial_values, grid_config, boundary_info, boundary_values);
  } // ... create(...)

  static ConfigType default_config(const ConfigType grid_config,
                                   const Dune::QuadratureRule<double, 3>& quadrature,
                                   const Polyhedron_3& poly,
                                   const RangeFieldType psi_vac = 5e-9)
  {
    ConfigType config;
    config.add(grid_config, "grid");
    config.add(default_boundary_info_config(), "boundary_info");
    config.add(create_flux_config(quadrature, poly), "flux");
    config.add(create_rhs_config(grid_config, quadrature, poly), "rhs");
    config.add(create_initial_value_config(grid_config, quadrature, poly, psi_vac), "initial_values");
    config.add(create_boundary_value_config(quadrature, poly, psi_vac), "boundary_values");
    return config;
  } // ... default_config(...)

  template <class... Args>
  PointSourcePnPartialMoments(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  // flux matrix A_i,nm = <Omega_i h_n h_m>
  static ConfigType create_flux_config(const Dune::QuadratureRule<double, 3>& quadrature, const Polyhedron_3& poly)
  {
    ConfigType flux_config;
    flux_config["type"] = DefaultFluxType::static_id();
    MatrixType A_0(0), A_1(0), A_2(0);
    for (const auto& quad_point : quadrature) {
      const auto point = quad_point.position();
      const auto basis_evaluated = evaluate_linear_partial_basis<RangeType, DomainType, Polyhedron_3>(point, poly);
      const auto weight = quad_point.weight();
      for (int nn = 0; nn < dimRange; ++nn) {
        for (int mm = 0; mm < dimRange; ++mm) {
          A_0[nn][mm] += basis_evaluated[nn] * basis_evaluated[mm] * point[0] * weight;
          A_1[nn][mm] += basis_evaluated[nn] * basis_evaluated[mm] * point[1] * weight;
          A_2[nn][mm] += basis_evaluated[nn] * basis_evaluated[mm] * point[2] * weight;
        } // mm
      } // nn
    } // quadrature
    flux_config["A.0"] = XT::Common::to_string(A_0, precision);
    flux_config["A.1"] = XT::Common::to_string(A_1, precision);
    flux_config["A.2"] = XT::Common::to_string(A_2, precision);
    flux_config["b"] = Dune::XT::Common::to_string(FluxRangeType(0));
    return flux_config;
  } // ... create_flux_config(...)

  // returns <b>, where b is the basis functions vector
  static RangeType basisfunctions_integrated(const Dune::QuadratureRule<double, 3>& quadrature,
                                             const Polyhedron_3& poly)
  {
    RangeType ret(0);
    for (const auto& quad_point : quadrature) {
      const auto point = quad_point.position();
      const auto basis_evaluated = evaluate_linear_partial_basis<RangeType, DomainType, Polyhedron_3>(point, poly);
      const auto weight = quad_point.weight();
      for (size_t nn = 0; nn < dimRange; ++nn)
        ret[nn] += basis_evaluated[nn] * weight;
    } // quadrature
    return ret;
  }

  // n-th component of RHS is -sigma_t u_n + sigma_s/(4 pi) <psi><b_n> + Q<b_n>
  // For this test case (sigma_t = sigma_s + sigma_a),
  // sigma_a = 0, sigma_s = 1, Q = 0
  // As <psi> = sum_i u_i, the rhs becomes
  // -u_n + 1/(4 pi) sum_i u_i <b_n>
  static ConfigType create_rhs_config(const ConfigType grid_config,
                                      const Dune::QuadratureRule<double, 3>& quadrature,
                                      const Polyhedron_3& poly)
  {
    const auto basis_integrated = basisfunctions_integrated(quadrature, poly);
    ConfigType rhs_config;
    rhs_config["lower_left"] = grid_config["lower_left"];
    rhs_config["upper_right"] = grid_config["upper_right"];
    rhs_config["num_elements"] = "[1 1 1]";
    rhs_config["name"] = DefaultRHSType::static_id();
    MatrixType A(1);
    A *= 1. / (4 * M_PI);
    RangeType b(0);
    for (size_t nn = 0; nn < dimRange; ++nn) {
      A[nn] *= basis_integrated[nn];
      A[nn][nn] -= 1.;
    }
    rhs_config["A.0"] = XT::Common::to_string(A);
    rhs_config["b.0"] = XT::Common::to_string(b);
    return rhs_config;
  } // ... create_rhs_config()

  // Initial value of the kinetic equation is psi_vac + 1/(8 pi sigma^2) * exp(-|x|^2/(2*sigma^2)).
  // Thus the initial value for the n-th moment is base_integrated_n * (psi_vac + 1/(8 pi sigma^2) *
  // exp(-|x|^2/(2*sigma^2))).
  static ConfigType create_initial_value_config(const ConfigType& grid_config,
                                                const Dune::QuadratureRule<double, 3>& quadrature,
                                                const Polyhedron_3& poly,
                                                const double psi_vac = 5e-9)
  {
    static const double sigma = 0.03;
    ConfigType initial_value_config;
    initial_value_config["lower_left"] = grid_config["lower_left"];
    initial_value_config["upper_right"] = grid_config["upper_right"];
    initial_value_config["num_elements"] = "[1 1 1]";
    initial_value_config["variable"] = "x";
    initial_value_config["name"] = DefaultInitialValueType::static_id();
    const auto basis_integrated = basisfunctions_integrated(quadrature, poly);
    std::string str = "[";
    for (size_t rr = 0; rr < dimRange; ++rr) {
      if (rr > 0)
        str += " ";
      str += XT::Common::to_string(basis_integrated[rr], precision) + "*(" + XT::Common::to_string(psi_vac, precision)
             + "+" + XT::Common::to_string(1. / (8. * M_PI * sigma * sigma), precision)
             + "*exp((x[0]*x[0]+x[1]*x[1]+x[2]*x[2])/(" + XT::Common::to_string(-2. * sigma * sigma, precision) + ")))";
    } // rr
    str += "]";
    initial_value_config["values.0"] = str;
    initial_value_config["order.0"] = "20";
    return initial_value_config;
  } // ... create_initial_value_config(...)

  // Initial value of the kinetic equation is psi_vac + 1/(8 pi sigma^2) * exp(-|x|^2/(2*sigma^2)).
  // Thus the initial value for the n-th moment is base_integrated_n * (psi_vac + 1/(8 pi sigma^2) *
  // exp(-|x|^2/(2*sigma^2))).
  static std::vector<std::function<RangeType(DomainType)>> create_initial_value_lambda(
      const Dune::QuadratureRule<double, 3>& quadrature, const Polyhedron_3& poly, const double psi_vac = 5e-9)
  {
    std::vector<std::function<RangeType(DomainType)>> ret;
    static const double sigma = 0.03;
    const auto basis_integrated = basisfunctions_integrated(quadrature, poly);
    ret.push_back([basis_integrated, psi_vac](const DomainType& x) {
      RangeType result = basis_integrated;
      result *= psi_vac + 1. / (8. * M_PI * sigma * sigma) * std::exp(-1. * x.two_norm() / (2. * sigma * sigma));
      return result;
    });
    return ret;
  } // ... create_initial_value_lambda(...)

  // Boundary value of kinetic equation is psi_vac at both boundaries
  // so n-th component of boundary value has to be \psi_{vac}*base_integrated_n at both boundaries.
  // Modell with constant function.
  static ConfigType create_boundary_value_config(const Dune::QuadratureRule<double, 3>& quadrature,
                                                 const Polyhedron_3& poly,
                                                 const double psi_vac = 5e-9)
  {
    ConfigType boundary_value_config;
    boundary_value_config["type"] = DefaultBoundaryValueType::static_id();
    boundary_value_config["variable"] = "x";
    boundary_value_config["order"] = "1";
    const RangeType integrated_basis = basisfunctions_integrated(quadrature, poly);
    std::string str = "[";
    for (size_t rr = 0; rr < dimRange; ++rr) {
      if (rr > 0)
        str += " ";
      str += XT::Common::to_string(integrated_basis[rr] * psi_vac, precision);
    }
    str += "]";
    boundary_value_config["expression"] = str;
    return boundary_value_config;
  } // ... create_boundary_value_config(...)

}; // class PointSourcePnPartialMoments


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_POINTSOURCE_HH
