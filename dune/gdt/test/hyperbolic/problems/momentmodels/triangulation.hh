// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_TRIANGULATION_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_TRIANGULATION_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/common/typetraits.hh>

#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/string.hh>

#include <dune/xt/la/container.hh>

#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/affine.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/l2.hh>
#include <dune/gdt/spaces/cg.hh>

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


template <class FieldType, size_t dimDomain>
class Vertex
{
public:
  typedef XT::Common::FieldVector<FieldType, dimDomain> DomainType;

  Vertex() = default;

  Vertex(DomainType pos, size_t index)
    : position_(pos)
    , index_(index)
  {
  }

  const DomainType& position() const
  {
    return position_;
  }

  DomainType& position()
  {
    return position_;
  }

  size_t index() const
  {
    return index_;
  }

  void set_child_with(const std::shared_ptr<Vertex>& parent, std::shared_ptr<Vertex>& child)
  {
    children_.insert({parent->index(), child});
  }

  bool has_child_with(const std::shared_ptr<Vertex>& parent)
  {
    return children_.count(parent->index());
  }

  const std::shared_ptr<Vertex>& child(const std::shared_ptr<Vertex>& parent)
  {
    return children_[parent->index()];
  }

private:
  DomainType position_;
  size_t index_;
  std::map<size_t, std::shared_ptr<Vertex>> children_;
}; // class Vertex

template <class RangeFieldImp = double>
class SphericalTriangle
{
  typedef SphericalTriangle<RangeFieldImp> ThisType;

public:
  typedef RangeFieldImp RangeFieldType;
  typedef Vertex<RangeFieldType, 3> VertexType;
  typedef typename VertexType::DomainType DomainType;
  typedef typename std::vector<std::shared_ptr<VertexType>> TriangulationVerticesVectorType;
  typedef typename Dune::FieldVector<std::shared_ptr<VertexType>, 3> VertexVectorType;
  typedef typename Dune::FieldVector<std::shared_ptr<ThisType>, 4> SubtrianglesVectorType;
  typedef typename Dune::QuadraturePoint<RangeFieldType, 3> QuadraturePointType;
  typedef typename Dune::QuadratureRule<RangeFieldType, 3> QuadratureRuleType;
  typedef XT::Common::FieldVector<RangeFieldType, 3> FieldVectorType;

  SphericalTriangle() = default;
  SphericalTriangle(ThisType&& other) = default;
  ThisType& operator=(ThisType&& other) = default;

  SphericalTriangle(TriangulationVerticesVectorType& triangulation_vertices,
                    const VertexVectorType& vertices,
                    std::atomic<size_t>* current_face_index,
                    std::atomic<size_t>* current_vertex_index,
                    std::mutex* triangulation_vertices_mutex,
                    const Dune::QuadratureRule<RangeFieldType, 2>& reference_quadrature_rule)
    : triangulation_vertices_(triangulation_vertices)
    , vertices_(vertices)
    , current_face_index_(current_face_index)
    , current_vertex_index_(current_vertex_index)
    , triangulation_vertices_mutex_(triangulation_vertices_mutex)
    , index_((*current_face_index_)++)
    , subtriangle_once_flag_(new std::once_flag)
    , reference_quadrature_rule_(reference_quadrature_rule)
  {
    calculate_quadrature_rule();
  }

  SphericalTriangle(TriangulationVerticesVectorType& triangulation_vertices,
                    const std::shared_ptr<VertexType> vertex_1,
                    const std::shared_ptr<VertexType> vertex_2,
                    const std::shared_ptr<VertexType> vertex_3,
                    std::atomic<size_t>* current_face_index,
                    std::atomic<size_t>* current_vertex_index,
                    std::mutex* triangulation_vertices_mutex,
                    const Dune::QuadratureRule<RangeFieldType, 2>& reference_quadrature_rule)
    : triangulation_vertices_(triangulation_vertices)
    , vertices_{vertex_1, vertex_2, vertex_3}
    , current_face_index_(current_face_index)
    , current_vertex_index_(current_vertex_index)
    , triangulation_vertices_mutex_(triangulation_vertices_mutex)
    , index_((*current_face_index_)++)
    , subtriangle_once_flag_(new std::once_flag)
    , reference_quadrature_rule_(reference_quadrature_rule)
  {
    calculate_quadrature_rule();
  }

  const SubtrianglesVectorType& subtriangles() const
  {
    std::call_once(*subtriangle_once_flag_, [&] { initialize_subtriangles(); });
    return subtriangles_;
  }

  SubtrianglesVectorType& subtriangles()
  {
    std::call_once(*subtriangle_once_flag_, [&] { initialize_subtriangles(); });
    return subtriangles_;
  }

  const VertexVectorType& vertices()
  {
    return vertices_;
  }

  //  const QuadraturePointType& quadrature_point() const
  //  {
  //    return *quadrature_point_;
  //  }

  const QuadratureRuleType& quadrature_rule() const
  {
    return *quadrature_rule_;
  }

  size_t index()
  {
    return index_;
  }

private:
  void initialize_subtriangles() const
  {
    std::lock_guard<std::mutex> vertices_lock(*triangulation_vertices_mutex_);
    VertexVectorType midpoints;
    for (size_t ii = 0; ii < 3; ++ii) {
      auto& vertex1 = vertices_[ii];
      auto& vertex2 = vertices_[(1 + ii) % 3];
      auto midpoint_position = vertex1->position() + vertex2->position();
      midpoint_position /= midpoint_position.two_norm();
      if (vertex1->has_child_with(vertex2)) {
        midpoints[ii] = vertex1->child(vertex2);
      } else {
        triangulation_vertices_.emplace_back(
            std::make_shared<VertexType>(midpoint_position, (*current_vertex_index_)++));
        vertex1->set_child_with(vertex2, triangulation_vertices_.back());
        vertex2->set_child_with(vertex1, triangulation_vertices_.back());
        midpoints[ii] = triangulation_vertices_.back();
      }
    } // ii
    subtriangles_[0] = std::make_shared<ThisType>(triangulation_vertices_,
                                                  vertices_[0],
                                                  midpoints[0],
                                                  midpoints[2],
                                                  current_face_index_,
                                                  current_vertex_index_,
                                                  triangulation_vertices_mutex_,
                                                  reference_quadrature_rule_);
    subtriangles_[1] = std::make_shared<ThisType>(triangulation_vertices_,
                                                  vertices_[1],
                                                  midpoints[1],
                                                  midpoints[0],
                                                  current_face_index_,
                                                  current_vertex_index_,
                                                  triangulation_vertices_mutex_,
                                                  reference_quadrature_rule_);
    subtriangles_[2] = std::make_shared<ThisType>(triangulation_vertices_,
                                                  vertices_[2],
                                                  midpoints[2],
                                                  midpoints[1],
                                                  current_face_index_,
                                                  current_vertex_index_,
                                                  triangulation_vertices_mutex_,
                                                  reference_quadrature_rule_);
    subtriangles_[3] = std::make_shared<ThisType>(triangulation_vertices_,
                                                  midpoints[0],
                                                  midpoints[1],
                                                  midpoints[2],
                                                  current_face_index_,
                                                  current_vertex_index_,
                                                  triangulation_vertices_mutex_,
                                                  reference_quadrature_rule_);
  } // initialize_subtriangles()

  void calculate_quadrature_rule() const
  {
    quadrature_rule_ = std::make_unique<QuadratureRuleType>();
    for (const auto& quad_point : reference_quadrature_rule_) {
      const auto& ref_pos = quad_point.position();
      const auto& ref_weight = quad_point.weight();
      // map point to spherical triangle
      const FieldVectorType vertices_1_minus_0 = vertices_[1]->position() - vertices_[0]->position();
      const FieldVectorType vertices_2_minus_0 = vertices_[2]->position() - vertices_[0]->position();
      FieldVectorType ff =
          FieldVectorType(vertices_[0]->position() + ref_pos[0] * vertices_1_minus_0 + ref_pos[1] * vertices_2_minus_0);
      const RangeFieldType norm_ff = ff.two_norm();
      const FieldVectorType pos = ff / norm_ff;
      const RangeFieldType norm_ff_3 = std::pow(norm_ff, 3);
      const FieldVectorType partial_s_gg = vertices_2_minus_0 / norm_ff - ff * (ff * vertices_2_minus_0) / norm_ff_3;
      const FieldVectorType partial_t_gg = vertices_1_minus_0 / norm_ff - ff * (ff * vertices_1_minus_0) / norm_ff_3;
      const RangeFieldType weight = XT::Common::cross_product(partial_s_gg, partial_t_gg).two_norm() * ref_weight;
      quadrature_rule_->emplace_back(Dune::QuadraturePoint<RangeFieldType, 3>(pos, weight));
    }
  }

  TriangulationVerticesVectorType& triangulation_vertices_;
  const VertexVectorType vertices_;
  mutable std::unique_ptr<QuadratureRuleType> quadrature_rule_;
  mutable SubtrianglesVectorType subtriangles_;
  mutable std::atomic<size_t>* current_face_index_;
  mutable std::atomic<size_t>* current_vertex_index_;
  mutable std::mutex* triangulation_vertices_mutex_;
  size_t index_;
  std::unique_ptr<std::once_flag> subtriangle_once_flag_;
  const Dune::QuadratureRule<RangeFieldType, 2>& reference_quadrature_rule_;
}; // class SphericalTriangle<...>

template <class RangeFieldImp = double>
class SphericalTriangulation
{
public:
  typedef SphericalTriangle<RangeFieldImp> TriangleType;
  typedef std::vector<std::shared_ptr<TriangleType>> TriangleVectorType;
  typedef typename TriangleType::TriangulationVerticesVectorType VertexVectorType;
  typedef typename TriangleType::VertexType VertexType;
  typedef typename VertexType::DomainType DomainType;

  static QuadratureRule<RangeFieldImp, 2> barycentre_rule()
  {
    Dune::QuadratureRule<RangeFieldImp, 2> ret;
    ret.push_back(Dune::QuadraturePoint<RangeFieldImp, 2>({1. / 3., 1. / 3.}, 0.5));
    return ret;
  }

  SphericalTriangulation(const std::vector<DomainType>& initial_points,
                         size_t num_refinements = 0,
                         const QuadratureRule<RangeFieldImp, 2>& reference_quadrature_rule = barycentre_rule())
    : current_face_index_(0)
    , current_vertex_index_(0)
    , reference_quadrature_rule_(reference_quadrature_rule)
  {
    calculate_faces(initial_points);
    refine(num_refinements);
  }

  const TriangleVectorType& faces() const
  {
    return faces_;
  }

  TriangleVectorType& faces()
  {
    return faces_;
  }

  const VertexVectorType& vertices() const
  {
    return vertices_;
  }

  VertexVectorType& vertices()
  {
    return vertices_;
  }

  // get indices of all faces that contain point
  // indices are unordered
  std::vector<size_t> get_face_indices(const DomainType& v) const
  {
    std::vector<size_t> face_indices;
    DXT_ASSERT(XT::Common::FloatCmp::eq(v * v, 1.));
    FieldMatrix<RangeFieldImp, 3, 3> vertices_matrix;
    FieldMatrix<RangeFieldImp, 3, 3> determinant_matrix;
    for (const auto& face : faces()) {
      bool v_in_this_face = false;
      bool second_check = true;
      const auto& vertices = face->vertices();
      for (size_t ii = 0; ii < 3; ++ii) {
        // if v is not on the same octant of the sphere as the vertices, return false
        // assumes the triangulation is fine enough that vertices[ii]*vertices[jj] >= 0 for all triangles
        const auto scalar_prod = v * vertices[ii]->position();
        if (XT::Common::FloatCmp::lt(scalar_prod, 0.)) {
          second_check = false;
          break;
        } else if (XT::Common::FloatCmp::eq(scalar_prod, 1.)) {
          v_in_this_face = true;
          second_check = false;
          break;
        }
        vertices_matrix[ii] = vertices[ii]->position();
      } // ii

      if (second_check) {
        // Vertices are ordered counterclockwise, so if the point is inside the spherical triangle,
        // the coordinate system formed by two adjacent vertices and v is always right-handed, i.e.
        // the triple product is positive.
        // The triple products that need to be positive are the determinants of the matrices (v1, v2, v), (v2, v3, v),
        // (v3, v1, v), where vi is the ith vertex. Swapping two columns changes the sign of det, the matrices used
        // below all have an even number of column swaps.
        // The determinant is 0 iff v is on the same plane as the two vertices. Then, to be on the edge of the
        // current face, v has be in between the two vertices, i.e.  v = x * v1 + y * v2 with x, y >= 0. This is
        // equivalent to v * v1 >= (v*v2)v2 * v1 && v * v2 >= (v*v1)v1 * v2 (the projection of v to v1 has to be
        // greater than the projection of its v2-projection to v1).
        v_in_this_face = true;
        for (size_t ii = 0; ii < 3; ++ii) {
          determinant_matrix = vertices_matrix;
          determinant_matrix[ii] = v;
          auto det = determinant_matrix.determinant();
          if (XT::Common::FloatCmp::eq(det, 0.)) {
            const auto& v1 = vertices_matrix[(ii + 1) % 3];
            const auto& v2 = vertices_matrix[(ii + 2) % 3];
            v_in_this_face = v * v1 > (v * v2) * (v2 * v1) && v * v2 > (v * v1) * (v1 * v2) ? true : false;
            break;
          } else if (det < 0.) {
            v_in_this_face = false;
            break;
          }
        } // ii
      } // if (second_check)
      if (v_in_this_face)
        face_indices.push_back(face->index());
    } // faces
    DXT_ASSERT(face_indices.size());
    return face_indices;
  }

  // 3D quadrature on sphere (from http://www.unizar.es/galdeano/actas_pau/PDFVIII/pp61-69.pdf)
  Dune::QuadratureRule<RangeFieldImp, 3> quadrature_rule(size_t refinements = 0) const
  {
    TriangleVectorType quadrature_faces = get_subtriangles(refinements);
    Dune::QuadratureRule<RangeFieldImp, 3> ret;
    for (size_t ii = 0; ii < quadrature_faces.size(); ++ii) {
      const auto& quad_rule = quadrature_faces[ii]->quadrature_rule();
      for (const auto& quad_point : quad_rule)
        ret.emplace_back(quad_point);
    }
    return ret;
  }

private:
  void refine(size_t times = 1)
  {
    set_faces_to_subtriangles(times);
    vertices_ = all_vertices_;
  } // void refine(...)

  TriangleVectorType get_subtriangles(size_t refinements = 1) const
  {
    TriangleVectorType subtriangles = faces_;
    while (refinements-- > 0)
      get_subtriangles(subtriangles);
    return subtriangles;
  } // ... get_subtriangles(...)

  void get_subtriangles(TriangleVectorType& subtriangles) const
  {
    const size_t old_size = subtriangles.size();
    subtriangles.resize(4. * old_size);
    for (size_t ii = 0; ii < old_size; ++ii) {
      const auto& local_subtriangles = subtriangles[ii]->subtriangles();
      for (size_t jj = 0; jj < 4; ++jj)
        subtriangles[(3 - jj) * old_size + ii] = local_subtriangles[jj];
    }
  }

  void set_faces_to_subtriangles(size_t refinements = 1)
  {
    while (refinements-- > 0) {
      current_face_index_ = 0;
      get_subtriangles(faces_);
    }
  }

  void calculate_faces(const std::vector<DomainType>& points0)
  {
    for (const auto& point : points0)
      all_vertices_.emplace_back(std::make_shared<VertexType>(point, current_vertex_index_++));
    vertices_ = all_vertices_;
    const auto all_vertices = all_vertices_;
    auto vertices0 = all_vertices_;
    while (vertices0.size() > 0) {
      const auto v0 = vertices0.back();
      vertices0.pop_back();
      auto vertices1 = vertices0;
      while (vertices1.size() > 0) {
        const auto v1 = vertices1.back();
        vertices1.pop_back();
        for (const auto& v2 : vertices1) {
          // calculate plane equation defined by three points
          const DomainType v0v1 = v1->position() - v0->position();
          const DomainType v0v2 = v2->position() - v0->position();
          const DomainType normal = XT::Common::cross_product(v0v1, v0v2);
          if (XT::Common::FloatCmp::ne(normal.two_norm2(), 0.)) {
            bool is_face = true;
            double max_value = std::numeric_limits<double>::lowest();
            double min_value = std::numeric_limits<double>::max();
            for (const auto& v3 : all_vertices) {
              const auto v0v3 = v3->position() - v0->position();
              const auto value = normal * v0v3;
              max_value = std::max(max_value, value);
              min_value = std::min(min_value, value);
              if (XT::Common::FloatCmp::lt(min_value * max_value, 0.)) {
                is_face = false;
                break;
              }
            } // p3
            if (is_face) {
              // if max_value is <= 0, all values are less or equal zero,
              // i.e the normal points outwards and thus p0, p1, p2 are oriented counterclockwise, which is what we want
              if (XT::Common::FloatCmp::le(max_value, 0.))
                faces_.emplace_back(std::make_shared<TriangleType>(all_vertices_,
                                                                   v0,
                                                                   v1,
                                                                   v2,
                                                                   &current_face_index_,
                                                                   &current_vertex_index_,
                                                                   &all_vertices_mutex_,
                                                                   reference_quadrature_rule_));
              else
                faces_.emplace_back(std::make_shared<TriangleType>(all_vertices_,
                                                                   v0,
                                                                   v2,
                                                                   v1,
                                                                   &current_face_index_,
                                                                   &current_vertex_index_,
                                                                   &all_vertices_mutex_,
                                                                   reference_quadrature_rule_));
            } // if (is_face)
          } // check if points define a plane
        } // p2
      } // p1
    } // p0
  } // void calculate_faces(...)

  TriangleVectorType faces_;
  VertexVectorType vertices_;
  VertexVectorType all_vertices_;
  mutable std::mutex all_vertices_mutex_;
  std::atomic<size_t> current_face_index_;
  std::atomic<size_t> current_vertex_index_;
  const QuadratureRule<RangeFieldImp, 2> reference_quadrature_rule_;
}; // class SphericalTriangulation<...>


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_TRIANGULATION_HH
