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

#ifndef DUNE_GDT_MOMENTMODELS_TRIANGULATION_HH
#define DUNE_GDT_MOMENTMODELS_TRIANGULATION_HH

#include <memory>
#include <vector>
#include <string>
#include <mutex>

#include <dune/common/typetraits.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/common/fvector.hh>

namespace Dune {
namespace GDT {


template <class FieldType, size_t dimDomain>
class Vertex
{
public:
  using DomainType = XT::Common::FieldVector<FieldType, dimDomain>;

  Vertex() = default;

  Vertex(const DomainType& pos, const size_t index)
    : position_(pos)
    , index_(index)
  {}

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
  using ThisType = SphericalTriangle;

public:
  using RangeFieldType = RangeFieldImp;
  using VertexType = Vertex<RangeFieldType, 3>;
  using DomainType = typename VertexType::DomainType;
  using TriangulationVerticesVectorType = typename std::vector<std::shared_ptr<VertexType>>;
  using VertexVectorType = typename XT::Common::FieldVector<std::shared_ptr<VertexType>, 3>;
  using SubtrianglesVectorType = typename XT::Common::FieldVector<std::shared_ptr<ThisType>, 4>;
  using QuadraturePointType = QuadraturePoint<RangeFieldType, 3>;
  using QuadratureRuleType = QuadratureRule<RangeFieldType, 3>;
  using FieldVectorType = XT::Common::FieldVector<RangeFieldType, 3>;

  SphericalTriangle() = default;
  SphericalTriangle(ThisType&& other) = default;
  ThisType& operator=(ThisType&& other) = default;

  SphericalTriangle(TriangulationVerticesVectorType& triangulation_vertices,
                    const VertexVectorType& vertices,
                    std::shared_ptr<size_t> current_vertex_index,
                    std::shared_ptr<std::mutex> triangulation_vertices_mutex)
    : triangulation_vertices_(triangulation_vertices)
    , vertices_(vertices)
    , current_vertex_index_(current_vertex_index)
    , triangulation_vertices_mutex_(triangulation_vertices_mutex)
  {}

  SphericalTriangle(TriangulationVerticesVectorType& triangulation_vertices,
                    const std::shared_ptr<VertexType> vertex_1,
                    const std::shared_ptr<VertexType> vertex_2,
                    const std::shared_ptr<VertexType> vertex_3,
                    std::shared_ptr<size_t> current_vertex_index,
                    std::shared_ptr<std::mutex> triangulation_vertices_mutex)
    : triangulation_vertices_(triangulation_vertices)
    , vertices_{vertex_1, vertex_2, vertex_3}
    , current_vertex_index_(current_vertex_index)
    , triangulation_vertices_mutex_(triangulation_vertices_mutex)
  {}

  const SubtrianglesVectorType& subtriangles() const
  {
    return subtriangles_;
  }

  SubtrianglesVectorType& subtriangles()
  {
    return subtriangles_;
  }

  const VertexVectorType& vertices()
  {
    return vertices_;
  }

  DomainType center() const
  {
    DomainType center = (vertices_[0]->position() + vertices_[1]->position() + vertices_[2]->position()) / 3.;
    center /= center.two_norm();
    return center;
  }

  QuadratureRuleType quadrature_rule(const QuadratureRule<RangeFieldType, 2>& reference_quadrature_rule) const
  {
    QuadratureRuleType quadrature_rule;
    quadrature_rule.reserve(reference_quadrature_rule.size());
    const FieldVectorType vertices_1_minus_0 = vertices_[1]->position() - vertices_[0]->position();
    const FieldVectorType vertices_2_minus_0 = vertices_[2]->position() - vertices_[0]->position();
    FieldVectorType ff, partial_s_gg, partial_t_gg;
    for (const auto& quad_point : reference_quadrature_rule) {
      const auto& ref_pos = quad_point.position();
      const auto& ref_weight = quad_point.weight();
      // map point to spherical triangle
      ff = vertices_[0]->position() + ref_pos[0] * vertices_1_minus_0 + ref_pos[1] * vertices_2_minus_0;
      const RangeFieldType norm_ff = ff.two_norm();
      const RangeFieldType norm_ff_3 = std::pow(norm_ff, 3);
      partial_s_gg = vertices_2_minus_0 / norm_ff - ff * ((ff * vertices_2_minus_0) / norm_ff_3);
      partial_t_gg = vertices_1_minus_0 / norm_ff - ff * ((ff * vertices_1_minus_0) / norm_ff_3);
      const RangeFieldType weight = XT::Common::cross_product(partial_s_gg, partial_t_gg).two_norm() * ref_weight;
      quadrature_rule.emplace_back(ff / norm_ff, weight);
    }
    return quadrature_rule;
  }

  void initialize_subtriangles()
  {
    if (!(subtriangles_[0])) {
      std::lock_guard<std::mutex> vertices_lock(*triangulation_vertices_mutex_);
      VertexVectorType midpoints;
      for (size_t ii = 0; ii < 3; ++ii) {
        auto& vertex1 = vertices_[ii];
        auto& vertex2 = vertices_[(1 + ii) % 3];
        DomainType midpoint_position = vertex1->position() + vertex2->position();
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
                                                    current_vertex_index_,
                                                    triangulation_vertices_mutex_);
      subtriangles_[1] = std::make_shared<ThisType>(triangulation_vertices_,
                                                    vertices_[1],
                                                    midpoints[1],
                                                    midpoints[0],
                                                    current_vertex_index_,
                                                    triangulation_vertices_mutex_);
      subtriangles_[2] = std::make_shared<ThisType>(triangulation_vertices_,
                                                    vertices_[2],
                                                    midpoints[2],
                                                    midpoints[1],
                                                    current_vertex_index_,
                                                    triangulation_vertices_mutex_);
      subtriangles_[3] = std::make_shared<ThisType>(triangulation_vertices_,
                                                    midpoints[0],
                                                    midpoints[1],
                                                    midpoints[2],
                                                    current_vertex_index_,
                                                    triangulation_vertices_mutex_);
    }
  } // initialize_subtriangles()

private:
  TriangulationVerticesVectorType& triangulation_vertices_;
  const VertexVectorType vertices_;
  SubtrianglesVectorType subtriangles_;
  std::shared_ptr<size_t> current_vertex_index_;
  std::shared_ptr<std::mutex> triangulation_vertices_mutex_;
}; // class SphericalTriangle<...>

template <class RangeFieldImp = double>
class SphericalTriangulation
{
public:
  using TriangleType = SphericalTriangle<RangeFieldImp>;
  using TriangleVectorType = std::vector<std::shared_ptr<TriangleType>>;
  using VertexVectorType = typename TriangleType::TriangulationVerticesVectorType;
  using VertexType = typename TriangleType::VertexType;
  using DomainType = typename VertexType::DomainType;

  static QuadratureRule<RangeFieldImp, 2> barycentre_rule()
  {
    QuadratureRule<RangeFieldImp, 2> ret;
    ret.push_back(QuadraturePoint<RangeFieldImp, 2>({1. / 3., 1. / 3.}, 0.5));
    return ret;
  }

  SphericalTriangulation() {}

  SphericalTriangulation(size_t num_refinements,
                         const std::vector<DomainType>& initial_points =
                             {{1., 0., 0.}, {-1., 0., 0.}, {0., 1., 0.}, {0., -1., 0.}, {0., 0., 1.}, {0., 0., -1.}})
    : current_vertex_index_(std::make_shared<size_t>(0))
    , all_vertices_mutex_(std::make_shared<std::mutex>())
  {
    calculate_faces(initial_points);
    refine(num_refinements);
    calculate_neighbors();
  }

  // Do not allow copying as the triangles hold references to the vectors in this class.
  SphericalTriangulation(const SphericalTriangulation& other) = delete;
  SphericalTriangulation(SphericalTriangulation&& other) = delete;
  SphericalTriangulation& operator=(const SphericalTriangulation& other) = delete;
  SphericalTriangulation& operator=(SphericalTriangulation&& other) = delete;

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

  const std::vector<std::set<size_t>>& neighbors() const
  {
    return neighbors_;
  }

  const std::set<size_t>& neighbors(const VertexType& vertex) const
  {
    return neighbors_[vertex.index()];
  }

  // get indices of all faces that contain point
  std::vector<size_t> get_face_indices(const DomainType& v) const
  {
    std::vector<size_t> face_indices;
    assert(XT::Common::FloatCmp::eq(v * v, 1.));
    FieldMatrix<RangeFieldImp, 3, 3> vertices_matrix;
    FieldMatrix<RangeFieldImp, 3, 3> determinant_matrix;
    for (size_t kk = 0; kk < faces().size(); ++kk) {
      const auto& face = faces_[kk];
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
        face_indices.push_back(kk);
    } // faces
    assert(face_indices.size());
    return face_indices;
  }

  std::vector<QuadratureRule<RangeFieldImp, 3>>
  quadrature_rules(size_t refinements = 0,
                   const QuadratureRule<RangeFieldImp, 2>& reference_quadrature_rule = barycentre_rule()) const
  {
    std::vector<TriangleVectorType> quadrature_faces = get_subtriangles(refinements);
    std::vector<QuadratureRule<RangeFieldImp, 3>> ret(faces_.size());
    for (size_t jj = 0; jj < faces_.size(); ++jj) {
      for (size_t ii = 0; ii < quadrature_faces[jj].size(); ++ii) {
        const auto quad_rule = quadrature_faces[jj][ii]->quadrature_rule(reference_quadrature_rule);
        ret[jj].insert(ret[jj].end(), quad_rule.begin(), quad_rule.end());
      }
    }
    return ret;
  }

  // This returns a vector, which contains for each face kk a FieldVector<size_t, 3> of the
  // indices of the faces that correspond to face kk when it is reflected in direction ii,
  // i.e. reflected_face_indices()[kk][0] is the index of the face which is the result of
  // reflecting face kk through the y-z Plane.
  const std::vector<XT::Common::FieldVector<size_t, 3>>& reflected_face_indices() const
  {
    return reflected_face_indices_;
  }

private:
  void refine(size_t times = 1)
  {
    set_faces_to_subtriangles(times);
    vertices_ = all_vertices_;
    calculate_reflected_faces();
  } // void refine(...)

  std::vector<TriangleVectorType> get_subtriangles(const size_t refinements = 1) const
  {
    std::vector<TriangleVectorType> ret(faces_.size());
    for (size_t jj = 0; jj < faces_.size(); ++jj) {
      TriangleVectorType subtriangles(1, faces_[jj]);
      size_t refs = refinements;
      while (refs-- > 0)
        get_subtriangles(subtriangles);
      ret[jj] = subtriangles;
    } // jj
    return ret;
  } // ... get_subtriangles(...)

  void get_subtriangles(TriangleVectorType& subtriangles) const
  {
    const size_t old_size = subtriangles.size();
    subtriangles.resize(4. * old_size);
    for (size_t ii = 0; ii < old_size; ++ii) {
      subtriangles[ii]->initialize_subtriangles();
      const auto& local_subtriangles = subtriangles[ii]->subtriangles();
      for (size_t jj = 0; jj < 4; ++jj)
        subtriangles[(3 - jj) * old_size + ii] = local_subtriangles[jj];
    }
  }

  void set_faces_to_subtriangles(size_t refinements = 1)
  {
    while (refinements-- > 0)
      get_subtriangles(faces_);
  }

  void calculate_faces(const std::vector<DomainType>& points0)
  {
    for (const auto& point : points0)
      all_vertices_.emplace_back(std::make_shared<VertexType>(point, (*current_vertex_index_)++));
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
                faces_.emplace_back(std::make_shared<TriangleType>(
                    all_vertices_, v0, v1, v2, current_vertex_index_, all_vertices_mutex_));
              else
                faces_.emplace_back(std::make_shared<TriangleType>(
                    all_vertices_, v0, v2, v1, current_vertex_index_, all_vertices_mutex_));
            } // if (is_face)
          } // check if points define a plane
        } // p2
      } // p1
    } // p0
    calculate_reflected_faces();
  } // void calculate_faces(...)

  void calculate_reflected_faces()
  {
    reflected_face_indices_.resize(faces_.size());
    for (size_t kk = 0; kk < faces_.size(); ++kk) {
      const auto midpoint_rule = faces_[kk]->quadrature_rule(barycentre_rule());
      assert(midpoint_rule.size() == 1);
      for (size_t ii = 0; ii < 3; ++ii) {
        auto midpoint_reflected_in_dir_ii = midpoint_rule[0].position();
        midpoint_reflected_in_dir_ii[ii] = -midpoint_reflected_in_dir_ii[ii];
        auto reflected_index = get_face_indices(midpoint_reflected_in_dir_ii);
        assert(reflected_index.size() == 1);
        reflected_face_indices_[kk][ii] = reflected_index[0];
      }
    }
  }

  void calculate_neighbors()
  {
    neighbors_.resize(vertices_.size());
    for (const auto& face : faces_) {
      for (const auto& vertex1 : face->vertices()) {
        for (const auto& vertex2 : face->vertices()) {
          if (vertex1->index() != vertex2->index()) {
            neighbors_[vertex1->index()].insert(vertex2->index());
            neighbors_[vertex2->index()].insert(vertex1->index());
          }
        }
      }
    }
  }

  TriangleVectorType faces_;
  VertexVectorType vertices_;
  std::vector<std::set<size_t>> neighbors_;
  VertexVectorType all_vertices_;
  std::shared_ptr<size_t> current_vertex_index_;
  std::vector<XT::Common::FieldVector<size_t, 3>> reflected_face_indices_;
  std::shared_ptr<std::mutex> all_vertices_mutex_;
}; // class SphericalTriangulation<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_MOMENTMODELS_TRIANGULATION_HH
