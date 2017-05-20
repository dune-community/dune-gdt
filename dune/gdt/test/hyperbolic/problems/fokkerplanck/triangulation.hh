// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_BASISFUNCTIONS_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_BASISFUNCTIONS_HH

#include <memory>
#include <vector>
#include <string>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/l2.hh>
#include <dune/gdt/spaces/cg.hh>

#include <dune/xt/common/string.hh>
#include <dune/xt/functions/affine.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/la/container.hh>

#include "../default.hh"

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

private:
  DomainType position_;
  size_t index_;
}; // class Vertex

template <class RangeFieldImp = double>
class SphericalTriangle
{
  typedef SphericalTriangle<RangeFieldImp> ThisType;

public:
  typedef Vertex<RangeFieldImp, 3> VertexType;
  typedef typename VertexType::DomainType DomainType;
  typedef typename std::vector<std::shared_ptr<VertexType>> TriangulationVerticesVectorType;
  typedef typename Dune::FieldVector<std::shared_ptr<VertexType>, 3> VertexVectorType;
  typedef typename Dune::FieldVector<std::shared_ptr<ThisType>, 4> SubtrianglesVectorType;
  typedef typename Dune::QuadraturePoint<RangeFieldImp, 3> QuadraturePointType;

  SphericalTriangle() = default;
  SphericalTriangle(ThisType&& other) = default;
  ThisType& operator=(ThisType&& other) = default;

  SphericalTriangle(TriangulationVerticesVectorType& triangulation_vertices,
                    const VertexVectorType& vertices,
                    std::atomic<size_t>* current_face_index,
                    std::atomic<size_t>* current_vertex_index,
                    std::mutex* triangulation_vertices_mutex)
    : triangulation_vertices_(triangulation_vertices)
    , vertices_(vertices)
    , current_face_index_(current_face_index)
    , current_vertex_index_(current_vertex_index)
    , triangulation_vertices_mutex_(triangulation_vertices_mutex)
    , index_(current_face_index_++)
    , subtriangle_once_flag_(new std::once_flag)
  {
    calculate_barycentre_rule();
  }

  SphericalTriangle(TriangulationVerticesVectorType& triangulation_vertices,
                    const std::shared_ptr<VertexType> vertex_1,
                    const std::shared_ptr<VertexType> vertex_2,
                    const std::shared_ptr<VertexType> vertex_3,
                    std::atomic<size_t>* current_face_index,
                    std::atomic<size_t>* current_vertex_index,
                    std::mutex* triangulation_vertices_mutex)
    : triangulation_vertices_(triangulation_vertices)
    , vertices_{vertex_1, vertex_2, vertex_3}
    , current_face_index_(current_face_index)
    , current_vertex_index_(current_vertex_index)
    , triangulation_vertices_mutex_(triangulation_vertices_mutex)
    , index_((*current_face_index_)++)
    , subtriangle_once_flag_(new std::once_flag)
  {
    calculate_barycentre_rule();
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

  const QuadraturePointType& quadrature_point() const
  {
    return *quadrature_point_;
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
      auto midpoint_position = vertices_[ii]->position() + vertices_[(1 + ii) % 3]->position();
      midpoint_position /= midpoint_position.two_norm();
      const auto midpoint_iterator =
          std::find_if(triangulation_vertices_.begin(),
                       triangulation_vertices_.end(),
                       [&](const std::shared_ptr<VertexType>& vertex) {
                         return XT::Common::FloatCmp::eq(vertex->position(), midpoint_position);
                       });
      if (midpoint_iterator != triangulation_vertices_.end()) {
        midpoints[ii] = *midpoint_iterator;
      } else {
        triangulation_vertices_.emplace_back(
            std::make_shared<VertexType>(midpoint_position, (*current_vertex_index_)++));
        midpoints[ii] = triangulation_vertices_.back();
      }
    }
    subtriangles_[0] = std::make_shared<ThisType>(triangulation_vertices_,
                                                  vertices_[0],
                                                  midpoints[0],
                                                  midpoints[2],
                                                  current_face_index_,
                                                  current_vertex_index_,
                                                  triangulation_vertices_mutex_);
    subtriangles_[1] = std::make_shared<ThisType>(triangulation_vertices_,
                                                  vertices_[1],
                                                  midpoints[1],
                                                  midpoints[0],
                                                  current_face_index_,
                                                  current_vertex_index_,
                                                  triangulation_vertices_mutex_);
    subtriangles_[2] = std::make_shared<ThisType>(triangulation_vertices_,
                                                  vertices_[2],
                                                  midpoints[2],
                                                  midpoints[1],
                                                  current_face_index_,
                                                  current_vertex_index_,
                                                  triangulation_vertices_mutex_);
    subtriangles_[3] = std::make_shared<ThisType>(triangulation_vertices_,
                                                  midpoints[0],
                                                  midpoints[1],
                                                  midpoints[2],
                                                  current_face_index_,
                                                  current_vertex_index_,
                                                  triangulation_vertices_mutex_);
  } // initialize_subtriangles()

  void calculate_barycentre_rule() const
  {
    const auto ff = (vertices_[0]->position() + vertices_[1]->position() + vertices_[2]->position()) / 3.;
    const auto norm_ff = ff.two_norm();
    const auto bb = ff / norm_ff;
    const auto norm_ff_3 = std::pow(norm_ff, 3);
    const auto vertices_1_minus_0 = vertices_[1]->position() - vertices_[0]->position();
    const auto vertices_2_minus_0 = vertices_[2]->position() - vertices_[0]->position();
    const auto partial_s_gg = vertices_2_minus_0 / norm_ff - ff * (ff * vertices_2_minus_0 / norm_ff_3);
    const auto partial_t_gg = vertices_1_minus_0 / norm_ff - ff * (ff * vertices_1_minus_0 / norm_ff_3);
    const auto weight = Dune::PDELab::crossproduct(partial_s_gg, partial_t_gg).two_norm() / 2.;
    quadrature_point_ = std::make_unique<QuadraturePointType>(bb, weight);
  }

  TriangulationVerticesVectorType& triangulation_vertices_;
  const VertexVectorType vertices_;
  mutable std::unique_ptr<QuadraturePointType> quadrature_point_;
  mutable SubtrianglesVectorType subtriangles_;
  mutable std::atomic<size_t>* current_face_index_;
  mutable std::atomic<size_t>* current_vertex_index_;
  mutable std::mutex* triangulation_vertices_mutex_;
  size_t index_;
  std::unique_ptr<std::once_flag> subtriangle_once_flag_;
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

  SphericalTriangulation(const std::vector<DomainType>& initial_points, size_t num_refinements = 0)
    : current_face_index_(0)
    , current_vertex_index_(0)
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

  void refine(size_t times = 1)
  {
    faces_ = get_subtriangles(times);
  } // void refine(...)

  Dune::QuadratureRule<RangeFieldImp, 3> quadrature_rule(size_t refinements = 0) const
  {
    TriangleVectorType quadrature_faces = get_subtriangles(refinements);
    Dune::QuadratureRule<RangeFieldImp, 3> ret(quadrature_faces.size());
    for (size_t ii = 0; ii < quadrature_faces.size(); ++ii)
      ret[ii] = quadrature_faces[ii]->quadrature_point();
    return ret;
  }

private:
  TriangleVectorType get_subtriangles(size_t refinements = 1) const
  {
    TriangleVectorType subtriangles = faces_;
    while (refinements-- > 0) {
      const size_t old_size = subtriangles.size();
      subtriangles.resize(4. * old_size);
      for (size_t ii = 0; ii < old_size; ++ii) {
        const auto& local_subtriangles = subtriangles[ii]->subtriangles();
        for (size_t jj = 0; jj < 4; ++jj)
          subtriangles[(3 - jj) * old_size + ii] = local_subtriangles[jj];
      } // faces
    } // refinements
    return subtriangles;
  } // ... get_subtriangles(...)

  void calculate_faces(std::vector<DomainType> points0)
  {
    for (const auto& point : points0)
      vertices_.emplace_back(std::make_shared<VertexType>(point, current_vertex_index_++));
    const auto all_vertices = vertices_;
    auto vertices0 = vertices_;
    while (vertices0.size() > 0) {
      const auto v0 = vertices0.back();
      vertices0.pop_back();
      auto vertices1 = vertices0;
      while (vertices1.size() > 0) {
        const auto v1 = vertices1.back();
        vertices1.pop_back();
        for (const auto& v2 : vertices1) {
          // calculate plane equation defined by three points
          const auto v0v1 = v1->position() - v0->position();
          const auto v0v2 = v2->position() - v0->position();
          const auto normal = PDELab::crossproduct(v0v1, v0v2);
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
                    vertices_, v0, v1, v2, &current_face_index_, &current_vertex_index_, &vertices_mutex_));
              else
                faces_.emplace_back(std::make_shared<TriangleType>(
                    vertices_, v0, v2, v1, &current_face_index_, &current_vertex_index_, &vertices_mutex_));
            } // if (is_face)
          } // check if points define a plane
        } // p2
      } // p1
    } // p0
  } // void calculate_faces(...)

  TriangleVectorType faces_;
  VertexVectorType vertices_;
  mutable std::mutex vertices_mutex_;
  std::atomic<size_t> current_face_index_;
  std::atomic<size_t> current_vertex_index_;
}; // class SphericalTriangulation<...>


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_BASISFUNCTIONS_HH
