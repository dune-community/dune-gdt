// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_ASSEMBLER_BOUNDARY_HH
#define DUNE_GDT_ASSEMBLER_BOUNDARY_HH

#warning This header is deprecated, #include <dune/gdt/assembler/system.hh> instead (03.05.2017)!

#include <dune/common/deprecated.hh>

#include <dune/xt/common/parallel/threadstorage.hh>

#include "wrapper.hh"

namespace Dune {
namespace GDT {


/**
 * \todo merge with SystemAssembler
 */
template <class LocalSpaceType, class BoundaryGridPartType>
class DUNE_DEPRECATED_MSG("Use SystemAssembler instead (03.05.2017)!") BoundaryAssembler
{
  typedef BoundaryAssembler<LocalSpaceType, BoundaryGridPartType> ThisType;

public:
  typedef LocalSpaceType TestSpaceType;
  typedef LocalSpaceType AnsatzSpaceType;
  typedef BoundaryGridPartType GridLayerType;

  BoundaryAssembler(const BoundaryGridPartType& boundary_grid_part,
                    const LocalSpaceType& test_space,
                    const LocalSpaceType& ansatz_space)
    : boundary_grid_part_(boundary_grid_part)
    , test_space_(test_space)
    , ansatz_space_(ansatz_space)
  {
  }

  template <class V, class M, class R>
  ThisType& append(const LocalBoundaryTwoFormAssembler<V>& local_assembler,
                   XT::LA::MatrixInterface<M, R>& matrix,
                   const XT::Grid::ApplyOn::WhichIntersection<BoundaryGridPartType>* where =
                       new XT::Grid::ApplyOn::AllIntersections<BoundaryGridPartType>())
  {
    assert(matrix.rows() == test_space_->mapper().size());
    assert(matrix.cols() == ansatz_space_->mapper().size());
    typedef internal::LocalBoundaryTwoFormMatrixAssemblerWrapper<ThisType,
                                                                 LocalBoundaryTwoFormAssembler<V>,
                                                                 typename M::derived_type>
        WrapperType;
    codim1_functors_.emplace_back(new WrapperType(test_space_, ansatz_space_, where, local_assembler, matrix.as_imp()));
    return *this;
  } // ... append(...)

  template <class L, class V, class R>
  ThisType& append(const LocalFaceFunctionalAssembler<L>& local_assembler,
                   XT::LA::VectorInterface<V, R>& vector,
                   const XT::Grid::ApplyOn::WhichIntersection<BoundaryGridPartType>* where =
                       new XT::Grid::ApplyOn::AllIntersections<BoundaryGridPartType>())
  {
    assert(vector.size() == test_space_->mapper().size());
    typedef internal::LocalFaceFunctionalVectorAssemblerWrapper<ThisType,
                                                                LocalFaceFunctionalAssembler<L>,
                                                                typename V::derived_type>
        WrapperType;
    codim1_functors_.emplace_back(new WrapperType(test_space_, where, local_assembler, vector.as_imp()));
    return *this;
  } // ... append(...)

  void assemble()
  {
    for (auto&& entity : elements(boundary_grid_part_)) {
      // only walk the intersections, if there are codim1 functors present
      if (codim1_functors_.size() > 0) {
        // walk the intersections
        // not using intersections(...) on purpose here, would break DD::SubdomainBoundaryGridPart
        const auto intersection_it_end = boundary_grid_part_.iend(entity);
        for (auto intersection_it = boundary_grid_part_.ibegin(entity); intersection_it != intersection_it_end;
             ++intersection_it) {
          const auto& intersection = *intersection_it;

          // apply codim1 functors
          if (intersection.neighbor()) {
            const auto neighbor = intersection.outside();
            for (const auto& functor : codim1_functors_)
              if (functor->apply_on(boundary_grid_part_, intersection))
                functor->apply_local(intersection, entity, neighbor);
          } else
            for (const auto& functor : codim1_functors_)
              if (functor->apply_on(boundary_grid_part_, intersection))
                functor->apply_local(intersection, entity, entity);
        } // walk the intersections
      } // only walk the intersections, if there are codim1 functors present
    }
  } // ... assemble(...)

private:
  const BoundaryGridPartType& boundary_grid_part_;
  const Dune::XT::Common::PerThreadValue<const LocalSpaceType> test_space_;
  const Dune::XT::Common::PerThreadValue<const LocalSpaceType> ansatz_space_;

  std::vector<std::unique_ptr<XT::Grid::internal::Codim1Object<BoundaryGridPartType>>> codim1_functors_;
}; // class BoundaryAssembler


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_BOUNDARY_HH
