#ifndef DUNE_GDT_ASSEMBLER_COUPLING_HH
#define DUNE_GDT_ASSEMBLER_COUPLING_HH

#include <dune/xt/common/parallel/threadstorage.hh>

#include <dune/grid/common/rangegenerators.hh>

#include "wrapper.hh"

namespace Dune {
namespace GDT {


/**
 * \todo merge with SystemAssembler
 */
template <class LocalSpaceType, class CouplingGridPartType>
class CouplingAssembler
{
  typedef CouplingAssembler<LocalSpaceType, CouplingGridPartType> ThisType;

public:
  typedef LocalSpaceType TestSpaceType;
  typedef LocalSpaceType AnsatzSpaceType;
  typedef CouplingGridPartType GridLayerType;

  CouplingAssembler(const CouplingGridPartType& coupling_grid_part,
                    const LocalSpaceType& inner_test_space,
                    const LocalSpaceType& inner_ansatz_space,
                    const LocalSpaceType& outer_test_space,
                    const LocalSpaceType& outer_ansatz_space)
    : coupling_grid_part_(coupling_grid_part)
    , inner_test_space_(inner_test_space)
    , inner_ansatz_space_(inner_ansatz_space)
    , outer_test_space_(outer_test_space)
    , outer_ansatz_space_(outer_ansatz_space)
  {
  }

  template <class V, class M, class R>
  ThisType& append(const LocalCouplingTwoFormAssembler<V>& local_assembler,
                   XT::LA::MatrixInterface<M, R>& in_in_matrix,
                   XT::LA::MatrixInterface<M, R>& out_out_matrix,
                   XT::LA::MatrixInterface<M, R>& in_out_matrix,
                   XT::LA::MatrixInterface<M, R>& out_in_matrix)
  {
    assert(in_in_matrix.rows() == inner_test_space_->mapper().size());
    assert(in_in_matrix.cols() == inner_ansatz_space_->mapper().size());
    assert(out_out_matrix.rows() == outer_test_space_->mapper().size());
    assert(out_out_matrix.cols() == outer_ansatz_space_->mapper().size());
    assert(in_out_matrix.rows() == inner_test_space_->mapper().size());
    assert(in_out_matrix.cols() == outer_ansatz_space_->mapper().size());
    assert(out_in_matrix.rows() == outer_test_space_->mapper().size());
    assert(out_in_matrix.cols() == inner_ansatz_space_->mapper().size());
    typedef internal::LocalCouplingTwoFormMatrixAssemblerWrapper<ThisType,
                                                                 LocalCouplingTwoFormAssembler<V>,
                                                                 typename M::derived_type>
        WrapperType;
    codim1_functors_.emplace_back(new WrapperType(inner_test_space_,
                                                  inner_ansatz_space_,
                                                  outer_test_space_,
                                                  outer_ansatz_space_,
                                                  new XT::Grid::ApplyOn::AllIntersections<CouplingGridPartType>(),
                                                  local_assembler,
                                                  in_in_matrix.as_imp(),
                                                  out_out_matrix.as_imp(),
                                                  in_out_matrix.as_imp(),
                                                  out_in_matrix.as_imp()));
    return *this;
  } // ... append(...)

  void assemble()
  {
    for (auto&& entity : elements(coupling_grid_part_)) {
      // only walk the intersections, if there are codim1 functors present
      if (codim1_functors_.size() > 0) {
        // walk the intersections
        // not using intersections(...) on purpose here, would break DD::SubdomainBoundaryGridPart
        const auto intersection_it_end = coupling_grid_part_.iend(entity);
        for (auto intersection_it = coupling_grid_part_.ibegin(entity); intersection_it != intersection_it_end;
             ++intersection_it) {
          const auto& intersection = *intersection_it;
          // apply codim1 functors
          const auto neighbor = intersection.outside();
          for (const auto& functor : codim1_functors_)
            if (functor->apply_on(coupling_grid_part_, intersection))
              functor->apply_local(intersection, entity, neighbor);
        } // walk the intersections
      } // only walk the intersections, if there are codim1 functors present
    }
  } // ... assemble(...)

private:
  const CouplingGridPartType& coupling_grid_part_;
  const Dune::XT::Common::PerThreadValue<const LocalSpaceType> inner_test_space_;
  const Dune::XT::Common::PerThreadValue<const LocalSpaceType> inner_ansatz_space_;
  const Dune::XT::Common::PerThreadValue<const LocalSpaceType> outer_test_space_;
  const Dune::XT::Common::PerThreadValue<const LocalSpaceType> outer_ansatz_space_;

  std::vector<std::unique_ptr<XT::Grid::internal::Codim1Object<CouplingGridPartType>>> codim1_functors_;
}; // class CouplingAssembler


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_COUPLING_HH
