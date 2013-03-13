#ifndef DUNE_DETIALED_DISCRETIZATIONS_ASSEMBLER_MULTISCALE_SYSTEM_HH
#define DUNE_DETIALED_DISCRETIZATIONS_ASSEMBLER_MULTISCALE_SYSTEM_HH

#include <memory>
#include <vector>

#include <dune/grid/multiscale/default.hh>

#include <dune/detailed/discretizations/assembler/multiscale/coupling.hh>

#if HAVE_DUNE_DETAILED_SOLVERS
#include <dune/detailed/solvers/stationary/linear/elliptic/continuousgalerkin/dune-detailed-discretizations.hh>
#endif // HAVE_DUNE_DETAILED_SOLVERS

namespace Dune {

namespace Detailed {

namespace Discretizations {

namespace Assembler {

namespace Multiscale {

namespace System {

template <class MsGridImp, class LocalSolverImp>
class Primal;

#if HAVE_DUNE_DETAILED_SOLVERS
template <class GridImp, class ModelImp, int polynomialOrder>
class Primal<Dune::grid::Multiscale::Default<GridImp>,
             Dune::Detailed::Solvers::Stationary::Linear::Elliptic::ContinuousGalerkin::
                 DuneDetailedDiscretizations<ModelImp,
                                             typename Dune::grid::Multiscale::Default<GridImp>::LocalGridPartType,
                                             polynomialOrder>>
{
public:
  typedef Dune::grid::Multiscale::Default<GridImp> MsGridType;

  typedef typename MsGridType::LocalGridPartType LocalGridPartType;

  typedef typename Dune::Detailed::Solvers::Stationary::Linear::Elliptic::ContinuousGalerkin::
      DuneDetailedDiscretizations<ModelImp, LocalGridPartType, polynomialOrder> LocalSolverType;

  Primal(const MsGridType& msGrid, const std::vector<std::shared_ptr<LocalSolverType>> localSolvers)
    : msGrid_(msGrid)
    , localSolvers_(localSolvers)
  {
  }

  template <class CouplingMatrixType>
  void assembleCoupling(std::vector<std::map<unsigned int, CouplingMatrixType>>& couplingMatricesMaps) const
  {
    // loop over all subdomains
    for (unsigned int subdomain = 0; subdomain < msGrid_.size(); ++subdomain) {
      // get subdomain matrices map
      std::map<unsigned int, CouplingMatrixType>& subdomainCouplingMatricesMap = couplingMatricesMaps[subdomain];
      // loop over all neighbors
      typedef typename MsGridType::NeighboringSubdomainsSetType NeighboringSubdomainsSetType;
      const std::shared_ptr<const NeighboringSubdomainsSetType> neighbors = msGrid_.neighborsOf(subdomain);
      for (const unsigned int neighbor : *neighbors) {
        // ensure primal assembling
        if (subdomain < neighbor) {
          // get coupling grid part
          typedef typename MsGridType::CouplingGridPartType CouplingGridPartType;
          const std::shared_ptr<const CouplingGridPartType> couplingGridPart =
              msGrid_.couplingGridPart(subdomain, neighbor);
          // get neighbor matrices map
          std::map<unsigned int, CouplingMatrixType>& neighborCouplingMatricesMap = couplingMatricesMaps[neighbor];
          // get subdomain/neighbor coupling matrix
          typename std::map<unsigned int, CouplingMatrixType>::iterator subdomainNeighborMatrixResult =
              subdomainCouplingMatricesMap.find(neighbor);
          assert(subdomainNeighborMatrixResult != subdomainCouplingMatricesMap.end());
          CouplingMatrixType& subdomainNeighborMatrix = subdomainNeighborMatrixResult->second;
          // get neighbor/subdomain coupling matrix
          typename std::map<unsigned int, CouplingMatrixType>::iterator neighborSubdomainMatrixResult =
              neighborCouplingMatricesMap.find(subdomain);
          assert(neighborSubdomainMatrixResult != neighborCouplingMatricesMap.end());
          CouplingMatrixType& neighborSubdomainMatrix = neighborSubdomainMatrixResult->second;
          // create coupling assembler
          typedef typename Dune::Detailed::Discretizations::Assembler::Multiscale::Coupling::
              Primal<CouplingGridPartType,
                     typename LocalSolverType::AnsatzSpaceType,
                     typename LocalSolverType::TestSpaceType,
                     typename LocalSolverType::AnsatzSpaceType,
                     typename LocalSolverType::TestSpaceType> CouplingAssemblerType;
          const CouplingAssemblerType couplingAssembler(*couplingGridPart,
                                                        localSolvers_[subdomain]->ansatzSpace(),
                                                        localSolvers_[subdomain]->testSpace(),
                                                        localSolvers_[neighbor]->ansatzSpace(),
                                                        localSolvers_[neighbor]->testSpace());
        } // ensure primal assembling
      } // loop over all neighbors
    } // loop over all subdomains
  } // void assembleCoupling()

private:
  const MsGridType& msGrid_;
  const std::vector<std::shared_ptr<LocalSolverType>> localSolvers_;
}; // class Primal
#endif // HAVE_DUNE_DETAILED_SOLVERS

} // namespace System

} // namespace Multiscale

} // namespace Assembler

} // namespace Discretizations

} // namespace Detailed

} // namespace Dune

#endif // DUNE_DETIALED_DISCRETIZATIONS_ASSEMBLER_MULTISCALE_SYSTEM_HH
