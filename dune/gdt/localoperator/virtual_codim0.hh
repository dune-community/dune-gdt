// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LOCALOPERATOR_VIRTUAL_INTEGRAL_HH
#define DUNE_GDT_LOCALOPERATOR_VIRTUAL_INTEGRAL_HH

#include <vector>
#include <utility>
#include <type_traits>
#include <limits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/densematrix.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/stuff/functions/interfaces.hh>

#include "../localevaluation/interface.hh"
#include "interface.hh"

#include <dune/patches/finiteelementmap/refinedq1fem.hh>
#include <dune/patches/common/utility.hh>
#include <dune/patches/localfunction/refinedq1.hh>

namespace Dune {
namespace GDT {
namespace LocalOperator {


// forward, to be used in the traits
template< class BinaryEvaluationImp >
class VirtualRefinedCodim0Integral;


namespace internal {


template< class BinaryEvaluationImp >
class VirtualRefinedCodim0IntegralTraits
{
  static_assert(std::is_base_of< LocalEvaluation::Codim0Interface< typename BinaryEvaluationImp::Traits, 2 >,
                                 BinaryEvaluationImp >::value,
                "BinaryEvaluationImp has to be derived from LocalEvaluation::Codim0Interface< ..., 2 >!");
public:
  typedef VirtualRefinedCodim0Integral< BinaryEvaluationImp > derived_type;
};


} // namespace internal


template< class BinaryEvaluationType >
class VirtualRefinedCodim0Integral
  : public LocalOperator::Codim0Interface< internal::VirtualRefinedCodim0IntegralTraits< BinaryEvaluationType > >
{
  static const size_t numTmpObjectsRequired_ = 1;
public:
  typedef internal::VirtualRefinedCodim0IntegralTraits< BinaryEvaluationType > Traits;

  template< class... Args >
  explicit VirtualRefinedCodim0Integral(Args&& ...args)
    : integrand_(std::forward< Args >(args)...)
    , over_integrate_(0)
  {}

  template< class... Args >
  explicit VirtualRefinedCodim0Integral(const int over_integrate, Args&& ...args)
    : integrand_(std::forward< Args >(args)...)
    , over_integrate_(boost::numeric_cast< size_t >(over_integrate))
  {}

  template< class... Args >
  explicit VirtualRefinedCodim0Integral(const size_t over_integrate, Args&& ...args)
    : integrand_(std::forward< Args >(args)...)
    , over_integrate_(over_integrate)
  {}

  size_t numTmpObjectsRequired() const
  {
    return numTmpObjectsRequired_;
  }

  template< class E, class D, size_t d, class R, size_t rT, size_t rCT, size_t rA, size_t rCA >
  void apply(const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
             const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase,
             Dune::DynamicMatrix< R >& ret,
             std::vector< Dune::DynamicMatrix< R > >& tmpLocalMatrices) const
  {
    const auto& original_entity = ansatzBase.entity();

    {
      using Patch = Dune::Patches::Cube::Unconnected::Patch<R, d>;

      auto& gfs = testBase.backend().gfs();
      typedef typename std::remove_const<decltype(gfs)>::type GridView;
      using Factory = Dune::Patches::GridViewPatchFactory<GridView>;

      Factory factory(gfs.gridView());

      using LFS = Dune::PDELab::LocalFunctionSpace<typename std::remove_const<decltype(gfs)>::type>;
      LFS lfs(gfs);

      using LFSCache = Dune::PDELab::LFSIndexCache<LFS>;
      LFSCache lfsCache(lfs);

      //    using LocalView = typename V::template LocalView<LFSCache>;
      //    LocalView localView(v);

      auto patchp = factory.template create<Patch>(makeIteratorRange(&original_entity, &original_entity + 1));
      unsigned level = 0;
      for (auto subdiv = lfs.finiteElement().subDivisions(); subdiv > 1; subdiv /= 2)
        ++level;
      if (1u << level != lfs.finiteElement().subDivisions())
        DUNE_THROW(Dune::NotImplemented, "Refinements with non-power-of-2 subdivisions");
      auto pv = patchp->levelGridView(level);
      auto is = pv.indexSet();

      const auto original_localFunctions = integrand_.localFunctions(original_entity);
      const size_t integrand_order = integrand_.order(original_localFunctions, ansatzBase, testBase) + over_integrate_;
      const auto& volumeQuadrature = QuadratureRules< D, d >::rule(original_entity.type(),
                                                                   boost::numeric_cast< int >(integrand_order));
      for (auto& entity : pv) {
        const auto localFunctions = integrand_.localFunctions(entity);


        const size_t rows = lfs.size();
        const size_t cols = ansatzBase.size();
      }

    } // ... assembleLocal(...)


//    // check matrix and tmp storage
//    const size_t rows = testBase.size();
//    const size_t cols = ansatzBase.size();
//    ret *= 0.0;
//    assert(ret.rows() >= rows);
//    assert(ret.cols() >= cols);
//    assert(tmpLocalMatrices.size() >= numTmpObjectsRequired_);
//    auto& evaluationResult = tmpLocalMatrices[0];
//    // loop over all quadrature points
//    for (const auto& quadPoint : volumeQuadrature) {
//      const auto x = quadPoint.position();
//      // integration factors
//      const auto integrationFactor = entity.geometry().integrationElement(x);
//      const auto quadratureWeight = quadPoint.weight();
//      // evaluate the integrand
//      integrand_.evaluate(localFunctions, ansatzBase, testBase, x, evaluationResult);
//      // compute integral
//      for (size_t ii = 0; ii < rows; ++ii) {
//        auto& retRow = ret[ii];
//        const auto& evaluationResultRow = evaluationResult[ii];
//        for (size_t jj = 0; jj < cols; ++jj)
//          retRow[jj] += evaluationResultRow[jj] * integrationFactor * quadratureWeight;
//      } // compute integral
//    } // loop over all quadrature points
  } // ... apply(...)

private:
  const BinaryEvaluationType integrand_;
  const size_t over_integrate_;
}; // class VirtualRefinedCodim0Integral


} // namespace LocalOperator
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALOPERATOR_INTEGRAL_HH
