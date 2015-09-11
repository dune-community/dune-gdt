// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#warning This header is deprecated, include <dune/gdt/localoperator/interfaces.hh> instead (08.09.2015)!

#ifndef DUNE_GDT_LOCALOPERATOR_INTERFACE_HH
#define DUNE_GDT_LOCALOPERATOR_INTERFACE_HH

#include <vector>

#include <dune/common/deprecated.hh>
#include <dune/common/dynmatrix.hh>

#include <dune/stuff/common/crtp.hh>

#include <dune/gdt/basefunctionset/interface.hh>

namespace Dune {
namespace GDT {
namespace LocalOperator {


template< class Traits >
class
  DUNE_DEPRECATED_MSG("Use LocalOperatorInterface or LocalVolumeTwoFormInterface instead (08.09.2015)!")
      Codim0Interface
  : public Stuff::CRTPInterface< Codim0Interface< Traits >, Traits >
{
public:
  typedef typename Traits::derived_type derived_type;

  size_t numTmpObjectsRequired() const
  {
    CHECK_CRTP(this->as_imp().numTmpObjectsRequired());
    return this->as_imp().numTmpObjectsRequired();
  }

  /**
   *  \brief      Applies the local operator.
   *  \tparam T       Traits of the test BaseFunctionSetInterface implementation
   *  \tparam A       Traits of the ansatz BaseFunctionSetInterface implementation
   *  \tparam D       DomainFieldType
   *  \tparam d       dimDomain
   *  \tparam R       RangeFieldType
   *  \tparam r{T,A}  dimRange of the of the {testBase,ansatzBase}
   *  \tparam rC{T,a} dimRangeCols of the {testBase,ansatzBase}
   */
  template< class T, class A, class D, size_t d, class R, size_t rT, size_t rCT, size_t rA, size_t rCA >
  void apply(const BaseFunctionSetInterface< T, D, d, R, rT, rCT >& testBase,
             const BaseFunctionSetInterface< A, D, d, R, rA, rCA >& ansatzBase,
             Dune::DynamicMatrix< R >& ret,
             std::vector< Dune::DynamicMatrix< R > >& tmpLocalMatrices) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().apply(testBase, ansatzBase, ret, tmpLocalMatrices));
  }
}; // class Codim0Interface


template< class Traits >
class
  DUNE_DEPRECATED_MSG("Use LocalCouplingTwoFormInterface instead (08.09.2015)!")
      Codim1CouplingInterface
  : public Stuff::CRTPInterface< Codim1CouplingInterface< Traits >, Traits >
{
public:
  typedef typename Traits::derived_type derived_type;

  size_t numTmpObjectsRequired() const
  {
    CHECK_CRTP(this->as_imp().numTmpObjectsRequired());
    return this->as_imp().numTmpObjectsRequired();
  }

  /**
   *  \brief Applies the local operator.
   *  \tparam TE      Traits of the entity test BaseFunctionSetInterface implementation
   *  \tparam AE      Traits of the entity ansatz BaseFunctionSetInterface implementation
   *  \tparam TN      Traits of the neighbor test BaseFunctionSetInterface implementation
   *  \tparam AN      Traits of the neighbor ansatz BaseFunctionSetInterface implementation
   *  \tparam IntersectionType
   *  \tparam D       DomainFieldType
   *  \tparam d       dimDomain
   *  \tparam R       RangeFieldType
   *  \tparam r{T,A}  dimRange of the of the {testBase*,ansatzBase*}
   *  \tparam rC{T,a} dimRangeCols of the {testBase*,ansatzBase*}
   */
  template< class TE, class AE, class TN, class AN,
            class IntersectionType,
            class D, size_t d, class R, size_t rT, size_t rCT, size_t rA, size_t rCA >
  void apply(const BaseFunctionSetInterface< TE, D, d, R, rT, rCT >& testBaseEntity,
             const BaseFunctionSetInterface< AE, D, d, R, rA, rCA >& ansatzBaseEntity,
             const BaseFunctionSetInterface< TN, D, d, R, rT, rCT >& testBaseNeighbor,
             const BaseFunctionSetInterface< AN, D, d, R, rA, rCA >& ansatzBaseNeighbor,
             const IntersectionType& intersection,
             Dune::DynamicMatrix< R >& entityEntityRet,
             Dune::DynamicMatrix< R >& neighborNeighborRet,
             Dune::DynamicMatrix< R >& entityNeighborRet,
             Dune::DynamicMatrix< R >& neighborEntityRet,
             std::vector< Dune::DynamicMatrix< R > >& tmpLocalMatrices) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().apply(testBaseEntity, ansatzBaseEntity,
                                             testBaseNeighbor, ansatzBaseNeighbor,
                                             intersection,
                                             entityEntityRet, neighborNeighborRet,
                                             entityNeighborRet, neighborEntityRet,
                                             tmpLocalMatrices));
  }
}; // class Codim1CouplingInterface


template< class Traits >
class
  DUNE_DEPRECATED_MSG("Use LocalBoundaryTwoFormInterface instead (08.09.2015)!")
      Codim1BoundaryInterface
  : public Stuff::CRTPInterface< Codim1BoundaryInterface< Traits >, Traits >
{
public:
  typedef typename Traits::derived_type derived_type;

  size_t numTmpObjectsRequired() const
  {
    CHECK_CRTP(this->as_imp().numTmpObjectsRequired());
    return this->as_imp().numTmpObjectsRequired();
  }

  /**
   *  \brief Applies the local operator.
   *  \tparam T       Traits of the test BaseFunctionSetInterface implementation
   *  \tparam A       Traits of the ansatz BaseFunctionSetInterface implementation
   *  \tparam IntersectionType
   *  \tparam D       DomainFieldType
   *  \tparam d       dimDomain
   *  \tparam R       RangeFieldType
   *  \tparam r{T,A}  dimRange of the of the {testBase,ansatzBase}
   *  \tparam rC{T,a} dimRangeCols of the {testBase,ansatzBase}
   */
  template< class T, class A, class IntersectionType, class D, size_t d, class R, size_t rT, size_t rCT, size_t rA, size_t rCA >
  void apply(const BaseFunctionSetInterface< T, D, d, R, rT, rCT >& testBase,
             const BaseFunctionSetInterface< A, D, d, R, rA, rCA >& ansatzBase,
             const IntersectionType& intersection,
             Dune::DynamicMatrix< R >& ret,
             std::vector< Dune::DynamicMatrix< R > >& tmpLocalMatrices) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().apply(testBase, ansatzBase, intersection, ret, tmpLocalMatrices));
  }
}; // class Codim1BoundaryInterface


} // namespace LocalOperator
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALOPERATOR_INTERFACE_HH
