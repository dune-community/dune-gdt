// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_OPERATORS_BASE_HH
#define DUNE_GDT_OPERATORS_BASE_HH

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/grid/walker/apply-on.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/la/container/pattern.hh>
#include <dune/xt/functions/interfaces.hh>

#include <dune/gdt/local/assembler.hh>
#include <dune/gdt/assembler/wrapper.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/operators/interfaces.hh>
#include <dune/gdt/spaces/interface.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forward, required for the traits
template <class M,
          class RS,
          class GV = typename RS::GridViewType,
          class SS = RS,
          class F = typename M::RealType,
          ChoosePattern pt = ChoosePattern::face_and_volume>
class MatrixOperatorBase;


namespace internal {


template <class MatrixImp,
          class RangeSpaceImp,
          class GridViewImp,
          class SourceSpaceImp,
          class FieldImp,
          ChoosePattern pt>
class MatrixOperatorBaseTraits
{
  static_assert(XT::LA::is_matrix<MatrixImp>::value, "MatrixType has to be derived from XT::LA::MatrixInterface!");
  static_assert(is_space<RangeSpaceImp>::value, "RangeSpaceType has to be derived from SpaceInterface!");
  static_assert(is_space<SourceSpaceImp>::value, "SourceSpaceType has to be derived from SpaceInterface!");
  static_assert(std::is_same<typename RangeSpaceImp::GridViewType::template Codim<0>::Entity,
                             typename GridViewImp::template Codim<0>::Entity>::value,
                "RangeSpaceType and GridViewType have to match!");
  static_assert(std::is_same<typename SourceSpaceImp::GridViewType::template Codim<0>::Entity,
                             typename GridViewImp::template Codim<0>::Entity>::value,
                "SourceSpaceType and GridViewType have to match!");

public:
  typedef MatrixOperatorBase<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp, FieldImp, pt> derived_type;
  typedef FieldImp FieldType;
};


} // namespace internal


/**
 * \todo Check parallel case: there is probably/definitely communication missing in apply2!
 */
template <class GridViewImp,
          class RangeImp,
          class SourceImp = RangeImp,
          class FieldImp = typename RangeImp::RangeFieldType>
class LocalizableProductBase : public XT::Grid::Walker<GridViewImp>
{
  typedef LocalizableProductBase<GridViewImp, RangeImp, SourceImp> ThisType;
  typedef XT::Grid::Walker<GridViewImp> BaseType;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::EntityType;
  typedef RangeImp RangeType;
  typedef SourceImp SourceType;
  typedef FieldImp FieldType;

private:
  static_assert(XT::Functions::is_localizable_function<SourceType>::value,
                "SourceType has to be derived from XT::Functions::LocalizableFunctionInterface!");
  static_assert(XT::Functions::is_localizable_function<RangeType>::value,
                "RangeType has to be derived from XT::Functions::LocalizableFunctionInterface!");
  static_assert(std::is_same<typename SourceType::EntityType, EntityType>::value,
                "The EntityType of SourceType and GridViewType have to match!");
  static_assert(std::is_same<typename RangeType::EntityType, EntityType>::value,
                "The EntityType of RangeType and GridViewType have to match!");
  static_assert(std::is_same<typename SourceType::DomainFieldType, typename GridViewType::ctype>::value,
                "The DomainFieldType of SourceType and GridViewType have to match!");
  static_assert(std::is_same<typename RangeType::DomainFieldType, typename GridViewType::ctype>::value,
                "The DomainFieldType of RangeType and GridViewType have to match!");
  static_assert(SourceType::dimDomain == GridViewType::dimension,
                "The dimDomain of SourceType and GridViewType have to match!");
  static_assert(RangeType::dimDomain == GridViewType::dimension,
                "The dimDomain of RangeType and GridViewType have to match!");

public:
  LocalizableProductBase(GridViewType grd_vw, const RangeType& rng, const SourceType& src)
    : BaseType(grd_vw)
    , range_(rng)
    , source_(src)
    , walked_(false)
  {
  }

  LocalizableProductBase(GridViewType grd_vw, const RangeType& rng)
    : BaseType(grd_vw)
    , range_(rng)
    , source_(rng)
    , walked_(false)
  {
  }

  const SourceType& source() const
  {
    return source_;
  }

  const RangeType& range() const
  {
    return range_;
  }

  using BaseType::grid_view;
  using BaseType::append;

  template <class V>
  ThisType&
  append(const LocalVolumeTwoFormInterface<V>& local_volume_twoform,
         const XT::Grid::ApplyOn::WhichEntity<GridViewType>* where = new XT::Grid::ApplyOn::AllEntities<GridViewType>())
  {
    typedef LocalVolumeTwoFormAccumulator<GridViewType,
                                          typename LocalVolumeTwoFormInterface<V>::derived_type,
                                          RangeType,
                                          SourceType,
                                          FieldType>
        AccumulateFunctor;
    local_volume_twoforms_.emplace_back(
        new AccumulateFunctor(grid_view(), local_volume_twoform.as_imp(), range_, source_, *where));
    BaseType::append(*local_volume_twoforms_.back(), where);
    return *this;
  } // ... append(...)

  using BaseType::add;

  template <class... Args>
  DUNE_DEPRECATED_MSG("Use append() instead (since 11.01.2017)!")
  ThisType& add(Args&&... args)
  {
    return append(std::forward<Args>(args)...);
  }

  FieldType compute_locally(const EntityType& entity) const
  {
    FieldType local_result = 0.;
    for (const auto& local_volume_twoform : local_volume_twoforms_)
      local_result += local_volume_twoform->compute_locally(entity);
    return local_result;
  }

  FieldType apply2()
  {
    if (!walked_) {
      this->walk();
      walked_ = true;
    }
    FieldType result = 0.;
    for (const auto& local_volume_twoform : local_volume_twoforms_)
      result += local_volume_twoform->result();
    return result;
  }

protected:
  const RangeType& range_;
  const SourceType& source_;
  std::vector<std::unique_ptr<XT::Grid::internal::Codim0ReturnObject<GridViewType, FieldType>>> local_volume_twoforms_;
  bool walked_;
}; // class LocalizableProductBase


/**
 * \todo add static checks of dimensions
 * \note Does a const_cast in apply() and apply2(), not sure yet if this is fine.
 */
template <class MatrixImp,
          class RangeSpaceImp,
          class GridViewImp,
          class SourceSpaceImp,
          class FieldImp,
          ChoosePattern pt>
class MatrixOperatorBase : public OperatorInterface<internal::MatrixOperatorBaseTraits<MatrixImp,
                                                                                       RangeSpaceImp,
                                                                                       GridViewImp,
                                                                                       SourceSpaceImp,
                                                                                       FieldImp,
                                                                                       pt>>,
                           public SystemAssembler<RangeSpaceImp, GridViewImp, SourceSpaceImp>
{
  typedef OperatorInterface<internal::MatrixOperatorBaseTraits<MatrixImp,
                                                               RangeSpaceImp,
                                                               GridViewImp,
                                                               SourceSpaceImp,
                                                               FieldImp,
                                                               pt>>
      BaseOperatorType;
  typedef SystemAssembler<RangeSpaceImp, GridViewImp, SourceSpaceImp> BaseAssemblerType;
  typedef MatrixOperatorBase<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp, FieldImp, pt> ThisType;

public:
  typedef internal::MatrixOperatorBaseTraits<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp, FieldImp, pt>
      Traits;
  typedef typename BaseAssemblerType::TestSpaceType RangeSpaceType;
  typedef typename BaseAssemblerType::AnsatzSpaceType SourceSpaceType;
  typedef XT::LA::SparsityPatternDefault PatternType;
  typedef MatrixImp MatrixType;
  using typename BaseOperatorType::FieldType;
  using typename BaseOperatorType::derived_type;
  using typename BaseAssemblerType::GridViewType;

private:
  typedef XT::LA::Solver<MatrixType, typename SourceSpaceType::CommunicatorType> LinearSolverType;

  template <ChoosePattern pp = ChoosePattern::face_and_volume, bool anything = true>
  struct Compute
  {
    static PatternType
    pattern(const RangeSpaceType& rng_spc, const SourceSpaceType& src_spc, const GridViewType& grd_vw)
    {
      return rng_spc.compute_face_and_volume_pattern(grd_vw, src_spc);
    }
  };

  template <bool anything>
  struct Compute<ChoosePattern::volume, anything>
  {
    static PatternType
    pattern(const RangeSpaceType& rng_spc, const SourceSpaceType& src_spc, const GridViewType& grd_vw)
    {
      return rng_spc.compute_volume_pattern(grd_vw, src_spc);
    }
  };

  template <bool anything>
  struct Compute<ChoosePattern::face, anything>
  {
    static PatternType
    pattern(const RangeSpaceType& rng_spc, const SourceSpaceType& src_spc, const GridViewType& grd_vw)
    {
      return rng_spc.compute_face_pattern(grd_vw, src_spc);
    }
  };

public:
  static PatternType pattern(const RangeSpaceType& rng_spc, const SourceSpaceType& src_spc, const GridViewType& grd_vw)
  {
    return Compute<pt>::pattern(rng_spc, src_spc, grd_vw);
  }

  static PatternType pattern(const RangeSpaceType& rng_spc)
  {
    return pattern(rng_spc, rng_spc);
  }

  static PatternType pattern(const RangeSpaceType& rng_spc, const SourceSpaceType& src_spc)
  {
    return pattern(rng_spc, src_spc, rng_spc.grid_view());
  }

  static PatternType pattern(const RangeSpaceType& rng_spc, const GridViewType& grd_vw)
  {
    return pattern(rng_spc, rng_spc, grd_vw);
  }

  template <class... Args>
  explicit MatrixOperatorBase(MatrixType& mtrx, Args&&... args)
    : BaseAssemblerType(std::forward<Args>(args)...)
    , matrix_(mtrx)
  {
    if (matrix_.access().rows() != this->range_space().mapper().size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "matrix.rows(): " << matrix_.access().rows() << "\n"
                                   << "range_space().mapper().size(): "
                                   << this->range_space().mapper().size());
    if (matrix_.access().cols() != this->source_space().mapper().size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "matrix.cols(): " << matrix_.access().cols() << "\n"
                                   << "source_space().mapper().size(): "
                                   << this->source_space().mapper().size());
  } // MatrixOperatorBase(...)

  /// \todo Guard against copy and move ctor (Args = ThisType)!
  template <class... Args>
  explicit MatrixOperatorBase(Args&&... args)
    : BaseAssemblerType(std::forward<Args>(args)...)
    , matrix_(new MatrixType(this->range_space().mapper().size(),
                             this->source_space().mapper().size(),
                             pattern(this->range_space(), this->source_space(), this->grid_view())))
  {
  }

  /// \sa SystemAssembler
  MatrixOperatorBase(const ThisType& other) = delete;
  MatrixOperatorBase(ThisType&& source) = delete;
  MatrixOperatorBase(ThisType& other) = delete; // <- b.c. of the too perfect forwarding ctor

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  const MatrixType& matrix() const
  {
    return matrix_.access();
  }

  MatrixType& matrix()
  {
    return matrix_.access();
  }

  const SourceSpaceType& source_space() const
  {
    return this->ansatz_space();
  }

  const RangeSpaceType& range_space() const
  {
    return this->test_space();
  }

  using BaseAssemblerType::append;

  template <class V>
  ThisType&
  append(const LocalVolumeTwoFormInterface<V>& local_volume_twoform,
         const XT::Grid::ApplyOn::WhichEntity<GridViewType>* where = new XT::Grid::ApplyOn::AllEntities<GridViewType>())
  {
    typedef internal::LocalVolumeTwoFormWrapper<ThisType,
                                                typename LocalVolumeTwoFormInterface<V>::derived_type,
                                                MatrixType>
        WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(
        this->test_space_, this->ansatz_space_, where, local_volume_twoform.as_imp(), matrix_.access()));
    return *this;
  } // ... append(...)

  template <class C>
  ThisType& append(const LocalCouplingTwoFormInterface<C>& local_coupling_twoform,
                   const XT::Grid::ApplyOn::WhichIntersection<GridViewType>* where =
                       new XT::Grid::ApplyOn::InnerIntersectionsPrimally<GridViewType>())
  {
    typedef internal::LocalCouplingTwoFormWrapper<ThisType,
                                                  typename LocalCouplingTwoFormInterface<C>::derived_type,
                                                  MatrixType>
        WrapperType;
    this->codim1_functors_.emplace_back(new WrapperType(
        this->test_space_, this->ansatz_space_, where, local_coupling_twoform.as_imp(), matrix_.access()));
    return *this;
  } // ... append(...)

  template <class B>
  ThisType& append(const LocalBoundaryTwoFormInterface<B>& local_boundary_twoform,
                   const XT::Grid::ApplyOn::WhichIntersection<GridViewType>* where =
                       new XT::Grid::ApplyOn::InnerIntersectionsPrimally<GridViewType>())
  {
    typedef internal::LocalBoundaryTwoFormWrapper<ThisType,
                                                  typename LocalBoundaryTwoFormInterface<B>::derived_type,
                                                  MatrixType>
        WrapperType;
    this->codim1_functors_.emplace_back(new WrapperType(
        this->test_space_, this->ansatz_space_, where, local_boundary_twoform.as_imp(), matrix_.access()));
    return *this;
  } // ... append(...)

  using BaseAssemblerType::add;

  template <class... Args>
  DUNE_DEPRECATED_MSG("Use append() instead (since 11.01.2017)!")
  ThisType& add(Args&&... args)
  {
    return append(std::forward<Args>(args)...);
  }

  template <class S, class R>
  void apply(const XT::LA::VectorInterface<S>& source, XT::LA::VectorInterface<R>& range) const
  {
    const_cast<ThisType&>(*this).assemble();
    matrix().mv(source.as_imp(), range.as_imp());
  }

  template <class S, class R>
  void apply(const ConstDiscreteFunction<SourceSpaceType, S>& source, DiscreteFunction<RangeSpaceType, R>& range) const
  {
    apply(source.vector(), range.vector());
  }

  template <class R, class S>
  FieldType apply2(const XT::LA::VectorInterface<R>& range, const XT::LA::VectorInterface<S>& source) const
  {
    const_cast<ThisType&>(*this).assemble();
    auto tmp = range.copy();
    matrix().mv(source.as_imp(source), tmp);
    return range.dot(tmp);
  }

  template <class R, class S>
  FieldType apply2(const ConstDiscreteFunction<RangeSpaceType, R>& range,
                   const ConstDiscreteFunction<SourceSpaceType, S>& source) const
  {
    return apply2(range.vector(), source.vector());
  }

  //! \todo Implement a base for matrix operators that only gets a matrix and handles the apply, apply2, apply_inverse,
  //! etc.
  //  template< class SourceType >
  //  JacobianType jacobian(const SourceType& /*source*/) const
  //  {
  //    return JacobianType(matrix());
  //  }

  using BaseOperatorType::apply_inverse;

  template <class R, class S>
  void apply_inverse(const XT::LA::VectorInterface<R>& range,
                     XT::LA::VectorInterface<S>& source,
                     const XT::Common::Configuration& opts) const
  {
    this->assemble();
    LinearSolverType(matrix(), source_space().communicator()).apply(range.as_imp(), source.as_imp(), opts);
  }

  template <class R, class S>
  void apply_inverse(const ConstDiscreteFunction<RangeSpaceType, R>& range,
                     ConstDiscreteFunction<SourceSpaceType, S>& source,
                     const XT::Common::Configuration& opts) const
  {
    apply_inverse(range.vector(), source.vector(), opts);
  }

  std::vector<std::string> invert_options() const
  {
    return LinearSolverType::types();
  }

  XT::Common::Configuration invert_options(const std::string& type) const
  {
    return LinearSolverType::options(type);
  }

protected:
  using BaseAssemblerType::codim0_functors_;
  using BaseAssemblerType::codim1_functors_;

private:
  Dune::XT::Common::StorageProvider<MatrixType> matrix_;
}; // class MatrixOperatorBase


template <class GridViewImp, class SourceImp, class RangeImp>
class LocalizableOperatorBase : public XT::Grid::Walker<GridViewImp>
{
  typedef LocalizableOperatorBase<GridViewImp, SourceImp, RangeImp> ThisType;
  typedef XT::Grid::Walker<GridViewImp> BaseType;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::EntityType;
  typedef SourceImp SourceType;
  typedef RangeImp RangeType;

private:
  static_assert(XT::Functions::is_localizable_function<SourceType>::value,
                "SourceType has to be derived from XT::Functions::LocalizableFunctionInterface!");
  static_assert(is_discrete_function<RangeType>::value, "RangeType has to be a DiscreteFunctionDefault!");
  static_assert(std::is_same<typename SourceType::EntityType, EntityType>::value, "Have to match!");
  static_assert(std::is_same<typename RangeType::EntityType, EntityType>::value, "Have to match!");
  static_assert(std::is_same<typename SourceType::DomainFieldType, typename GridViewType::ctype>::value,
                "Have to match!");
  static_assert(std::is_same<typename RangeType::DomainFieldType, typename GridViewType::ctype>::value,
                "Have to match!");
  static_assert(SourceType::dimDomain == GridViewType::dimension, "Have to match!");
  static_assert(RangeType::dimDomain == GridViewType::dimension, "Have to match!");

public:
  LocalizableOperatorBase(GridViewType grd_vw, const SourceType& src, RangeType& rng)
    : BaseType(grd_vw)
    , source_(src)
    , range_(rng)
    , walked_(false)
  {
  }

  const SourceType& source() const
  {
    return source_;
  }

  const RangeType& range() const
  {
    return range_;
  }

  RangeType& range()
  {
    return range_;
  }

  using BaseType::grid_view;
  using BaseType::append;

  template <class L>
  ThisType&
  append(const LocalOperatorInterface<L>& local_operator,
         const XT::Grid::ApplyOn::WhichEntity<GridViewType>* where = new XT::Grid::ApplyOn::AllEntities<GridViewType>())
  {
    typedef LocalOperatorApplicator<GridViewType,
                                    typename LocalOperatorInterface<L>::derived_type,
                                    SourceType,
                                    RangeType>
        Applicator;
    local_operators_codim_0.emplace_back(new Applicator(grid_view(), local_operator.as_imp(), source_, range_, *where));
    BaseType::append(*local_operators_codim_0.back(), where);
  } // ... append(...)

  template <class T>
  ThisType& append(const LocalCouplingOperatorInterface<T>& local_operator,
                   const XT::Grid::ApplyOn::WhichIntersection<GridViewType>* where =
                       new XT::Grid::ApplyOn::InnerIntersections<GridViewType>())
  {
    typedef LocalCouplingOperatorApplicator<GridViewType,
                                            typename LocalCouplingOperatorInterface<T>::derived_type,
                                            SourceType,
                                            RangeType>
        Applicator;
    local_operators_codim_1.emplace_back(new Applicator(grid_view(), local_operator.as_imp(), source_, range_, *where));
    BaseType::append(*local_operators_codim_1.back(), where);
  } // ... append(...)

  template <class T>
  ThisType& append(const LocalBoundaryOperatorInterface<T>& local_operator,
                   const XT::Grid::ApplyOn::WhichIntersection<GridViewType>* where =
                       new XT::Grid::ApplyOn::BoundaryIntersections<GridViewType>())
  {
    typedef LocalBoundaryOperatorApplicator<GridViewType,
                                            typename LocalBoundaryOperatorInterface<T>::derived_type,
                                            SourceType,
                                            RangeType>
        Applicator;
    local_operators_codim_1.emplace_back(new Applicator(grid_view(), local_operator.as_imp(), source_, range_, *where));
    BaseType::append(*local_operators_codim_1.back(), where);
  } // ... append(...)

  using BaseType::add;

  template <class... Args>
  DUNE_DEPRECATED_MSG("Use append() instead (since 11.01.2017)!")
  ThisType& add(Args&&... args)
  {
    return append(std::forward<Args>(args)...);
  }

  void apply()
  {
    if (!walked_) {
      this->walk();
      walked_ = true;
    }
  } // ... apply(...)

protected:
  const SourceType& source_;
  RangeType& range_;
  std::vector<std::unique_ptr<XT::Grid::internal::Codim0Object<GridViewType>>> local_operators_codim_0;
  std::vector<std::unique_ptr<XT::Grid::internal::Codim1Object<GridViewType>>> local_operators_codim_1;
  bool walked_;
}; // class LocalizableOperatorBase


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_BASE_HH
