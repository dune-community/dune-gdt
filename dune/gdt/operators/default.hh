// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_DEFAULT_HH
#define DUNE_GDT_OPERATORS_DEFAULT_HH

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/grid/walker/apply-on.hh>
#include <dune/stuff/grid/walker.hh>
#include <dune/stuff/la/container/pattern.hh>

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
template <class M, class RS, class GV = typename RS::GridViewType, class SS = RS, class F = typename M::RealType,
          ChoosePattern pt = ChoosePattern::face_and_volume>
class MatrixOperatorDefault;


namespace internal {


template <class MatrixImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp, class FieldImp,
          ChoosePattern pt>
class MatrixOperatorDefaultTraits
{
  static_assert(Stuff::LA::is_matrix<MatrixImp>::value,
                "MatrixType has to be derived from Stuff::LA::MatrixInterface!");
  static_assert(is_space<RangeSpaceImp>::value, "RangeSpaceType has to be derived from SpaceInterface!");
  static_assert(is_space<SourceSpaceImp>::value, "SourceSpaceType has to be derived from SpaceInterface!");
  static_assert(std::is_same<typename RangeSpaceImp::GridViewType::template Codim<0>::Entity,
                             typename GridViewImp::template Codim<0>::Entity>::value,
                "RangeSpaceType and GridViewType have to match!");
  static_assert(std::is_same<typename SourceSpaceImp::GridViewType::template Codim<0>::Entity,
                             typename GridViewImp::template Codim<0>::Entity>::value,
                "SourceSpaceType and GridViewType have to match!");

public:
  typedef MatrixOperatorDefault<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp, FieldImp, pt> derived_type;
  typedef FieldImp FieldType;
};


} // namespace internal


/**
 * \todo Check parallel case: there is probably/definitely communication missing in apply2!
 */
template <class GridViewImp, class RangeImp, class SourceImp = RangeImp,
          class FieldImp                                     = typename RangeImp::RangeFieldType>
class LocalizableProductDefault : public Stuff::Grid::Walker<GridViewImp>
{
  typedef Stuff::Grid::Walker<GridViewImp> BaseType;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::EntityType;
  typedef RangeImp RangeType;
  typedef SourceImp SourceType;
  typedef FieldImp FieldType;

private:
  static_assert(Stuff::is_localizable_function<SourceType>::value,
                "SourceType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(Stuff::is_localizable_function<RangeType>::value,
                "RangeType has to be derived from Stuff::LocalizableFunctionInterface!");
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
  LocalizableProductDefault(GridViewType grd_vw, const RangeType& rng, const SourceType& src)
    : BaseType(grd_vw)
    , range_(rng)
    , source_(src)
    , walked_(false)
  {
  }

  LocalizableProductDefault(GridViewType grd_vw, const RangeType& rng)
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

  template <class V>
  void add(const LocalVolumeTwoFormInterface<V>& local_volume_twoform,
           const DSG::ApplyOn::WhichEntity<GridViewType>* where = new DSG::ApplyOn::AllEntities<GridViewType>())
  {
    typedef LocalVolumeTwoFormAccumulator<GridViewType,
                                          typename LocalVolumeTwoFormInterface<V>::derived_type,
                                          RangeType,
                                          SourceType,
                                          FieldType> AccumulateFunctor;
    local_volume_twoforms_.emplace_back(
        new AccumulateFunctor(grid_view(), local_volume_twoform.as_imp(), range_, source_, *where));
    BaseType::add(*local_volume_twoforms_.back(), where);
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
  std::vector<std::unique_ptr<DSG::internal::Codim0ReturnObject<GridViewType, FieldType>>> local_volume_twoforms_;
  bool walked_;
}; // class LocalizableProductDefault


/**
 * \todo add static checks of dimensions
 * \note Does a const_cast in apply() and apply2(), not sure yet if this is fine.
 */
template <class MatrixImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp, class FieldImp,
          ChoosePattern pt>
class MatrixOperatorDefault
    : public OperatorInterface<internal::MatrixOperatorDefaultTraits<MatrixImp, RangeSpaceImp, GridViewImp,
                                                                     SourceSpaceImp, FieldImp, pt>>,
      public SystemAssembler<RangeSpaceImp, GridViewImp, SourceSpaceImp>
{
  typedef OperatorInterface<internal::MatrixOperatorDefaultTraits<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp,
                                                                  FieldImp, pt>> BaseOperatorType;
  typedef SystemAssembler<RangeSpaceImp, GridViewImp, SourceSpaceImp> BaseAssemblerType;
  typedef MatrixOperatorDefault<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp, FieldImp, pt> ThisType;

public:
  typedef internal::MatrixOperatorDefaultTraits<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp, FieldImp, pt>
      Traits;
  typedef typename BaseAssemblerType::TestSpaceType RangeSpaceType;
  typedef typename BaseAssemblerType::AnsatzSpaceType SourceSpaceType;
  typedef Stuff::LA::SparsityPatternDefault PatternType;
  typedef MatrixImp MatrixType;
  using typename BaseOperatorType::FieldType;
  using typename BaseOperatorType::derived_type;
  using typename BaseAssemblerType::GridViewType;

private:
  typedef Stuff::LA::Solver<MatrixType, typename SourceSpaceType::CommunicatorType> LinearSolverType;

  template <ChoosePattern pp = ChoosePattern::face_and_volume, bool anything = true>
  struct Compute
  {
    static PatternType pattern(const RangeSpaceType& rng_spc, const SourceSpaceType& src_spc,
                               const GridViewType& grd_vw)
    {
      return rng_spc.compute_face_and_volume_pattern(grd_vw, src_spc);
    }
  };

  template <bool anything>
  struct Compute<ChoosePattern::volume, anything>
  {
    static PatternType pattern(const RangeSpaceType& rng_spc, const SourceSpaceType& src_spc,
                               const GridViewType& grd_vw)
    {
      return rng_spc.compute_volume_pattern(grd_vw, src_spc);
    }
  };

  template <bool anything>
  struct Compute<ChoosePattern::face, anything>
  {
    static PatternType pattern(const RangeSpaceType& rng_spc, const SourceSpaceType& src_spc,
                               const GridViewType& grd_vw)
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
  explicit MatrixOperatorDefault(MatrixType& mtrx, Args&&... args)
    : BaseAssemblerType(std::forward<Args>(args)...)
    , matrix_(mtrx)
  {
    if (matrix_.access().rows() != this->range_space().mapper().size())
      DUNE_THROW(Stuff::Exceptions::shapes_do_not_match,
                 "matrix.rows(): " << matrix_.access().rows() << "\n"
                                   << "range_space().mapper().size(): "
                                   << this->range_space().mapper().size());
    if (matrix_.access().cols() != this->source_space().mapper().size())
      DUNE_THROW(Stuff::Exceptions::shapes_do_not_match,
                 "matrix.cols(): " << matrix_.access().cols() << "\n"
                                   << "source_space().mapper().size(): "
                                   << this->source_space().mapper().size());
  } // MatrixOperatorDefault(...)

  template <class... Args>
  explicit MatrixOperatorDefault(Args&&... args)
    : BaseAssemblerType(std::forward<Args>(args)...)
    , matrix_(new MatrixType(this->range_space().mapper().size(), this->source_space().mapper().size(),
                             pattern(this->range_space(), this->source_space(), this->grid_view())))
  {
  }

  MatrixOperatorDefault(ThisType&& source) = default;

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

  using BaseAssemblerType::add;

  template <class V>
  void add(const LocalVolumeTwoFormInterface<V>& local_volume_twoform,
           const DSG::ApplyOn::WhichEntity<GridViewType>* where = new DSG::ApplyOn::AllEntities<GridViewType>())
  {
    typedef internal::LocalVolumeTwoFormWrapper<ThisType,
                                                typename LocalVolumeTwoFormInterface<V>::derived_type,
                                                MatrixType> WrapperType;
    this->codim0_functors_.emplace_back(new WrapperType(
        this->test_space_, this->ansatz_space_, where, local_volume_twoform.as_imp(), matrix_.access()));
  }

  template <class C>
  void add(const LocalCouplingTwoFormInterface<C>& local_coupling_twoform,
           const DSG::ApplyOn::WhichIntersection<GridViewType>* where =
               new DSG::ApplyOn::InnerIntersectionsPrimally<GridViewType>())
  {
    typedef internal::LocalCouplingTwoFormWrapper<ThisType,
                                                  typename LocalCouplingTwoFormInterface<C>::derived_type,
                                                  MatrixType> WrapperType;
    this->codim1_functors_.emplace_back(new WrapperType(
        this->test_space_, this->ansatz_space_, where, local_coupling_twoform.as_imp(), matrix_.access()));
  }

  template <class B>
  void add(const LocalBoundaryTwoFormInterface<B>& local_boundary_twoform,
           const DSG::ApplyOn::WhichIntersection<GridViewType>* where =
               new DSG::ApplyOn::InnerIntersectionsPrimally<GridViewType>())
  {
    typedef internal::LocalBoundaryTwoFormWrapper<ThisType,
                                                  typename LocalBoundaryTwoFormInterface<B>::derived_type,
                                                  MatrixType> WrapperType;
    this->codim1_functors_.emplace_back(new WrapperType(
        this->test_space_, this->ansatz_space_, where, local_boundary_twoform.as_imp(), matrix_.access()));
  }

  template <class S, class R>
  void apply(const Stuff::LA::VectorInterface<S>& source, Stuff::LA::VectorInterface<R>& range) const
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
  FieldType apply2(const Stuff::LA::VectorInterface<R>& range, const Stuff::LA::VectorInterface<S>& source) const
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
  void apply_inverse(const Stuff::LA::VectorInterface<R>& range, Stuff::LA::VectorInterface<S>& source,
                     const Stuff::Common::Configuration& opts) const
  {
    this->assemble();
    LinearSolverType(matrix(), source_space().communicator()).apply(range.as_imp(), source.as_imp(), opts);
  }

  template <class R, class S>
  void apply_inverse(const ConstDiscreteFunction<SourceSpaceType, R>& range,
                     ConstDiscreteFunction<RangeSpaceType, S>& source, const Stuff::Common::Configuration& opts) const
  {
    apply_inverse(range.vector(), source.vector(), opts);
  }

  std::vector<std::string> invert_options() const
  {
    return LinearSolverType::types();
  }

  Stuff::Common::Configuration invert_options(const std::string& type) const
  {
    return LinearSolverType::options(type);
  }

protected:
  using BaseAssemblerType::codim0_functors_;
  using BaseAssemblerType::codim1_functors_;

private:
  DSC::StorageProvider<MatrixType> matrix_;
}; // class MatrixOperatorDefault


template <class GridViewImp, class SourceImp, class RangeImp>
class LocalizableOperatorDefault : public Stuff::Grid::Walker<GridViewImp>
{
  typedef Stuff::Grid::Walker<GridViewImp> BaseType;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::EntityType;
  typedef SourceImp SourceType;
  typedef RangeImp RangeType;

private:
  static_assert(Stuff::is_localizable_function<SourceType>::value,
                "SourceType has to be derived from Stuff::LocalizableFunctionInterface!");
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
  LocalizableOperatorDefault(GridViewType grd_vw, const SourceType& src, RangeType& rng)
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

  using BaseType::add;
  using BaseType::grid_view;

  template <class L>
  void add(const LocalOperatorInterface<L>& local_operator,
           const DSG::ApplyOn::WhichEntity<GridViewType>* where = new DSG::ApplyOn::AllEntities<GridViewType>())
  {
    typedef LocalOperatorApplicator<GridViewType,
                                    typename LocalOperatorInterface<L>::derived_type,
                                    SourceType,
                                    RangeType> Applicator;
    local_operators_codim_0.emplace_back(new Applicator(grid_view(), local_operator.as_imp(), source_, range_, *where));
    BaseType::add(*local_operators_codim_0.back(), where);
  } // ... add(...)

  template <class T>
  void
  add(const LocalCouplingOperatorInterface<T>& local_operator,
      const DSG::ApplyOn::WhichIntersection<GridViewType>* where = new DSG::ApplyOn::InnerIntersections<GridViewType>())
  {
    typedef LocalCouplingOperatorApplicator<GridViewType,
                                            typename LocalCouplingOperatorInterface<T>::derived_type,
                                            SourceType,
                                            RangeType> Applicator;
    local_operators_codim_1.emplace_back(new Applicator(grid_view(), local_operator.as_imp(), source_, range_, *where));
    BaseType::add(*local_operators_codim_1.back(), where);
  } // ... add(...)

  template <class T>
  void add(const LocalBoundaryOperatorInterface<T>& local_operator,
           const DSG::ApplyOn::WhichIntersection<GridViewType>* where =
               new DSG::ApplyOn::BoundaryIntersections<GridViewType>())
  {
    typedef LocalBoundaryOperatorApplicator<GridViewType,
                                            typename LocalBoundaryOperatorInterface<T>::derived_type,
                                            SourceType,
                                            RangeType> Applicator;
    local_operators_codim_1.emplace_back(new Applicator(grid_view(), local_operator.as_imp(), source_, range_, *where));
    BaseType::add(*local_operators_codim_1.back(), where);
  } // ... add(...)


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
  std::vector<std::unique_ptr<DSG::internal::Codim0Object<GridViewType>>> local_operators_codim_0;
  std::vector<std::unique_ptr<DSG::internal::Codim1Object<GridViewType>>> local_operators_codim_1;
  bool walked_;
}; // class LocalizableOperatorDefault


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_DEFAULT_HH
