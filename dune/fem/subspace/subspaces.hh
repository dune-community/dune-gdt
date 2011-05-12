#ifndef DUNE_FEM_FUNCTIONALS_SUBSPACE_SUBSPACES_HH
#define DUNE_FEM_FUNCTIONALS_SUBSPACE_SUBSPACES_HH


namespace Dune {

namespace Functionals {

/**
 * @addtogroup SubSpaces Introduction to Subspaces
 * Let @f$V_h@f$ be a discrete function space and @f$C=\{C_1,\ldots,C_M\}@f$ be a set of
 * constraints on @f$V_h@f$. Then we call
 * @f[V_{h,C}:=\{v_h\in V_h\mid C_i[v_h]=0\quad i\in\{1,\ldots,M\}\}@f]
 * a @link Subspace::Linear @e linear @e subspace @endlink
 * of @f$V_h@f$ with respect to @f$C@f$.
 *
 * Let @f$\hat{g}\in V_h@f$ be a discrete function.
 * Then we call
 * @f[V_{\hat{g}}:=\{v_h+\hat{g}\mid v_h\in V_{h,C}\}@f]
 * an @link Subspace::Affine @e affine @e subspace @endlink
 * with respect to @f$\hat{g}@f$ and @f$V_{h,C}@f$.
 *
 *  See @link Constraints::Dirichlet "dirichlet constraints" @endlink to
 * understand how constraints work.
 */

//! Contains various subspaces.
namespace Subspace {

/**
 * @brief Represents a linear discrete function space with constrained doFs.
 *
 * Let @f$V_h@f$ be a discrete function space and @f$C=\{C_1,\ldots,C_M\}@f$ be a set of
 * constraints on @f$V_h@f$. Then we call
 * @f[V_{h,C}:=\{v_h\in V_h\mid C_i[v_h]=0\quad i\in\{1,\ldots,M\}\}@f]
 * a @e linear @e subspace of @f$V_h@f$ with respect to @f$C@f$.
 *
 * See @link Constraints::Dirichlet "dirichlet constraints" @endlink to
 * understand how constraints work.
 *
 * @note We are expecting an linear @b discrete function space at the moment.
 * It could be generalized in future.
 *
 * @tparam DiscreteFunctionSpaceImp Type of the discrete function space.
 * @tparam ConstraintsImp Type of the constraints.
 *
 * @ingroup SubSpaces
 *
 * @todo should prepare sparsity patterns and such things!
 */
template <class DiscreteFunctionSpaceImp, class ConstraintsImp>
class Linear : public DiscreteFunctionSpaceImp
{
public:
  //! Type of the discrete function space.
  typedef DiscreteFunctionSpaceImp DiscreteFunctionSpaceType;
  //! Type of the constraints.
  typedef ConstraintsImp ConstraintsType;

  //! Traits class.
  typedef typename DiscreteFunctionSpaceType::Traits Traits;
  //! Type of the analytical function space.
  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
  //! Type of the base function set.
  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
  //! Type of the grid part.
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  //! Type of the grid.
  typedef typename DiscreteFunctionSpaceType::GridType GridType;
  //! Type of the index set.
  typedef typename DiscreteFunctionSpaceType::IndexSetType IndexSetType;
  //! Type of the iterator.
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  //! Type of the entities.
  typedef typename DiscreteFunctionSpaceType::EntityType EntityType;
  //! Type of the dof manager.
  typedef typename DiscreteFunctionSpaceType::DofManagerType DofManagerType;
  //! Type of the communication manager.
  typedef typename DiscreteFunctionSpaceType::CommunicationManagerType CommunicationManagerType;
  //! Type of the mapper.
  typedef typename DiscreteFunctionSpaceType::MapperType MapperType;
  //! Type of the block mapper.
  typedef typename DiscreteFunctionSpaceType::BlockMapperType BlockMapperType;

private:
  typedef DiscreteFunctionSpaceType BaseType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;


public:
  //! The local block size.
  enum
  {
    localBlockSize = DiscreteFunctionSpaceType::localBlockSize
  };

  /**
   * @brief Constructor.
   *
   * @param space The discrete function space.
   * @param constraints The constraints.
   */
  Linear(DiscreteFunctionSpaceType& space, ConstraintsType& constraints)
    : DiscreteFunctionSpaceType(space.gridPart())
    , space_(space)
    , constraints_(constraints)
  {
  }

  /**
   * @brief Copy constructor.
   *
   * @param lin A reference to an existing linear subspace @f$V'_{h,C}@f$.
   */
  Linear(const Linear& lin)
    : DiscreteFunctionSpaceType(lin.space().gridPart())
    , space_(lin.space())
    , constraints_(lin.constraints())
  {
  }

  /**
   * @brief Returns the constraints.
   *
   * @return A reference to the constraints @f$C@f$.
   */
  ConstraintsType& constraints() const
  {
    return constraints_;
  }

  /**
   * @brief Returns the discrete function space.
   *
   * @return A reference to the discrete function space @f$V_h@f$.
   */
  DiscreteFunctionSpaceType& space() const
  {
    return space_;
  }

private:
  DiscreteFunctionSpaceType& space_;
  ConstraintsType& constraints_;
}; // end of class Linear


/**
 * @brief Represents an affine discrete function space.
 *
 * Let @f$V_h@f$ be a discrete function space and @f$C=\{C_1,\ldots,C_M\}@f$ be a set of
 * constraints on @f$V_h@f$. Let @f$V_{h,C}@f$, defined by
 * @f[V_{h,C}:=\{v_h\in V_h\mid C_i[v_h]=0\quad i\in\{1,\ldots,M\}\},@f]
 * a linear subspace and let @f$\hat{g}\in V_h@f$ be a discrete function.
 * Then we call
 * @f[V_{\hat{g}}:=\{v_h+\hat{g}\mid v_h\in V_{h,C}\}@f]
 * an @e affine @e subspace with respect to @f$\hat{g}@f$ and @f$V_{h,C}@f$.
 *
 * See @link Subspace::Linear "linear subspace" @endlink to understand how a linear
 * subspace should be created.
 * See @link Constraints::Dirichlet "dirichlet constraints" @endlink to
 * understand how constraints work.
 *
 * @note We are expecting an affine @b discrete function space at the moment.
 * It could be generalized in future.
 *
 * @tparam DiscreteFunctionSpaceImp Type of the discrete function space.
 * @tparam ConstraintsImp Type of the constraints.
 *
 * @ingroup SubSpaces
 *
 * @todo should provide a means to access the offset-function!
 */

template <class LinearSubspaceImp, class OffsetFunctionImp>
class Affine : public LinearSubspaceImp
{
public:
  //! Type of the linear subspace @f$V_{h,C}@f$.
  typedef LinearSubspaceImp LinearSubspaceType;
  //! Type of the offset function @f$\hat{g}@f$.
  typedef OffsetFunctionImp OffsetFunctionType;

  //! Type of the discrete function space @f$V_h@f$.
  typedef typename LinearSubspaceType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  //! Type of the constraints @f$C@f$ extracted from the linear subspace @f$V_{h,C}@f$.
  typedef typename LinearSubspaceType::ConstraintsType ConstraintsType;

  //! Traits class.
  typedef typename LinearSubspaceType::Traits Traits;

  //! Type of the analytical function space.
  typedef typename LinearSubspaceType::FunctionSpaceType FunctionSpaceType;
  //! Type of the basis function set.
  typedef typename LinearSubspaceType::BaseFunctionSetType BaseFunctionSetType;
  //! Type of the grid part.
  typedef typename LinearSubspaceType::GridPartType GridPartType;
  //! Type of the grid.
  typedef typename LinearSubspaceType::GridType GridType;
  //! Type of the index set.
  typedef typename LinearSubspaceType::IndexSetType IndexSetType;
  //! Type of the iterator.
  typedef typename LinearSubspaceType::IteratorType IteratorType;
  //! Type of the iterator iterating over intersections.
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  //! Type of the entity.
  typedef typename LinearSubspaceType::EntityType EntityType;
  //! Type of the dof manager.
  typedef typename LinearSubspaceType::DofManagerType DofManagerType;
  //! Type of the communication manager.
  typedef typename LinearSubspaceType::CommunicationManagerType CommunicationManagerType;
  //! Type of the mapper.
  typedef typename LinearSubspaceType::MapperType MapperType;
  //! Type of the block mapper.
  typedef typename LinearSubspaceType::BlockMapperType BlockMapperType;

private:
  typedef LinearSubspaceImp BaseType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;


public:
  //! The local block size.
  enum
  {
    localBlockSize = LinearSubspaceType::localBlockSize
  };

  /**
   * @brief Constructor.
   *
   * @param linear The linear subspace @f$V_{h,C}@f$.
   * @param offset The discrete function @f$\hat{g}@f$.
   */
  Affine(LinearSubspaceType& linear, OffsetFunctionType& offset)
    : LinearSubspaceType(linear)
    , linear_(linear)
    , offset_(offset)
  {
  }

  /**
   * @brief Copy constructor.
   *
   * @param aff A reference to an existing affine subspace @f$V'_{\hat{g}}@f$.
   */
  Affine(const Affine& aff)
    : LinearSubspaceType(aff.linearSpace())
    , linear_(aff.linearSpace())
    , offset_(aff.offset())
  {
  }

  /**
   * @brief Returns a reference to the linear subspace @f$V_{h,C}@f$.
   *
   * @return The linear subspace.
   */
  LinearSubspaceType& linearSpace() const
  {
    return linear_;
  }

  /**
   * @brief Returns a reference to the discrete function @f$\hat{g}@f$.
   */
  OffsetFunctionType& offset() const
  {
    return offset_;
  }

private:
  LinearSubspaceType& linear_;
  OffsetFunctionType& offset_;
}; // end of class Affine

} // end of namespace Subspace

} // end of namespace Functionals

} // end of namespace Dune


#endif /* end of include guard: DUNE_FEM_FUNCTIONALS_SUBSPACE_SUBSPACES_HH */
