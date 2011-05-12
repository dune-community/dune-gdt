#ifndef DUNE_FEM_FUNCTIONALS_SOLVER_FEMASSEMBLER_HH
#define DUNE_FEM_FUNCTIONALS_SOLVER_FEMASSEMBLER_HH

#include <dune/fem/common/localmatrix.hh>
#include <dune/fem/common/localvector.hh>

namespace Dune
{

namespace Functionals
{


/**
 * @mainpage Finite element solution of elliptic boundary value problems
 *
 *  @section Start How to start
 *
 *  A good start could be the documentation
 *  <a href="../../../doc/latex/dune-fem-functionals.pdf">dune-fem-functionals.pdf</a>.
 * 
 *  See <a href="modules.html">modules</a>
 *  to keep track of the main ingredients of @c dunefemfunctionals.
 *
 *  We consequently make use of namespaces to structure our concept,
 *  so you should have a look at our 
 *  <a href="namespaces.html">namespaces</a>.
 *
 *  To go in detail you should have a look at the 
 *  <a href="annotated.html">full class list</a>.
 *
 *  Otherwise you can also have a look at our doxygen introduction below...
 *
 *  @section introduction
 *
 *  Assuming standard notation, the following elliptic PDE is one 
 *  of the simplest sample problem we would like to solve
 *  with @c dunefemfunctionals.
 *
 *  @subsection EllipticBoundaryValueProblem Elliptic boundary value problem
 *  @anchor exampleintroductionelliptic_pde 
 *
 *  Let @f$\Omega \subset \mathbb{R}@f$ be a bounded connected lipschitz-domain 
 *  and let @f$a,f : \Omega \rightarrow \mathbb{R}@f$ and
 *  @f$g: \partial\Omega \rightarrow \mathbb{R}@f$ be given functions.
 *
 *  Find @f$u: \Omega \rightarrow \mathbb{R}@f$, such that
 *  @anchor equationintroductionelliptic_pde
 *  @f{align*}{
 *  - \nabla \cdot
 *    (a \nabla u )
 *    &= f &&\text{in } \Omega, \tag{1}\\
 *  u &= g &&\text{on } \partial\Omega.
 *  @f}
 *
 *  The following algorithm gives rise to the corresponding 
 *  constructs that we will need in @c dunefemfunctionals.
 *
 *  @subsection WeakFormulation Weak Formulation
 *  @anchor definitionintroductionweak_formulation
 *  Let @f$H^1@f$ and @f$H^1_0@f$ be given as usual and 
 *  let the affine subspace @f$H^1_g@f$ be defined as
 *  @f{align*}{
 *    H^1_g := \{ v \in H^1 \mid  v = v_0 + \hat{g} \text{ for a } v_0 \in H^1_0  \},
 *  @f}
 *  where @f${\hat{g} \in H^{1}}@f$ is a @f${H^{1}}@f$ representation of @f$g@f$.
 *  The weak formulation of the
 *  @ref equationintroductionelliptic_pde "Elliptic boundary value problem" 
 *  then reads as follows. 
 *
 *  Find @f${u \in H^1_g}@f$, such that
 *  @anchor equationintroductionweak_formulation
 *  @f{align*}{
 *    \int\limits_{\Omega}
 *      a \nabla u \nabla v
 *    \mathrm{d}x
 *    =
 *    \int\limits_{\Omega}
 *      f v
 *    \mathrm{d}x
 *    &&\text{for all } v \in H^1_0.
 *  @f}
 *  
 *  Using the fact, that the solution lies in @f${H_g^1}@f$,
 *  we can decomposed the solution as @f${u = u_0 + \hat{g}}@f$, where
 *  @f${u_0 \in H_0^1}@f$ is the solution to problem 
 *  @ref equationintroductionelliptic_pde "elliptic boundary value problem" 
 *  with homogeneous boundary conditions and a modified right hand side. 
 *  We can thus rewrite the
 *  @ref equationintroductionweak_formulation "weak formulation 
 *  as follows.
 * 
 *  @anchor remarkintroductionweak_formulation 
 *  With the notation from definition @ref definitionintroductionweak_formulation 
 *  "weak formulation", find @f${u_0 \in H_0^1}@f$, such that
 *  @f{align*}{
 *    \int\limits_{\Omega} a\nabla u_0 \nabla v \mathrm{d}x
 *    = \int\limits_{\Omega} f v \mathrm{d}x -
 *    \int\limits_{\Omega} a \nabla \hat{g} \nabla v  \mathrm{d}x
 *    &&\text{for all } v \in H^1_0.
 *  @f}
 *  The weak solution from definition @ref definitionintroductionweak_formulation 
 *  "weak formulation" can then be obtained by
 *  @f{align*}{
 *    u = u_0 + \hat{g}.
 *  @f}
 * 
 *  @subsection Finite_Element_Discretization Finite element discretization
 *  @anchor definitionintroductionfinite_element_discretization
 *  With the notation from definition @ref definitionintroductionweak_formulation
 *  "weak formulation", let @f$ \tau_h@f$ be a conform admissable
 *  triangulation of the domain @f$\Omega@f$ with codim 0 elements 
 *  @f$T \in \tau_h@f$. The usual finite element lagrange spaces are then 
 *  given by
 *
 *  @anchor equationintroductionlagrange_space
 *  @f{align*}{
 *    S_h^k &:=\{ v_h \in C^0(\Omega) \mid 
 *      v_h|_T \in \mathbb{P}^k(T) \quad\forall T \in \tau_h
 *      \},\\
 *    {S_h^k}_{0} &:=
 *      \{ v_h \in S_h^k \mid
 *        v_h = 0 \text{ on } \partial \Omega
 *      \},\\
 *    {S_h^k}_g &:=
 *      \{v_h \in S_h^k \mid
 *        {v_h}_0 + g_h \text{ for a } {v_h}_0 \in {S_h^k}_0
 *      \},
 *  @f}
 *  where @f$g_{h} \in S_{h}^{k}@f$ is the projection of @f$\hat{g}@f$ onto @f${S_{h}^{k}}@f$.
 *
 *  With these discrete function spaces at hand, we can define the finite element solution of the
 *  @ref equationintroductionelliptic_pde "elliptic boundary value problem".
 *
 *  @subsection finite_element_solution Finite element solution
 *  @anchor definitionintroductionfinite_element_solution
 *  With the notation from remark @ref remarkintroductionweak_formulation 
 *  "weak formulation" and definition 
 *  @ref definitionintroductionfinite_element_discretization 
 *  "finite element discretization",
 *  find @f${{u_{h}}_{0} \in {S_{h}^{1}}_0}@f$, such that
 *  @anchor equationintroductionfinite_element_solution_0
 *  @f{align*}{
 *    \int\limits_{\Omega}
 *      a \nabla {u_{h}}_{0} \nabla v_h
 *    \mathrm{d}x =
 *    \int\limits_{\Omega}
 *      f v_h
 *    \mathrm{d}x -
 *    \int\limits_{\Omega}
 *      a \nabla g_h \nabla v_h
 *    \mathrm{d}x
 *    &&\text{for all } v_{h} \in {S_{h}^{1}}_0.
 *  @f}
 *  The finite element solution @f${u_h \in {S_h^1}}@f$ of the 
 *  @ref equationintroductionelliptic_pde "elliptic boundary value problem"
 *  is then given by
 *  @anchor equationintroductionfinite_element_solution
 *  @f{align*}
 *      u_{h} := {u_{h}}_{0} + g_{h}.
 *  @f}
 *
 *  The @ref equationintroductionweak_formulation
 *  "weak formulation" and the @ref equationintroductionfinite_element_solution_0
 *  "finite element solution" gives rise to the introduction of
 *  functionals and operators. A rigorous mathematical definition 
 *  of these can be found in the next section.
 *  The following is only intended to give the basic idea.
 *
 *
 *  @subsection operator_and_functional Operator and functional
 *  @anchor definitionintroductionoperators_functionals
 *  The function @f$f@f$ from example @ref exampleintroductionelliptic_pde 
 *  "ellptic boundary value problem" induces a functional
 *  @f{align*}{
 *    F: S_h^k &\rightarrow \mathbb{R}
 *      \notag\\
 *    v &\mapsto F[v] :=
 *      \int\limits_{\Omega}
 *        f v
 *      \mathrm{d}x.
 *  @f}
 *  Accordingly the function @f$a@f$ from example @ref exampleintroductionelliptic_pde
 *  "elliptic boundary value problem" induces an operator
 *  @f{align*}{
 *    A: S_h^k &\rightarrow {S_h^k}^{-1}\\
 *    u &\mapsto A(u),
 *  @f}
 *  where @f$A(u)@f$ itself is a functional, defined by
 *  @f{align*}{
 *    A(u): S_h^k &\rightarrow \mathbb{R} \\
 *    v &\mapsto A(u)[v] :=
 *      \int\limits_{\Omega}
 *        a \nabla u \nabla v
 *      \mathrm{d}x.
 *  @f}
 *
 *  @subsection Finite_element_solution_using Finite element solution, using operators and functionals
 *  @anchor remarkintroductionvariational_formulation_functionals_operators
 *  With the notation from definition @ref definitionintroductionfinite_element_solution
 *  "finite element solution" let @f$A@f$ and @f$F@f$ be as in
 *  definition @ref definitionintroductionoperators_functionals
 *  "operator and functional". The finite element
 *  @ref equationintroductionweak_formulation "weak formulation"
 *  can be rewritten as follows.
 *
 *  Find @f${{u_h}_0 \in {S^1_h}_0}@f$, such that
 *  @anchor equationintroductionfinite_element_formulation_functionals_operators
 *  @f{align*}{
 *      A({u_h}_0)[v_h] = F[v_h] - A(g_h)[v_h] &&\text{for all } v_h \in {S^1_h}_0.
 *  @f}
 *  The finite element solution is then given as in definition 
 *  @ref definitionintroductionweak_formulation "weak formulation".
 *
 *  @section Algorithm_For_Solving Algorithm to solve the elliptic boundary value problem
 *
 *  @subsection Step1 Grid
 *  First we have to determine a domain @f$\Omega\subset\mathbb{R}^d@f$,
 *  (lets say @f$d=3@f$)
 *  and a discretization @f$\tau_h@f$ of @f$\Omega@f$.
 *  @code
 *  typedef Dune::GridSelector::GridType
 *    HGridType;
 *
 *  typedef Dune::AdaptiveLeafGridPart< HGridType >
 *    GridPartType;
 *
 *  Dune::GridPtr< HGridType > gridPtr( "macrogrids/unitcube2.dgf" );
 *
 *  GridPartType gridPart( *gridPtr );
 *  @endcode
 *
 *  @subsection Step2 Function space
 *  Define an analytical function space @f$V:\mathbb{R}^d\rightarrow\mathbb{R}@f$
 *  and an analytical function @f$v\in\{\mathbb{R}\rightarrow\mathbb{R}\}@f$.
 *  @code
 *  typedef Dune::FunctionSpace< double, double, HGridType::dimension, 1 >
 *    FunctionSpaceType;
 *
 *  typedef Dune::Function< double, double >
 *    FunctionType;
 *  @endcode
 *
 *  @subsection Step3 Local operations
 *  Define a product operation @f$F:S_h^1\rightarrow\mathbb{R}@f$,
 *  defined by @f$v\mapsto fv@f$, and a finite element operator 
 *  @f${A: {S_{h}^{1}} \rightarrow {S_{h}^{1}}^{-1}}@f$,
 *  defined by @f$(u,v)\mapsto a\nabla u\nabla v@f$.
 *  @code
 *  typedef ProductOperation< FunctionSpaceType >
 *    ProductOperationType;
 *  ProductOperationType productOperation;
 *
 *  typedef EllipticOperation< FunctionSpaceType >
 *    EllipticOperationType;
 *  EllipticOperationType ellipticOperation;
 *  @endcode
 *
 *
 * @subsection Step4 Integration
 * Define an integrator @f$I@f$ for the product operation, which computes
 * @f[I(F(v)):=\int_e fv\mathrm{d}x@f] numerically,
 * and an integrator @f$I'@f$ for the elliptic integrator, which computes
 * @f[I(F(u,v)):= \int_e a\nabla u\nabla v\mathrm{d}x@f] numerically.
 * @code
 * typedef LocalOperation::Integrator::Codim0< FunctionSpaceType, ProductOperationType >
 *   ProductIntegratorType;
 * ProductIntegratorType productIntegrator( productOperation );
 *
 * typedef LocalOperation::Integrator::Codim0< FunctionSpaceType, EllipticOperationType >
 *   EllipticIntegratorType;
 * EllipticIntegratorType ellipticIntegrator( ellipticOperation );
 * @endcode
 *
 * @subsection Step5_1 Discrete function space
 * Define the discrete function space @f$S_h^1@f$.
 * @code
 * typedef Dune::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder >
 *   DiscreteH1Type;
 * DiscreteH1Type discreteH1( gridPart );
 * typedef Dune::AdaptiveDiscreteFunction< DiscreteH1Type >
 *   DiscreteFunctionType;
 * @endcode
 * @subsection Step5_2 Test space
 * Define constraints @f$C=\{C_1,\ldots,C_M\}@f$ and
 * the test space @f${S_h^1}_0\subset S_h^1@f$ as a linear subspace.
 *
 * Here we choose Dirichlet constraints, i.e.
 * @f[
 *  C_i[v_h]=0\quad \forall v_h\in S_h^1 \text{ on }\partial \Omega
 * @f]
 * @code
 * typedef Constraints::Dirichlet< DiscreteH1Type >
 *   DirichletConstraints;
 * DirichletConstraints dirichletConstraints( discreteH1 );
 *
 * typedef Subspace::Linear< DiscreteH1Type, DirichletConstraints >
 *   DiscreteH10Type;
 * DiscreteH10Type discreteH10( discreteH1, dirichletConstraints );
 * @endcode
 *
 * @subsection Step5_3 Ansatz space
 * Define the ansatz space @f${S_h^1}_g\subset S_h^1@f$
 * as an affine subspace.
 * @code
 * //not needed here.
 * @endcode
 *
 *
 * @subsection Step6 operator and functional
 * Define an operator and a functional.
 * @code
 * typedef Operator::Linear< EllipticIntegratorType, DiscreteH1Type >
 *   FEMellipticOperatorType;
 * FEMellipticOperatorType femEllipticOperator( ellipticIntegrator, discreteH1 );
 *
 * typedef Functional::FiniteElementLOP< DiscreteH1Type, ProductIntegratorType >
 *   FEMrhsFunctionalType;
 * FEMrhsFunctionalType femRhsFunctional( discreteH1, productIntegrator );
 * @endcode
 *
 *
 * @subsection Step7 Matrix, rhs and solution storage
 * Define storage objects for storing the dofs.
 * @code
 * typedef Dune::FieldMatrix< double, 1, 1 >
 *   FieldMatrixType;
 *
 * typedef Container::MatrixFactory< Dune::BCRSMatrix< FieldMatrixType > >
 *   MatrixFactoryType;
 *
 * typedef MatrixFactoryType::ContainerType
 *   MatrixContainer;
 *
 * typedef MatrixFactoryType::AutoPtrType
 *   MatrixContainerPtr;
 *
 * typedef Container::VectorFactory< Dune::BlockVector< Dune::FieldVector< double, 1 > > >
 *   VectorFactoryType;
 *
 * typedef VectorFactoryType::ContainerType
 *   VectorContainer;
 *
 * typedef VectorFactoryType::AutoPtrType
 *   VectorContainerPtr;
 *
 * MatrixContainerPtr A  = MatrixFactoryType::create( discreteH10 );
 * VectorContainerPtr F  = VectorFactoryType::create( discreteH10 );
 * VectorContainerPtr u0 = VectorFactoryType::create( discreteH10 );
 * @endcode
 *
 * @subsection Step8 Assembler
 * Define an assember and assemble a matrix with entries
 * @f{align*}{
 *   \label{equation::introduction::matrix_entries}
 *   (A)_{i,j} := A( \varphi_{i} )[ \psi_j ]
 * @f}
 * for all basefunctions of the ansatz-space @f${\varphi_i \in S_{h}^1}@f$ 
 * and all basefunctions of the test-space
 * @f${\psi_j \in S_h^1}@f$.
 * Constrain the algebraic system according to @f${S_h^1}_0@f$.
 * @code
 * typedef Assembler::FiniteElement< MatrixContainer, VectorContainer >
 *   Assembler;
 *
 * Assembler::assembleMatrix( femEllipticOperator, *A );
 * Assembler::applyMatrixConstraints( discreteH10, *A );
 * @endcode
 * Assemble a vector with entries
 * @f{align*}{
 *   (F)_j:=F[\psi_j]
 * @f}
 * for all basefunctions of the test space @f$\psi\in S_h^1@f$.
 * Constrain the vector according to @f${S_h^1}_0@f$.
 * @code
 * Assembler::assembleVector( femRhsFunctional, *F );
 * Assembler::applyVectorConstraints( discreteH10, *F );
 * @endcode
 * 
 *
 * @subsection Step9 Preconditioner and solver
 * Define a preconditioner and a solver.
 * @code
 * typedef Dune::MatrixAdapter< MatrixContainer, VectorContainer, VectorContainer >
 *   MatrixAdapterType;
 * MatrixAdapterType op( *A );
 *
 * typedef Dune::SeqILU0< MatrixContainer, VectorContainer, VectorContainer, 1 >
 *   SeqILU0Type;
 * SeqILU0Type prec( *A, 1.0 );
 *
 * typedef Dune::CGSolver< VectorContainer >
 *   CG;
 * @endcode
 * Solve the algebraic system:
 *
 * Find @f${u_h}_0\in S_h^1@f$ such that
 * @f{align*}{
 *   A({u_h}_0)&=F.
 * @f}
 * @code
 * CG cg( op, prec, 1e-4, 100, 2 );
 *
 * Dune::InverseOperatorResult res;
 * cg.apply( *u0, *F, res );
 * @endcode
 * 
 * Compute the solution @f$u_h\in{S_h^1}_g@f$ as 
 * @f[u_h:={u_h}_0+g_h@f].
 * @code
 * //not needed here
 * @endcode
 *
 *
 *
 *
 *
 *
 *
 *
 */

//! Contains several solvers.
namespace Solver
{

/**
 * @brief 
 *
 * @tparam MatrixImp Type of the matrix we want to assemble.
 * @tparam VectorImp Type of the Vector we want to assemble.
 */
template< class MatrixImp, class VectorImp >
class FEMAssembler
{
public:
  //! Type of the matrix.
  typedef MatrixImp
    MatrixType;

  //! Type of the vector.
  typedef VectorImp
    VectorType;

  //! Type of the vector, i.e. double.
  typedef typename VectorType::field_type
    FieldType;

  //! Type of the local matrix.
  typedef Dune::Functionals::Common::LocalMatrix< FieldType >
    LocalMatrixType;

  //! Type of the local vector.
  typedef Dune::Functionals::Common::LocalVector< FieldType >
    LocalVectorType;

  /**
   * @brief Assembles the system matrix.
   *
   * @param op The operator.
   * @param m The matrix.
   */
  template< class Operator >
  static void assembleMatrix( const Operator& op, MatrixType& m )
  {
    typedef typename Operator::DiscreteFunctionSpaceType
      DFSType;
    typedef typename DFSType::BaseFunctionSetType
      BFSType;
    typedef typename DFSType::IteratorType
      ItType;
    typedef typename ItType::Entity
      EntityType;

    const DFSType& space = op.space();
    // check that number of dofs in space is equal to matrix size
    // assert()
    ItType it = space.begin();
    for( ; it != space.end(); ++it )
    {
      const EntityType &en = *it;
      const BFSType& bfs = space.baseFunctionSet( en );
      LocalMatrixType localMatrix( bfs.numBaseFunctions(), bfs.numBaseFunctions() );
      localMatrix = op.applyLocal( en );

      addToMatrix( space, localMatrix, en, m );
    }
  }

  template< class Functional >
  static void assembleVector( const Functional& functional, VectorType& vec )
  {
    typedef typename Functional::DiscreteFunctionSpaceType
      DFSType;
//    typedef typename DFSType::BaseFunctionSetType
//      BFSType;
    typedef typename DFSType::IteratorType
      ItType;
    typedef typename ItType::Entity
      EntityType;

    const DFSType& space = functional.space();

    // check that number of dofs in space is equal to vector size
    assert( space.size() == (int)vec.size() );
    ItType it = space.begin();
    for( ; it != space.end(); ++it )
    {
      const EntityType &en = *it;
//      const BFSType& bfs = space.baseFunctionSet( en );

//      LocalVectorType localVector = functional.applyLocal( en, bfs );
      LocalVectorType localVector = functional.applyLocal( en );

      addToVector( space, localVector, en, vec );
    }
  }


  /// \todo merge later with assembleMatrix
  /// \todo implement a PrecompiledConstraints class which wraps an existing
  ///       Constraints class for efficiency at the cost of one grid walk
  template< class ConstrainedDFS >
  static void applyMatrixConstraints( const ConstrainedDFS& cSpace,
                                      MatrixType& m )
  {
    typedef typename ConstrainedDFS::ConstraintsType
      ConstraintsType;
    typedef typename ConstrainedDFS::DiscreteFunctionSpaceType
      DFSType;
    typedef typename ConstraintsType::LocalConstraintsType
      LocalConstraintsType;
    typedef typename DFSType::IteratorType
      ItType;
    typedef typename ItType::Entity
      EntityType;

    const ConstraintsType& constraints = cSpace.constraints();

    // check that number of dofs in space is equal to matrix size
    // assert()

    ItType it  = cSpace.space().begin();
    ItType end = cSpace.space().end();

    for( ; it != end; ++it )
    {
      const EntityType& en = *it;

      const LocalConstraintsType localConstraints = constraints.local( en );

      setLocalConstraintsInMatrix( cSpace, localConstraints, en, m );
    }

  }

  // @todo implementation
  template< class ConstrainedDFS >
  static void applyVectorConstraints( const ConstrainedDFS& cSpace,
                                      VectorType& vec )
  {
    typedef typename ConstrainedDFS::ConstraintsType
      ConstraintsType;
    typedef typename ConstrainedDFS::DiscreteFunctionSpaceType
      DFSType;
    typedef typename ConstraintsType::LocalConstraintsType
      LocalConstraintsType;
    typedef typename DFSType::IteratorType
      ItType;
    typedef typename ItType::Entity
      EntityType;

    const ConstraintsType & constraints = cSpace.constraints();

    ItType it  = cSpace.space().begin();
    ItType end = cSpace.space().end();

    for( ; it != end; ++it )
    {
      const EntityType& en = *it;

      const LocalConstraintsType localConstraints = constraints.local( en );

      for (unsigned int i = 0; i < localConstraints.rowDofsSize(); ++i) {
        vec[localConstraints.rowDofs(i)] = 0;
      }
    }
  }

public:
  /// \todo move to matrixContainer factory
  template< class DFSType,
            class Entity >
  static void addToMatrix( const DFSType& space,
                           const LocalMatrixType& localMatrix,
                           const Entity& en,
                           MatrixType &m )
  {
    for( unsigned int i = 0; i < localMatrix.N(); i++ )
    {
      for( unsigned int j = 0; j < localMatrix.M(); j++ )
      {
        const int globalI = space.mapToGlobal( en, i );
        const int globalJ = space.mapToGlobal( en, j );

        m[globalI][globalJ] += localMatrix.get( i, j );
      }
    }
  }

  template< class DFSType,
            class Entity >
  static void addToVector( const DFSType& space,
                           const LocalVectorType& localVector,
                           const Entity& entity,
                           VectorType& vector )
  {
    for( unsigned int ii = 0; ii < localVector.size(); ++ii )
    {
      const int globalI = space.mapToGlobal( entity, ii );

      vector[globalI] += localVector[ii];
    }
  }

  /// \todo move to Constraints class
  template< class ConstrainedDFS,
            class Entity >
  static void setLocalConstraintsInMatrix( const ConstrainedDFS& cSpace,
                                           const typename ConstrainedDFS::ConstraintsType::LocalConstraintsType& localConstraints,
                                           const Entity& en,
                                           MatrixType &m )
  {

    for( unsigned int i = 0; i < localConstraints.rowDofsSize(); i++ )
    {
      for( unsigned int j = 0; j < localConstraints.columnDofsSize(); j++ )
      {
        m[localConstraints.rowDofs(i)][localConstraints.columnDofs(j)]
          = localConstraints.localMatrix(i,j);
      }
    }
  }


}; // end of class FEMAssembler

} // end of namespace Solver

/**
 * @brief Contains various assemblers.
 *
 * @todo The namespace is inserted into the wrong directory/file.
 */
namespace Assembler
{

template< class MatrixImp, class VectorImp >
class FiniteElement
{
public:
  typedef MatrixImp
    MatrixType;

  typedef VectorImp
    VectorType;

  typedef typename VectorType::field_type
    FieldType;

  typedef Dune::Functionals::Common::LocalMatrix< FieldType >
    LocalMatrixType;

  typedef Dune::Functionals::Common::LocalVector< FieldType >
    LocalVectorType;

  template< class OperatorType >
  static void assembleMatrix( const OperatorType& op, MatrixType& matrix )
  {
    // some types
    typedef typename OperatorType::DiscreteAnsatzFunctionSpaceType
      DiscreteAnsatzFunctionSpaceType;
    typedef typename OperatorType::DiscreteTestFunctionSpaceType
      DiscreteTestFunctionSpaceType;
    typedef typename DiscreteAnsatzFunctionSpaceType::IteratorType
      EntityIteratorType;
    typedef typename EntityIteratorType::Entity
      EntityType;

    // discrete function spaces
    const DiscreteAnsatzFunctionSpaceType& ansatzSpace = op.ansatzSpace();
    const DiscreteTestFunctionSpaceType& testSpace = op.testSpace();

    // gridwalk
    const EntityIteratorType behindLastEntity = ansatzSpace.end();
    for( EntityIteratorType entityIterator = ansatzSpace.begin(); entityIterator != behindLastEntity; ++entityIterator )
    {
      const EntityType& entity = *entityIterator;
      LocalMatrixType localMatrix(  ansatzSpace.baseFunctionSet( entity ).numBaseFunctions(),
                                    testSpace.baseFunctionSet( entity ).numBaseFunctions() );

      localMatrix = op.applyLocal( entity );

      addToMatrix( ansatzSpace, testSpace, entity, localMatrix, matrix );
    }
  }

  template< class Functional >
  static void assembleVector( const Functional& functional, VectorType& vec )
  {
    typedef typename Functional::DiscreteFunctionSpaceType
      DFSType;
//    typedef typename DFSType::BaseFunctionSetType
//      BFSType;
    typedef typename DFSType::IteratorType
      ItType;
    typedef typename ItType::Entity
      EntityType;

    const DFSType& space = functional.space();

    // check that number of dofs in space is equal to vector size
    assert( space.size() == (int)vec.size() );
    ItType it = space.begin();
    for( ; it != space.end(); ++it )
    {
      const EntityType &en = *it;
//      const BFSType& bfs = space.baseFunctionSet( en );

//      LocalVectorType localVector = functional.applyLocal( en, bfs );
      LocalVectorType localVector = functional.applyLocal( en );

      addToVector( space, localVector, en, vec );
    }
  }


  /// \todo merge later with assembleMatrix
  /// \todo implement a PrecompiledConstraints class which wraps an existing
  ///       Constraints class for efficiency at the cost of one grid walk
  template< class ConstrainedDFS >
  static void applyMatrixConstraints( const ConstrainedDFS& cSpace,
                                      MatrixType& m )
  {
    typedef typename ConstrainedDFS::ConstraintsType
      ConstraintsType;
    typedef typename ConstrainedDFS::DiscreteFunctionSpaceType
      DFSType;
    typedef typename ConstraintsType::LocalConstraintsType
      LocalConstraintsType;
    typedef typename DFSType::IteratorType
      ItType;
    typedef typename ItType::Entity
      EntityType;

    const ConstraintsType& constraints = cSpace.constraints();

    // check that number of dofs in space is equal to matrix size
    // assert()

    ItType it  = cSpace.space().begin();
    ItType end = cSpace.space().end();

    for( ; it != end; ++it )
    {
      const EntityType& en = *it;

      const LocalConstraintsType localConstraints = constraints.local( en );

      setLocalConstraintsInMatrix( cSpace, localConstraints, en, m );
    }

  }

  // @todo implementation
  template< class ConstrainedDFS >
  static void applyVectorConstraints( const ConstrainedDFS& cSpace,
                                      VectorType& vec )
  {
    typedef typename ConstrainedDFS::ConstraintsType
      ConstraintsType;
    typedef typename ConstrainedDFS::DiscreteFunctionSpaceType
      DFSType;
    typedef typename ConstraintsType::LocalConstraintsType
      LocalConstraintsType;
    typedef typename DFSType::IteratorType
      ItType;
    typedef typename ItType::Entity
      EntityType;

    const ConstraintsType & constraints = cSpace.constraints();

    ItType it  = cSpace.space().begin();
    ItType end = cSpace.space().end();

    for( ; it != end; ++it )
    {
      const EntityType& en = *it;

      const LocalConstraintsType localConstraints = constraints.local( en );

      for (unsigned int i = 0; i < localConstraints.rowDofsSize(); ++i) {
        vec[localConstraints.rowDofs(i)] = 0;
      }
    }
  }

public:
  /// \todo move to matrixContainer factory
  template< class DiscreteAnsatzFunctionSpaceType,
            class DiscreteTestFunctionSpaceType,
            class EntityType >
  static void addToMatrix( const DiscreteAnsatzFunctionSpaceType& ansatzSpace,
                           const DiscreteTestFunctionSpaceType& testSpace,
                           const EntityType& en,
                           const LocalMatrixType& localMatrix,
                           MatrixType& matrix )
  {
    for( unsigned int i = 0; i < localMatrix.N(); i++ )
    {
      for( unsigned int j = 0; j < localMatrix.M(); j++ )
      {
        const int globalI = ansatzSpace.mapToGlobal( en, i );
        const int globalJ = testSpace.mapToGlobal( en, j );

        matrix[globalI][globalJ] += localMatrix[i][j];
      }
    }
  }

  template< class DFSType,
            class Entity >
  static void addToVector( const DFSType& space,
                           const LocalVectorType& localVector,
                           const Entity& entity,
                           VectorType& vector )
  {
    for( unsigned int ii = 0; ii < localVector.size(); ++ii )
    {
      const int globalI = space.mapToGlobal( entity, ii );

      vector[globalI] += localVector[ii];
    }
  }

  /// \todo move to Constraints class
  template< class ConstrainedDFS,
            class Entity >
  static void setLocalConstraintsInMatrix( const ConstrainedDFS& cSpace,
                                           const typename ConstrainedDFS::ConstraintsType::LocalConstraintsType& localConstraints,
                                           const Entity& en,
                                           MatrixType &m )
  {

    for( unsigned int i = 0; i < localConstraints.rowDofsSize(); i++ )
    {
      for( unsigned int j = 0; j < localConstraints.columnDofsSize(); j++ )
      {
        m[localConstraints.rowDofs(i)][localConstraints.columnDofs(j)]
          = localConstraints.localMatrix(i,j);
      }
    }
  }


}; // end of class FiniteElement

} // end of namespace Assembler

} // end of namespace Functionals

} // end of namespace Dune


#endif /* end of include guard: DUNE_FEM_FUNCTIONALS_SOLVER_FEMASSEMBLER_HH */
