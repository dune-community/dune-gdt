#ifndef DUNE_DETAILED_DISCRETIZATIONS_LA_BACKEND_SOLVER_EIGEN_HH
#define DUNE_DETAILED_DISCRETIZATIONS_LA_BACKEND_SOLVER_EIGEN_HH

//#ifdef HAVE_EIGEN

#include <Eigen/Eigen>

namespace Dune
{

namespace DetailedDiscretizations
{

namespace LA
{

namespace Solver
{

namespace Eigen
{

struct BicgstabIlut
{
public:
  static const std::string id;

  template< class MatrixType, class SolutionType, class RhsType >
  static void apply(MatrixType& matrix, SolutionType& solution, RhsType& rhs, unsigned int maxIter, double precision)
  {
    typedef typename MatrixType::EntryType EntryType;
    typedef typename MatrixType::StorageType EigenMatrixType;
    typedef ::Eigen::IncompleteLUT< EntryType > PreconditionerType;
    typedef ::Eigen::BiCGSTAB< EigenMatrixType, PreconditionerType > SolverType;
    SolverType solver;
    solver.setMaxIterations(maxIter);
    solver.setTolerance(precision);
    solver.compute(*(matrix.storage()));
    *(solution.storage()) = solver.solve(*(rhs.storage()));
  }
}; // struct BicgstabIlut

const std::string BicgstabIlut::id = "eigen.bicgstab.incompletelut";

struct CgDiagonalUp
{
public:
  static const std::string id;

  template< class MatrixType, class SolutionType, class RhsType >
  static void apply(MatrixType& matrix, SolutionType& solution, RhsType& rhs, unsigned int maxIter, double precision)
  {
    typedef typename MatrixType::EntryType EntryType;
    typedef typename MatrixType::StorageType EigenMatrixType;
    typedef ::Eigen::DiagonalPreconditioner< EntryType > PreconditionerType;
    typedef ::Eigen::ConjugateGradient< EigenMatrixType, ::Eigen::Upper, PreconditionerType > SolverType;
    SolverType solver;
    solver.setMaxIterations(maxIter);
    solver.setTolerance(precision);
    solver.compute(*(matrix.storage()));
    *(solution.storage()) = solver.solve(*(rhs.storage()));
  }
}; // struct CgDiagonalUp

const std::string CgDiagonalUp::id = "eigen.cg.diagonal.up";

} // namespace Eigen

} // namespace Solver

} // namespace LA

} // DetailedDiscretizations

} // namespace Dune

//#endif // HAVE_EIGEN

#endif // DUNE_DETAILED_DISCRETIZATIONS_LA_BACKEND_SOLVER_EIGEN_HH
