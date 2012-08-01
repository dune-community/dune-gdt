#ifndef DUNE_DETAILED_DISCRETIZATIONS_LA_BACKEND_SOLVER_EIGEN_HH
#define DUNE_DETAILED_DISCRETIZATIONS_LA_BACKEND_SOLVER_EIGEN_HH

//#ifdef HAVE_EIGEN

#include <Eigen/Eigen>

namespace Dune {

namespace Detailed {

namespace Discretizations {

namespace LA {

namespace Solver {

namespace Eigen {

struct BicgstabIlut
{
public:
  static const std::string id;

  template <class MatrixType, class SolutionType, class RhsType>
  static void apply(MatrixType& matrix, SolutionType& solution, RhsType& rhs, unsigned int maxIter = 5000,
                    double precision = 1e-12)
  {
    typedef typename MatrixType::EntryType EntryType;
    typedef typename MatrixType::StorageType EigenMatrixType;
    typedef ::Eigen::IncompleteLUT<EntryType> PreconditionerType;
    typedef ::Eigen::BiCGSTAB<EigenMatrixType, PreconditionerType> SolverType;
    SolverType solver;
    solver.setMaxIterations(maxIter);
    solver.setTolerance(precision);
    solver.compute(*(matrix.storage()));
    *(solution.storage()) = solver.solve(*(rhs.storage()));
  }
}; // struct BicgstabIlut

const std::string BicgstabIlut::id = "eigen.bicgstab.incompletelut";

struct BicgstabDiagonal
{
public:
  static const std::string id;

  template <class MatrixType, class SolutionType, class RhsType>
  static void apply(MatrixType& matrix, SolutionType& solution, RhsType& rhs, unsigned int maxIter = 5000,
                    double precision = 1e-12)
  {
    typedef typename MatrixType::EntryType EntryType;
    typedef typename MatrixType::StorageType EigenMatrixType;
    typedef ::Eigen::DiagonalPreconditioner<EntryType> PreconditionerType;
    typedef ::Eigen::BiCGSTAB<EigenMatrixType, PreconditionerType> SolverType;
    SolverType solver;
    solver.setMaxIterations(maxIter);
    solver.setTolerance(precision);
    solver.compute(*(matrix.storage()));
    *(solution.storage()) = solver.solve(*(rhs.storage()));
  }
}; // struct BicgstabDiagonal

const std::string BicgstabDiagonal::id = "eigen.bicgstab.diagonal";

struct CgDiagonalUpper
{
public:
  static const std::string id;

  template <class MatrixType, class SolutionType, class RhsType>
  static void apply(MatrixType& matrix, SolutionType& solution, RhsType& rhs, unsigned int maxIter = 5000,
                    double precision = 1e-12)
  {
    typedef typename MatrixType::EntryType EntryType;
    typedef typename MatrixType::StorageType EigenMatrixType;
    typedef ::Eigen::DiagonalPreconditioner<EntryType> PreconditionerType;
    typedef ::Eigen::ConjugateGradient<EigenMatrixType, ::Eigen::Upper, PreconditionerType> SolverType;
    SolverType solver;
    solver.setMaxIterations(maxIter);
    solver.setTolerance(precision);
    solver.compute(*(matrix.storage()));
    *(solution.storage()) = solver.solve(*(rhs.storage()));
  }
}; // struct CgDiagonalUpper

const std::string CgDiagonalUpper::id = "eigen.cg.diagonal.upper";

struct CgDiagonalLower
{
public:
  static const std::string id;

  template <class MatrixType, class SolutionType, class RhsType>
  static void apply(MatrixType& matrix, SolutionType& solution, RhsType& rhs, unsigned int maxIter = 5000,
                    double precision = 1e-12)
  {
    typedef typename MatrixType::EntryType EntryType;
    typedef typename MatrixType::StorageType EigenMatrixType;
    typedef ::Eigen::DiagonalPreconditioner<EntryType> PreconditionerType;
    typedef ::Eigen::ConjugateGradient<EigenMatrixType, ::Eigen::Lower, PreconditionerType> SolverType;
    SolverType solver;
    solver.setMaxIterations(maxIter);
    solver.setTolerance(precision);
    solver.compute(*(matrix.storage()));
    *(solution.storage()) = solver.solve(*(rhs.storage()));
  }
}; // struct CgDiagonalLower

const std::string CgDiagonalLower::id = "eigen.cg.diagonal.lower";

struct SimplicialcholeskyUpper
{
public:
  static const std::string id;

  template <class MatrixType, class SolutionType, class RhsType>
  static void apply(MatrixType& matrix, SolutionType& solution, RhsType& rhs, unsigned int maxIter = 5000,
                    double precision = 1e-12)
  {
    typedef typename MatrixType::EntryType EntryType;
    typedef typename MatrixType::StorageType EigenMatrixType;
    typedef ::Eigen::SimplicialCholesky<EigenMatrixType, ::Eigen::Upper> SolverType;
    SolverType solver;
    solver.compute(*(matrix.storage()));
    *(solution.storage()) = solver.solve(*(rhs.storage()));
  }
}; // struct SimplicialcholeskyUpper

const std::string SimplicialcholeskyUpper::id = "eigen.simplicialcholesky.upper";

struct SimplicialcholeskyLower
{
public:
  static const std::string id;

  template <class MatrixType, class SolutionType, class RhsType>
  static void apply(MatrixType& matrix, SolutionType& solution, RhsType& rhs, unsigned int maxIter = 5000,
                    double precision = 1e-12)
  {
    typedef typename MatrixType::EntryType EntryType;
    typedef typename MatrixType::StorageType EigenMatrixType;
    typedef ::Eigen::SimplicialCholesky<EigenMatrixType, ::Eigen::Lower> SolverType;
    SolverType solver;
    solver.compute(*(matrix.storage()));
    *(solution.storage()) = solver.solve(*(rhs.storage()));
  }
}; // struct SimplicialcholeskyLower

const std::string SimplicialcholeskyLower::id = "eigen.simplicialcholesky.lower";

} // namespace Eigen

} // namespace Solver

} // namespace LA

} // namespace Discretizations

} // namespace Detailed

} // namespace Dune

//#endif // HAVE_EIGEN

#endif // DUNE_DETAILED_DISCRETIZATIONS_LA_BACKEND_SOLVER_EIGEN_HH
