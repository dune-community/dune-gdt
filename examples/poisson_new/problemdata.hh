#ifndef PROBLEM_HH
#define PROBLEM_HH

#include <cassert>
#include <cmath>

#include <dune/common/array.hh>

#include <dune/fem-howto/probleminterfaces.hh>

namespace Dune {

//! sine product problem
template <int dim>
struct SinProduct
{
  typedef FunctionSpace<double, double, dim, 1> FunctionSpaceType;

  static const int dimDomain = FunctionSpaceType::dimDomain;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

  void evaluate(const DomainType& x, RangeType& ret) const
  {
    ret = RangeFieldType(1);
    for (int i = 0; i < dimDomain; ++i)
      ret *= sin(2. * M_PI * x[i]);
  }

  void jacobian(const DomainType& x, JacobianRangeType& grad) const
  {
    array<RangeFieldType, dimDomain> costbl, sintbl;
    for (int i = 0; i < dimDomain; ++i) {
      costbl[i] = cos(2. * M_PI * x[i]);
      sintbl[i] = sin(2. * M_PI * x[i]);
    }

    for (int i = 0; i < dimDomain; ++i) {
      grad[0][i] = 2. * M_PI;
      for (int l = 0; l < dimDomain; ++l)
        grad[0][i] *= (l == i ? costbl[l] : sintbl[l]);
    }
  }

  void hessian(const DomainType& x, HessianRangeType& hessian) const
  {
    array<RangeFieldType, dimDomain> costbl, sintbl;
    for (int i = 0; i < dimDomain; ++i) {
      costbl[i] = cos(2. * M_PI * x[i]);
      sintbl[i] = sin(2. * M_PI * x[i]);
    }

    for (int i = 0; i < dimDomain; ++i) {
      for (int j = 0; j < dimDomain; ++j) {
        if (i != j) {
          hessian[0][i][j] = 4. * M_PI * M_PI;
          for (int l = 0; l < dimDomain; ++l)
            hessian[0][i][j] *= ((l == i) || (l == j) ? costbl[l] : sintbl[l]);
        } else {
          hessian[0][i][j] = -4. * M_PI * M_PI;
          for (int l = 0; l < dimDomain; ++l)
            hessian[0][i][j] *= sintbl[l];
        }
      }
    }
  }
};


//! cosine product problem
template <int dim>
struct CosProduct
{
  typedef FunctionSpace<double, double, dim, 1> FunctionSpaceType;

  static const int dimDomain = FunctionSpaceType::dimDomain;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

  //! the exact solution
  void evaluate(const DomainType& x, RangeType& ret) const
  {
    ret = 1;
    for (int i = 0; i < dimDomain; ++i)
      ret *= cos(2. * M_PI * x[i]);
  }

  //! the gradient of the exact solution
  void jacobian(const DomainType& x, JacobianRangeType& grad) const
  {
    array<RangeFieldType, dimDomain> costbl;
    array<RangeFieldType, dimDomain> sintbl;
    for (int i = 0; i < dimDomain; ++i) {
      costbl[i] = cos(2. * M_PI * x[i]);
      sintbl[i] = sin(2. * M_PI * x[i]);
    }

    for (int i = 0; i < dimDomain; ++i) {
      grad[0][i] = -2. * M_PI;
      for (int j = 0; j < dimDomain; ++j)
        grad[0][i] *= (i == j ? sintbl[j] : costbl[j]);
    }
  }

  void hessian(const DomainType& x, HessianRangeType& hessian) const
  {
    array<RangeFieldType, dimDomain> costbl;
    array<RangeFieldType, dimDomain> sintbl;
    for (int i = 0; i < dimDomain; ++i) {
      costbl[i] = cos(2. * M_PI * x[i]);
      sintbl[i] = sin(2. * M_PI * x[i]);
    }

    for (int i = 0; i < dimDomain; ++i) {
      for (int j = 0; j < dimDomain; ++j) {
        if (i != j) {
          hessian[0][i][j] = 4. * M_PI * M_PI;
          for (int l = 0; l < dimDomain; ++l)
            hessian[0][i][j] *= ((l == i) || (l == j) ? sintbl[l] : costbl[l]);
        } else {
          hessian[0][i][j] = -4. * M_PI * M_PI;
          for (int l = 0; l < dimDomain; ++l)
            hessian[0][i][j] *= costbl[l];
        }
      }
    }
  }
};


//! test problem from ALBERTA
template <int dim>
struct AlbertaSolution
{
  typedef FunctionSpace<double, double, dim, 1> FunctionSpaceType;

  static const int dimDomain = FunctionSpaceType::dimDomain;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

  //! the exact solution
  void evaluate(const DomainType& x, RangeType& ret) const
  {
    // use scalar product implemented on FieldVector
    const RangeFieldType xsqr = x * x;
    ret                       = exp(-10.0 * xsqr);
  }

  //! the gradient of the exact solution
  void jacobian(const DomainType& x, JacobianRangeType& grad) const
  {
    // set components of gradient
    grad[0] = x;
    // multiply with factor
    grad[0] *= -20.0 * exp(-10.0 * (x * x));
  }

  void hessian(const DomainType& x, HessianRangeType& hessian) const
  {
    RangeType ux;
    evaluate(x, ux);

    array<RangeFieldType, dimDomain> dx;
    for (int i = 0; i < dimDomain; ++i) {
      dx[i] = -20. * x[i];
    }

    for (int i = 0; i < dimDomain; ++i) {
      for (int j = 0; j < dimDomain; ++j)
        hessian[0][i][j] = ux[0] * (dx[i] * dx[j]);
      hessian[0][i][i] -= 20.0 * ux[0];
    }
  }
};


//! The famous Corner problem (L-shape problem)
template <int dim>
struct CornerSolution
{
  typedef FunctionSpace<double, double, dim, 1> FunctionSpaceType;

  static const int dimDomain = FunctionSpaceType::dimDomain;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

  //! the exact solution
  void evaluate(const DomainType& x, RangeType& ret) const
  {
    const double r   = x.two_norm();
    const double fac = 2.0 / 3.0;

    if (r > 1e-8) {
      const double phi = (x[1] >= 0 ? acos(x[0] / r) : acos(-x[0] / r) + M_PI);

      ret = sin(fac * phi);

      ret *= pow(r, fac);
    } else
      ret = 0;
  }

  //! the gradient of the exact solution
  void jacobian(const DomainType& x, JacobianRangeType& grad) const
  {
    // only implemented for 2d
    assert(dimDomain == 2);

    const double r   = x.two_norm();
    const double fac = 2.0 / 3.0;

    const double phi = (x[1] >= 0 ? acos(x[0] / r) : acos(-x[0] / r) + M_PI);

    grad[0][0] = sin(fac * phi) * x[0] - cos(fac * phi) * x[1];
    grad[0][1] = sin(fac * phi) * x[1] + cos(fac * phi) * x[0];

    grad *= fac * pow(r, -2 * fac);
  }

  void hessian(const DomainType& x, HessianRangeType& hessian) const
  {
    // only implemented for 2d
    assert(dimDomain == 2);

    /*
    const double r = x.two_norm();
    const double fac = 2.0 / 3.0;

    const double phi
      = (x[1] >= 0 ? acos( x[0] / r ) : acos( -x[0] / r ) + M_PI);

    array< RangeFieldType, dimDomain > costbl, sintbl;
    for( int i = 0; i < dimDomain; ++i )
    {
      costbl[i] = cos( 2.0 * M_PI * x[i] );
      sintbl[i] = sin( 2.0 * M_PI * x[i] );
    }
    */

    for (int i = 0; i < dimDomain; ++i) {
      for (int j = 0; j < dimDomain; ++j) {
        hessian[0][i][j] = 0.;
        /*
         * should work for the momement, an explicite implementation of the
         * hessian matrix can be done later
        if( i != j )
        {
        }
        else
        {
        }
        */
      }
    }
  }
};


//! The rotating anisotropy problem
template <class GridImp, class U>
class RotatingAnisotropyProblem : public ProblemInterface<typename U::FunctionSpaceType>
{
  typedef ProblemInterface<typename U::FunctionSpaceType> BaseType;

public:
  typedef typename BaseType::FunctionSpaceType FunctionSpaceType;

  static const int dimRange  = FunctionSpaceType::dimRange;
  static const int dimDomain = FunctionSpaceType::dimDomain;

  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;
  typedef typename BaseType::DomainFieldType DomainFieldType;
  typedef typename BaseType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

  typedef typename BaseType::DiffusionMatrixType DiffusionMatrixType;

  typedef array<DiffusionMatrixType, dimDomain> DivergenceDiffunsionMatrixType;

  double delta_;

  explicit RotatingAnisotropyProblem(const U& u = U())
    : delta_(1E-3)
    , u_(u)
  {
  }

  //! the right hand side
  void f(const DomainType& x, RangeType& ret) const
  {
    if (dimDomain != 2)
      DUNE_THROW(InvalidStateException, "Problem 'rotating anisotropy' may only be used in 2D");

    HessianRangeType Hu;
    u_.hessian(x, Hu);

    DiffusionMatrixType Mk;
    K(x, Mk);

    DivergenceDiffunsionMatrixType divMk;
    DivK(x, divMk);

    JacobianRangeType grad;
    u_.jacobian(x, grad);

    for (int j = 0; j < dimRange; ++j) {
      ret[j] = RangeFieldType(0);
      for (int s = 0; s < dimDomain; ++s) {
        for (int l = 0; l < dimDomain; ++l) {
          ret[j] -= Hu[j][s][l] * Mk[s][l];
          ret[j] -= divMk[s][s][l] * grad[j][l];
        }
      }
    }
  }

  //! the exact solution
  void u(const DomainType& x, RangeType& ret) const
  {
    if (dimDomain != 2)
      DUNE_THROW(InvalidStateException, "Problem 'rotating anisotropy' may only be used in 2D");

    u_.evaluate(x, ret);
  }

  //! the diffusion matrix
  void K(const DomainType& arg, DiffusionMatrixType& k) const
  {
    double x  = arg[0];
    double y  = arg[1];
    double rt = x * x + y * y;
    k[0][0]   = (y * y + delta_ * x * x) / rt;
    k[1][1]   = (x * x + delta_ * y * y) / rt;
    k[1][0] = k[0][1] = -(1 - delta_) * x * y / rt;
  }

  //! divergence of the diffusion matrix
  void DivK(const DomainType& arg, DivergenceDiffunsionMatrixType& divK) const
  {
    double x         = arg[0];
    double y         = arg[1];
    double rt        = x * x + y * y;
    const double eps = 1 - delta_;

    divK[0][0][0] = eps * ((-2. * x * y * y) / (rt * rt));
    divK[0][1][1] = eps * ((2. * x * x * y) / (rt * rt));
    divK[0][0][1] = divK[0][1][0] = eps * ((-1. * y * (y * y - x * x)) / (rt * rt));
    divK[1][0][0] = eps * ((2. * x * y * y) / (rt * rt));
    divK[1][1][1] = eps * ((-2. * y * x * x) / (rt * rt));
    divK[1][0][1] = divK[1][1][0] = eps * ((-1. * x * (x * x - y * y)) / (rt * rt));
  }

  bool constantK() const
  {
    return false;
  }

  //! the gradient of the exact solution
  virtual void gradient(const DomainType& x, JacobianRangeType& grad) const
  {
    if (dimDomain != 2)
      DUNE_THROW(InvalidStateException, "Problem 'rotating anisotropy' may only be used in 2D");

    u_.jacobian(x, grad);
  }

private:
  U u_;
};


template <class GridImp, class U>
class LaplaceProblem : public ProblemInterface<typename U::FunctionSpaceType>
{
  typedef ProblemInterface<typename U::FunctionSpaceType> BaseType;

public:
  typedef typename BaseType::FunctionSpaceType FunctionSpaceType;

  static const int dimRange  = FunctionSpaceType::dimRange;
  static const int dimDomain = FunctionSpaceType::dimDomain;

  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;
  typedef typename BaseType::DomainFieldType DomainFieldType;
  typedef typename BaseType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

  typedef typename BaseType::DiffusionMatrixType DiffusionMatrixType;

  explicit LaplaceProblem(const U& u = U())
    : u_(u)
  {
  }

  //! the right hand side (i.e., the Laplace of u)
  void f(const DomainType& x, RangeType& ret) const
  {
    HessianRangeType Hu;
    u_.hessian(x, Hu);
    for (int j = 0; j < dimRange; ++j) {
      ret[j] = RangeFieldType(0);
      for (int k = 0; k < dimDomain; ++k)
        ret[j] -= Hu[j][k][k];
    }
  }

  //! the exact solution
  void u(const DomainType& x, RangeType& ret) const
  {
    u_.evaluate(x, ret);
  }

  //! the diffusion matrix
  void K(const DomainType& x, DiffusionMatrixType& m) const
  {
    m = 0;
    for (int i = 0; i < dimDomain; ++i)
      m[i][i] = 1;
  }

  bool constantK() const
  {
    return true;
  }

  //! the gradient of the exact solution
  virtual void gradient(const DomainType& x, JacobianRangeType& grad) const
  {
    u_.jacobian(x, grad);
  }

private:
  U u_;
};


template <class GridImp, class U>
class SphereProblem : public ProblemInterface<typename U::FunctionSpaceType>
{
  typedef ProblemInterface<typename U::FunctionSpaceType> BaseType;

public:
  typedef typename BaseType::FunctionSpaceType FunctionSpaceType;

  static const int dimRange  = FunctionSpaceType::dimRange;
  static const int dimDomain = FunctionSpaceType::dimDomain;

  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::JacobianRangeType JacobianRangeType;
  typedef typename BaseType::DomainFieldType DomainFieldType;
  typedef typename BaseType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

  typedef typename BaseType::DiffusionMatrixType DiffusionMatrixType;

  explicit SphereProblem(const U& u = U())
    : u_(u)
  {
  }

  //! the right hand side (i.e., the Laplace-Beltrami of u)
  void f(const DomainType& x, RangeType& ret) const
  {
    DomainType y = x;
    y /= x.two_norm();
    DomainType p = y;
    p[0] -= 1.; // lets have the origin on the sphere...

    JacobianRangeType Ju;
    u_.jacobian(p, Ju);
    Ju.mv(y, ret);
    ret *= dimDomain - 1;

    HessianRangeType Hu;
    u_.hessian(p, Hu);
    for (int i = 0; i < dimDomain; ++i) {
      DomainType ds = y;
      ds *= -y[i];
      ds[i] += 1;

      for (int j = 0; j < dimRange; ++j) {
        DomainType Hds;
        Hu[j].mv(ds, Hds);
        ret[j] -= ds * Hds;
      }
    }
  }

  //! the exact solution
  void u(const DomainType& x, RangeType& ret) const
  {
    DomainType y = x;
    y /= x.two_norm();
    DomainType p = y;
    p[0] -= 1.; // lets have the origin on the sphere...
    u_.evaluate(p, ret);
  }

  //! the diffusion matrix
  void K(const DomainType& x, DiffusionMatrixType& m) const
  {
    m = 0;
    for (int i = 0; i < dimDomain; ++i)
      m[i][i] = 1;
  }

  bool constantK() const
  {
    return true;
  }

  //! the gradient of the exact solution
  virtual void gradient(const DomainType& x, JacobianRangeType& grad) const
  {
    DomainType y = x;
    y /= x.two_norm();
    DomainType p = y;
    p[0] -= 1.; // lets have the origin on the sphere...

    JacobianRangeType Ju;
    u_.jacobian(p, Ju);

    for (int i = 0; i < dimDomain; ++i) {
      DomainType ds = y;
      ds *= -y[i];
      ds[i] += 1;

      RangeType Jds;
      Ju.mv(ds, Jds);
      for (int j = 0; j < dimRange; ++j)
        grad[j][i] = Jds[j];
    }
  }

private:
  U u_;
};


template <class GridImp, class U>
static ProblemInterface<FunctionSpace<double, double, GridImp::dimensionworld, 1>>*
createProblemFromFunction(const int problemType, const U& u = U())
{
  switch (problemType) {
    case 0:
      return new LaplaceProblem<GridImp, U>();
    case 1:
      return new RotatingAnisotropyProblem<GridImp, U>();
    case 2:
      return new SphereProblem<GridImp, U>();
    default:
      std::cerr << "Wrong value for problem type, bye bye!" << std::endl;
  }
  return 0;
}


template <class GridImp>
static ProblemInterface<FunctionSpace<double, double, GridImp::dimensionworld, 1>>* createProblem()
{
  const int dimWorld = GridImp::dimensionworld;

  int problemFunction                      = 0; // default value
  const std::string problemFunctionTable[] = {"sin", "cos", "alberta", "corner", "heter"};
  // read Function from parameter file
  problemFunction = Parameter::getEnum("femhowto.problemsolution", problemFunctionTable, problemFunction);

  // default value
  int problemType                      = 0;
  const std::string problemTypeTable[] = {"laplace", "rotatinganisotropy", "sphere"};
  // read problem type from parameter file
  problemType = Parameter::getEnum("femhowto.problemtype", problemTypeTable, problemType);

  // see problem.hh for available problems
  switch (problemFunction) {
    case 0:
      return createProblemFromFunction<GridImp>(problemType, SinProduct<dimWorld>());
    case 1:
      return createProblemFromFunction<GridImp>(problemType, CosProduct<dimWorld>());
    case 2:
      return createProblemFromFunction<GridImp>(problemType, AlbertaSolution<dimWorld>());
    case 3:
      return new LaplaceProblem<GridImp, CornerSolution<dimWorld>>();
    default:
      std::cerr << "Wrong problem value, bye, bye!" << std::endl;
      abort();
  }
  return 0;
}
}

#endif // #ifndef PROBLEM_HH
