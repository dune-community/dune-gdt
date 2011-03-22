#ifndef DUNE_FEM_FUNCTIONALS_CONTAINER_FACTORY_HH
#define DUNE_FEM_FUNCTIONALS_CONTAINER_FACTORY_HH

#include <dune/common/exceptions.hh>

namespace Dune {
namespace Functionals {
namespace Container {

/** @brief interface of static factory class for matrix and vector classes
 */
template <class ContainerImp>
class Factory
{
public:
  //! return type for create() method
  typedef NotImplemented AutoPtrType;
  //! wrapped container type
  typedef NotImplemented ContainerType;

private:
  // class NotImplemented
  //{
  //};

public:
  /** @brief creates a new matrix/vector object and returns an auto_ptr
   * pointing to the allocated object
   *
   * - Matrices have size @f$ H \times H @f$ where @f$H@f$ is the number
   * of degrees of freedom in the discrete function space @f$ { \cal X }_H @f$.
   * The matrices' sparsity pattern is determined by the discrete function
   * space's basefunction overlap.
   * - Vectors have @f$ H @f$ components
   *
   * @param dfs the discrete function space @f$ { \cal X }_H @f$.
   */
  template <class DiscFuncSpace>
  static AutoPtrType create(DiscFuncSpace& dfs)
  {
    dune_static_assert(false, "Not Implemented: Factory!");
  }

  /** @brief creates a new matrix/vector object and returns a pointer to the
   * allocated object
   *
   * - Matrices have size @f$ H \times H @f$ where @f$H@f$ is the number
   * of degrees of freedom in the discrete function space @f$ { \cal X }_H @f$.
   * The matrices' sparsity pattern is determined by the discrete function
   * space's basefunction overlap.
   * - Vectors have @f$ H @f$ components
   *
   * @param dfs the discrete function space @f$ { \cal X }_H @f$.
   */
  template <class DiscFuncSpace>
  static ContainerType* createPtr(DiscFuncSpace& dfs)
  {
    dune_static_assert(false, "Not Implemented: Factory!");
  }
}; // end of class Factory

/** @brief interface of static factory class for matrix classes
 */
template <class ContainerImp>
class MatrixFactory : public Factory<ContainerImp>
{
};

// specialization for BCRSMatrix
template <class T>
class MatrixFactory<Dune::BCRSMatrix<T>>
{
public:
  typedef Dune::BCRSMatrix<T> ContainerType;
  typedef std::auto_ptr<ContainerType> AutoPtrType;

public:
  /** @brief creates a new BCRSMatrix object and returns an auto_ptr pointing
   * to the allocated object
   *
   * Matrices have size @f$ H \times H @f$ where @f$H@f$ is the number
   * of degrees of freedom in the discrete function space @f$ { \cal X }_H @f$.
   * The matrices' sparsity pattern is determined by the discrete function
   * space's basefunction overlap.
   *
   * @param dfs the discrete function space @f$ { \cal X }_H @f$.
   */
  template <class DiscFuncSpace>
  static AutoPtrType create(DiscFuncSpace& dfs)
  {
    return AutoPtrType(createPtr(dfs));
  }

  /** @brief creates a new BCRSMatrix object and returns a pointer to the
   * allocated object
   *
   * Matrices have size @f$ H \times H @f$ where @f$H@f$ is the number
   * of degrees of freedom in the discrete function space @f$ { \cal X }_H @f$.
   * The matrices' sparsity pattern is determined by the discrete function
   * space's basefunction overlap.
   *
   * @param dfs the discrete function space @f$ { \cal X }_H @f$.
   */
  template <class DiscFuncSpace>
  static ContainerType* createPtr(DiscFuncSpace& dfs)
  {
    typedef typename DiscFuncSpace::BaseFunctionSet BFS;
    typedef typename DiscFuncSpace::IteratorType ItType;
    typedef typename ItType::Entity Entity;

    const unsigned int numDofs = dfs.size();
    typedef std::vector<std::set<unsigned int>> PatternType;
    PatternType sPattern(numDofs);

    ContainerType* matrix = new ContainerType(numDofs, numDofs, ContainerType::random);

    // compute sparsity pattern
    // \todo precompile this in linear subspace
    // \todo use constraints for sparsity pattern
    ItType it = dfs.begin();
    for (; it != dfs.end(); it++) {
      const Entity& en = *it;
      const BFS& bfs   = dfs.baseFunctionSet(en);

      for (unsigned int i = 0; i < bfs.numBaseFunctions(); i++) {
        unsigned int ii = dfs.mapToGlobal(en, i);
        for (unsigned int j = 0; j < bfs.numBaseFunctions(); j++) {
          unsigned int jj = dfs.mapToGlobal(en, j);
          sPattern[ii].insert(jj);
        }
      }
    }

    for (unsigned int i = 0; i < sPattern.size(); i++) {
      matrix->setRowSize(i, sPattern[i].size());
    }
    matrix->endrowsizes();

    for (unsigned int i = 0; i < sPattern.size(); ++i) {
      typedef std::set<unsigned int>::const_iterator SetIterator;
      SetIterator sit = sPattern[i].begin();
      for (; sit != sPattern[i].end(); sit++) {
        matrix->addindex(i, *sit);
      }
    }
    matrix->endindices();

    return matrix;
  }

}; // end of MatrixFactory<BCRSMatrix<T> >


/** @brief interface of static factory class for vector classes
 */
template <class ContainerImp>
class VectorFactory : public Factory<ContainerImp>
{
};

/** @brief static factory class for vector classes of type Dune::BlockVector
 */
template <class T>
class VectorFactory<Dune::BlockVector<T>>
{
public:
  //! \copydoc Factory::ContainerType
  typedef Dune::BlockVector<T> ContainerType;
  //! \copydoc Factory::AutoPtrType;
  typedef std::auto_ptr<ContainerType> AutoPtrType;

public:
  /** @brief creates a new vector object and returns an auto_ptr pointing to
   * the allocated object
   *
   * The vector has @f$ H @f$ components which is the number of degrees of
   * freedom of the given discrete function space @f$ {\cal X}_H @f$.
   *
   * @param dfs the discrete function space @f$ { \cal X }_H @f$.
   */
  template <class DFSType>
  static ContainerType* createPtr(DFSType& dfs)
  {
    const unsigned int numDofs = dfs.size();
    ContainerType* bv          = new ContainerType(numDofs);
    return bv;
  }

  template <class DFSType>
  static AutoPtrType create(DFSType& dfs)
  {
    return AutoPtrType(createPtr(dfs));
  }
}; // end of VectorFactory<BlockVector<T> >

} // end of namespace Container

} // end of namespace Functionals

} // end of namespace Dune


#endif /* end of include guard: DUNE_FEM_FUNCTIONALS_CONTAINER_FACTORY_HH */
