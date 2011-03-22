#ifndef FACTORY_FNGWXZ23
#define FACTORY_FNGWXZ23


namespace Dune {
namespace Functionals {
namespace Container {

template <class ContainerImp>
class Factory
{
public:
  typedef NotImplemented AutoPtrType;
  typedef NotImplemented ContainerType;

  class NotImplemented
  {
  };

public:
  template <class DiscFuncSpace>
  static AutoPtrType create(DiscFuncSpace& dfs)
  {
    dune_static_assert("Not Implemented: Factory!");
  }

  template <class DiscFuncSpace>
  static ContainerType* createPtr(DiscFuncSpace& dfs)
  {
    dune_static_assert("Not Implemented: Factory!");
  }
}; // end of class Factory

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
  template <class DiscFuncSpace>
  static AutoPtrType create(DiscFuncSpace& dfs)
  {
    return AutoPtrType(createPtr(dfs));
  }

  template <class DiscFuncSpace>
  static ContainerType* createPtr(DiscFuncSpace& dfs)
  {
    typedef DiscFuncSpace::BaseFunctionSet BFS;
    typedef DiscFuncSpace::IteratorType ItType;
    typedef ItType::Entity Entity;

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
        ii = dfs.mapToGlobal(en, i);
        for (unsigned int j = 0; j < bfs.numBaseFunctions(); j++) {
          jj = dfs.mapToGlobal(en, j);
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


template <class ContainerImp>
class VectorFactory : public Factory<ContainerImp>
{
};

// specialization for BlockVector<T>
template <class T>
class VectorFactory<Dune::BlockVector<T>>
{
public:
  typedef Dune::BlockVector<T> ContainerType;
  typedef std::auto_ptr<ContainerType> AutoPtrType;

public:
  template <class DFSType>
  static ContainerType* createPtr(DFSType& dfs)
  {
    const unsigned int numDofs = dfs.size();
    ContainerType* bv          = new ContainerType(numDofs);
    return bv;
  }

  static AutoPtrType create(DFSType& dfs)
  {
    return AutoPtrType(createPtr(dfs));
  }
}; // end of VectorFactory<BlockVector<T> >

} // end of namespace Container
} // end of namespace Functionals
} // end of namespace Dune


#endif /* end of include guard: FACTORY_FNGWXZ23 */
