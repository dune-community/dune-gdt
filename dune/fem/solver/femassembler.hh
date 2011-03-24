#ifndef DUNE_FEM_FUNCTIONALS_SOLVER_FEMASSEMBLER_HH
#define DUNE_FEM_FUNCTIONALS_SOLVER_FEMASSEMBLER_HH

#include <dune/fem/common/localmatrix.hh>
#include <dune/fem/common/localvector.hh>

namespace Dune {

namespace Functionals {

namespace Solver {

template <class MatrixImp, class VectorImp>
class FEMAssembler
{
public:
  typedef MatrixImp MatrixType;

  typedef VectorImp VectorType;

  typedef typename VectorType::field_type FieldType;

  typedef Dune::Functionals::Common::LocalMatrix<FieldType> LocalMatrixType;

  typedef Dune::Functionals::Common::LocalVector<FieldType> LocalVectorType;


public:
  template <class Operator>
  static void assembleMatrix(const Operator& op, MatrixType& m)
  {
    typedef typename Operator::DiscreteFunctionSpaceType DFSType;
    typedef typename DFSType::BaseFunctionSetType BFSType;
    typedef typename DFSType::IteratorType ItType;
    typedef typename ItType::Entity EntityType;

    const DFSType& space = op.space();
    // check that number of dofs in space is equal to matrix size
    // assert()
    ItType it = space.begin();
    for (; it != space.end(); ++it) {
      const EntityType& en = *it;
      const BFSType& bfs = space.baseFunctionSet(en);
      LocalMatrixType localMatrix(bfs.numBaseFunctions(), bfs.numBaseFunctions());
      localMatrix = op.applyLocal(en);

      addToMatrix(space, localMatrix, en, m);
    }
  }

  template <class Functional>
  static void assembleVector(const Functional& functional, VectorType& vec)
  {
    typedef typename Functional::DiscreteFunctionSpaceType DFSType;
    typedef typename DFSType::BaseFunctionSetType BFSType;
    typedef typename DFSType::IteratorType ItType;
    typedef typename ItType::Entity EntityType;

    const DFSType& space = functional.space();

    // check that number of dofs in space is equal to vector size
    assert(space.size() == vec.size());
    ItType it = space.begin();
    for (; it != space.end(); ++it) {
      const EntityType& en = *it;
      const BFSType& bfs   = space.baseFunctionSet(en);

      LocalVectorType localVector = functional.applyLocal(en, bfs);

      addToVector(space, localVector, en, vec);
    }
  }


  /// \todo merge later with assembleMatrix
  /// \todo implement a PrecompiledConstraints class which wraps an existing
  ///       Constraints class for efficiency at the cost of one grid walk
  template <class ConstrainedDFS>
  static void applyMatrixConstraints(const ConstrainedDFS& cSpace, MatrixType& m)
  {
    typedef typename ConstrainedDFS::ConstraintsType ConstraintsType;
    typedef typename ConstrainedDFS::DiscreteFunctionSpaceType DFSType;
    typedef typename ConstraintsType::LocalConstraintsType LocalConstraintsType;
    typedef typename DFSType::IteratorType ItType;
    typedef typename ItType::Entity EntityType;

    const ConstraintsType& constraints = cSpace.constraints();

    // check that number of dofs in space is equal to matrix size
    // assert()

    ItType it  = cSpace.space().begin();
    ItType end = cSpace.space().end();

    for (; it != end; ++it) {
      const EntityType& en = *it;

      const LocalConstraintsType localConstraints = constraints.local(en);

      setLocalConstraintsInMatrix(cSpace, localConstraints, en, m);
    }
  }

  // @todo implementation
  template <class ConstrainedDFS>
  static void applyVectorConstraints(const ConstrainedDFS& cSpace, VectorType& vec)
  {
    typedef typename ConstrainedDFS::ConstraintsType ConstraintsType;
    typedef typename ConstrainedDFS::DiscreteFunctionSpaceType DFSType;
    typedef typename ConstraintsType::LocalConstraintsType LocalConstraintsType;
    typedef typename DFSType::IteratorType ItType;
    typedef typename ItType::Entity EntityType;

    const ConstraintsType& constraints = cSpace.constraints();

    ItType it  = cSpace.space().begin();
    ItType end = cSpace.space().end();

    for (; it != end; ++it) {
      const EntityType& en = *it;

      const LocalConstraintsType localConstraints = constraints.local(en);

      for (unsigned int i = 0; i < localConstraints.rowDofsSize(); ++i) {
        vec[localConstraints.rowDofs(i)] = 0;
      }
    }
  }

public:
  /// \todo move to matrixContainer factory
  template <class DFSType, class Entity>
  static void addToMatrix(const DFSType& space, const LocalMatrixType& localMatrix, const Entity& en, MatrixType& m)
  {
    for (unsigned int i = 0; i < localMatrix.N(); i++) {
      for (unsigned int j = 0; j < localMatrix.M(); j++) {
        const int globalI = space.mapToGlobal(en, i);
        const int globalJ = space.mapToGlobal(en, j);

        m[globalI][globalJ] += localMatrix.get(i, j);
      }
    }
  }

  template <class DFSType, class Entity>
  static void addToVector(const DFSType& space, const LocalVectorType& localVector, const Entity& entity,
                          VectorType& vector)
  {
    for (unsigned int ii = 0; ii < localVector.size(); ++ii) {
      const int globalI = space.mapToGlobal(entity, ii);

      vector[globalI] += localVector[ii];
    }
  }

  /// \todo move to Constraints class
  template <class ConstrainedDFS, class Entity>
  static void
  setLocalConstraintsInMatrix(const ConstrainedDFS& cSpace,
                              const typename ConstrainedDFS::ConstraintsType::LocalConstraintsType& localConstraints,
                              const Entity& en, MatrixType& m)
  {

    for (unsigned int i = 0; i < localConstraints.rowDofsSize(); i++) {
      for (unsigned int j = 0; j < localConstraints.columnDofsSize(); j++) {
        m[localConstraints.rowDofs(i)][localConstraints.columnDofs(j)] = localConstraints.localMatrix(i, j);
      }
    }
  }


}; // end of class FEMAssembler

} // end of namespace Solver

} // end of namespace Functionals

} // end of namespace Dune


#endif /* end of include guard: DUNE_FEM_FUNCTIONALS_SOLVER_FEMASSEMBLER_HH */
