#ifndef FEMASSEMBLER_AZSG0L97
#define FEMASSEMBLER_AZSG0L97


namespace Dune {
namespace Fem {
namespace Functional {
namespace Solver {

template <class MatrixImp, class VectorImp>
class FEMAssembler
{
public:
  typedef MatrixImp MatrixType;

  typedef VectorImp VectorType;

  typedef VectorType::FieldType FieldType;

  typedef LocalMatrix<FieldType> LocalMatrixType;


public:
  template <class Operator>
  static void assembleMatrix(const Operator& op, MatrixType& m)
  {
    typedef typename Operator::DiscreteFunctionSpaceType DFS;
    typedef typename DFS::BaseFunctionSet BFS;
    typedef typename DFS::Iterator ItType;
    typedef typename ItType::Entity Entity;

    const DFS& space = op.space();
    // check that number of dofs in space is equal to matrix size
    // assert()
    const ItType it = space.begin();
    for (; it != space.end(); ++it) {
      const Entity& en = *it;
      const BFS& bfs = space.baseFunctionSet(en);
      LocalMatrixType localMatrix(bfs.numBaseFunctions(), bfs.numBaseFunctions);
      localMatrix = op.applyLocal(en);

      addToMatrix(space, localMatrix, en, m);
    }
  }

  /// \todo merge later with assembleMatrix
  /// \todo implement a PrecompiledConstraints class which wraps an existing
  ///       Constraints class for efficiency at the cost of one grid walk
  template <class ConstrainedDFS>
  static void applyMatrixConstraints(const ConstrainedDFS& cSpace, MatrixType& m)
  {
    typedef typename ConstrainedDFS::Constraints Constraints;
    typedef typename Constraints::LocalConstraints LocalConstraints;
    typedef typename ConstrainedDFS::BaseFunctionSet BFS;
    typedef typename DFS::Iterator ItType;
    typedef typename ItType::Entity Entity;

    const Constraints& constraints = cSpace.constraints();

    // check that number of dofs in space is equal to matrix size
    // assert()

    ItType it = cSpace.begin();
    for (; it != cSpace.end(); ++it) {
      const LocalConstraints localConstraints = constraints.local(en);

      setLocalConstraintsInMatrix(space, localConstraints, en, m);
    }
  }

private:
  /// \todo move to matrixContainer factory
  template <class DFSType, class Entity>
  static void addToMatrix(const DFSType& space, const LocalMatrix& localMatrix, const Entity& en, MatrixType& m)
  {
    for (unsigned int i = 0; i < localMatrix.getN(); i++) {
      for (unsigned int j = 0; j < localMatrix.getM(); j++) {
        const int globalI = space.mapToGlobal(en, i);
        const int globalJ = space.mapToGlobal(en, j);

        m[globalI][globalJ] = localMatrix[i][j];
      }
    }
  }

  /// \todo move to Constraints class
  template <class ConstrainedDFS, class Entity>
  static void setLocalConstraintsInMatrix(const ConstrainedDFS& cSpace, const LocalConstraints& localConstraints,
                                          const Entity& en, MatrixType& m)
  {

    for (unsigned int i = 0; i < localConstraints.rowDofsSize(); i++) {
      for (unsigned int j = 0; j < localConstraints.columnDofsSize(); j++) {
        m[localConstraints.rowDofs(i)][localConstraints.columnDofs(j)] = localConstraints.localMatrix(i, j);
      }
    }
  }
};

} // end of namespace Constraints
} // end of namespace Functional
} // end of namespace Fem
} // end of namespace Dune


#endif /* end of include guard: FEMASSEMBLER_AZSG0L97 */
