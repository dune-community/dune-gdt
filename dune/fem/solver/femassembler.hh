#ifndef FEMASSEMBLER_AZSG0L97
#define FEMASSEMBLER_AZSG0L97


namespace Dune
{
namespace Fem
{
namespace Functional
{
namespace Solver
{

  template<class MatrixImp, class VectorImp>
  class FEMAssembler
  {
  public:
    typedef MatrixImp
      MatrixType;

    typedef VectorImp
      VectorType;


  public:
    template<class Operator>
    static void assembleMatrix(Operator & op, MatrixType & m)
    {
      typedef typename Operator::DiscreteFunctionSpaceType
        DFS;
      typedef typename DFS::BaseFunctionSet
        BFS;
      typedef typename DFS::Iterator
        ItType;
      typedef typename ItType::Entity
        Entity;

      DFS & space = op.space();
      // check that number of dofs in space is equal to matrix size
      // assert()
      ItType it = space.begin();
      for(; it!=space.end(); ++it)
      {
        LocalMatrix localMatrix = op.applyLocal( en1 )

      }


    }

  };

} // end of namespace Constraints
} // end of namespace Functional
} // end of namespace Fem
} // end of namespace Dune


#endif /* end of include guard: FEMASSEMBLER_AZSG0L97 */
