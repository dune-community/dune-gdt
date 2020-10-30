from dune.xt.grid import Dim, Cube, make_cube_grid

grid = make_cube_grid(Dim(1), [0], [1], [2])
d = grid.dimension

from dune.xt.functions import ConstantFunction, GridFunction
from dune.xt.la import Istl
from dune.gdt import ContinuousLagrangeSpace, MatrixOperator, LocalElementIntegralBilinearForm, LocalLaplaceIntegrand

space = ContinuousLagrangeSpace(grid, 1)
op = MatrixOperator(grid, space, space, Istl())

func = ConstantFunction(Dim(1), Dim(1), [1])
op += LocalElementIntegralBilinearForm(LocalLaplaceIntegrand(GridFunction(grid, func, (Dim(d), Dim(d)))))
del func
op.assemble()
