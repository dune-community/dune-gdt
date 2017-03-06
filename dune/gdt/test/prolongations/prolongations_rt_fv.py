import spaces as sp
from grids import LevelGrids

fv = sp.FV(cache, base=LevelGrids, rdim=[1])
rt = sp.RT(cache, base=LevelGrids)
spaces = fv.spaces + rt.spaces
names = fv.names + rt.names
spaces_with_names = zip(spaces, names)
