import spaces as sp

cg = sp.CG(cache)
dg = sp.DG(cache)
fv = sp.FV(cache)
rt = sp.RT(cache)
spaces = cg.spaces + dg.spaces + fv.spaces + rt.spaces
names = cg.names + dg.names + fv.names + rt.names
spaces_with_names = zip(spaces, names)
