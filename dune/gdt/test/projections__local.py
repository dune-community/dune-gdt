import spaces as sp

dg = sp.DG(cache)
fv = sp.FV(cache, rdim=[1])
rt = sp.RT(cache)
spaces = dg.spaces + fv.spaces + rt.spaces
names = dg.names + fv.names + rt.names
spaces_with_names = zip(spaces, names)
