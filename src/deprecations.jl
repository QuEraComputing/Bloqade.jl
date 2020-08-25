@deprecate soft_misloss(args...) gibbs_loss(args...)
@deprecate lattice_atoms(n, ff, geometry) square_lattice(n, ff)
@deprecate emulate!(r, ts, hs, cache) emulate!(r, ts, hs; cache=cache)
