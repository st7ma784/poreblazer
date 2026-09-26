# Ambuild's fork of Poreblazer

This fork of [richardjgowers/poreblazer](https://github.com/richardjgowers/poreblazer)
(v3.0.5, commit `a753c72`, unchanged upstream since 2018) is maintained for
[Ambuild](https://github.com/st7ma784/ambuild). Changes are on the `ambuild` branch.

## Changes

- **OpenMP enabled and made correct.** Upstream's `lattice_calculations` carried
  OpenMP directives that the Makefile never enabled. Enabled, the loop had data races:
  the Lennard-Jones temporaries were shared between threads, and the cubelet lists
  were filled using a counter read after other threads could have incremented it.
  The shared temporaries also made it slower on 4 threads than on 1. All temporaries
  are now private, and the cubelet lists are built after the parallel loop, in
  cubelet order.
- **`pore_distribution` runs in parallel.** The 10,000 sample sites are drawn first,
  in the original order, then sampled in parallel; the distribution is accumulated in
  sample order afterwards.
- **Cell list in `lattice_calculations`.** Atoms are binned into cells about 2 A wide,
  and each grid cubelet checks only the atoms in cells that can lie within the cutoff,
  in ascending atom order, so the Lennard-Jones sum, the overlap test and the nearest
  atom are exactly those found by checking every atom. A cubelet whose nearest atom or
  nearest surface could lie beyond the cutoff (a pore wider than about twice the
  cutoff) checks every atom, as upstream does. Orthorhombic cells only; other cells
  check every atom.
- **Parallel sort** of the cubelets by pore radius before the pore size distribution.
  The order of equal radii can differ from the serial sort, which changes no result.
- **Faster cluster relabelling** in the percolation analysis: a lookup table instead of
  a search of every label so far for every site (quadratic in the number of clusters,
  and overflowing after 100,000 clusters). The labels are unchanged.
- **Faster `nitrogen_network.grd` output**: one write statement per plane instead of
  one per cubelet, about 30% faster, same file.
- The Makefile builds with `-fopenmp`. `OMP_NUM_THREADS` sets the number of threads.

## Opt-in exact cluster labelling

The percolation analysis's cluster labelling (`clusteranalysis`) records only one
level of cluster merges, so it can split one connected cluster into several labels
(on random 40^3 lattices, 54 of 60 were labelled differently from their true
components, worst near the percolation threshold). This affects the percolating
networks, the pore limiting diameter and the pore size distribution, and makes
spanning non-monotonic in probe radius, which misleads the limiting-diameter
bisection. Ambuild's white paper and `benchmarks/percolation_study/` have the details.

The fork keeps upstream's labelling by default, so its results match upstream's. To
use exact union-find labelling instead, add `1` after the visualisation option (line 7
of `defaults.dat`):

```
2, 1
```

The log states which labelling ran (`Percolation labelling: exact` or `poreblazer`).
A `defaults.dat` without the second value behaves exactly as before.
