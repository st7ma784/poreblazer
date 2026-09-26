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
- **Pore size distribution search by blocks.** Each sample's answer is the largest
  sphere (of a geometrically accessible cubelet) containing its point. Upstream found it
  by scanning every cubelet from the largest radius down. The fork groups cubelets into
  blocks, visits the blocks from the largest radius they hold down, stops when no block
  left can beat the best sphere found, and skips blocks farther from the point than their
  largest radius. It uses upstream's containment test, so it finds the same sphere.
- **Spanning test in one pass.** The percolation analysis finds the planes every
  candidate cluster occupies in one pass over the grid, instead of up to three passes
  per candidate. The per-cluster test and its order are unchanged.
- **Surface area in parallel, with a cell list.** Every trial draws its two random
  numbers before any test, so the numbers are drawn first, in upstream's order. The
  atoms then run in parallel, and their areas are summed in atom order. The overlap
  test checks only atoms in nearby cells (orthorhombic cells).
- **Vectorised distances in the lattice step.** Distances to the candidate atoms are
  computed in a loop of their own over contiguous arrays, with the arithmetic of
  `fundcell_snglMinImage`. `anint` is written in a form the compiler vectorises: it
  truncates |q|, adds 1 when the part dropped is at least 1/2, and restores the sign,
  which is exact for these arguments.
- **About half the memory** (34 instead of about 70 bytes per grid cube). The
  accessibility masks are 1 byte instead of 2. Each cube's indices are computed from its
  number instead of stored. The unused list of geometric cubes is gone. The helium volume
  reads the helium mask instead of a list (percolation leaves the same cubes in both),
  and the nitrogen list is sized to its count. The sorted copy for the pore size
  distribution holds only the radii: the limiting diameter finds the cubes above a
  radius on the grid. Storage only, so results are unchanged.
- The Makefile builds with `-fopenmp`, and `OMP_NUM_THREADS` sets the number of threads.
  It also uses `-ffp-contract=off` (no fused multiply-adds, so results do not depend on
  the CPU) and `-fvect-cost-model=dynamic` (vectorises loops of unknown length; the
  operations are elementwise, so results are unchanged).

## Testing and the published image

`.github/workflows/ci.yml` runs on every push and pull request to `ambuild`:

- `tests/test_percolation.f90` checks the exact labelling against an independent flood
  fill on 200 random periodic lattices, and that Poreblazer's labelling still splits the
  8-site counterexample as upstream does.
- `tests/run_tests.py` builds upstream 3.0.5 (`a753c72`) beside the fork and runs both on
  upstream's example frameworks (HKUST-1, IRMOF-1, MIL-47(V), and the hexagonal MOF-180)
  and on three cells built by Ambuild (20 Å; 30 Å with 576 atoms; a near-empty 40 Å cell
  whose pores are wider than the cutoff). With the default labelling the fork's log (less
  the line naming the labelling) and every file it writes, `nitrogen_network.grd`
  included, must match upstream's at 1 and 4 threads. With exact labelling the output must
  match `tests/reference_exact.json` at both thread counts; for the three Ambuild cells
  that reference also matches an independent build of the exact labelling.

Once the tests pass on a push to `ambuild`, the workflow publishes the `Dockerfile` as
`ghcr.io/st7ma784/poreblazer`, tagged `sha-<commit>` and `ambuild`. The executable is
`/opt/poreblazer/poreblazer.exe` (it needs `libgfortran5` and `libgomp1`). Ambuild
copies it from a pinned `sha-` tag.

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
