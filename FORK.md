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
  cutoff) checks every atom, as upstream does. In a non-orthorhombic cell the cells are
  boxes in its slanted coordinates, searched to the cutoff divided by sqrt(lambda), where
  lambda bounds the smallest eigenvalue of U^T U for the unslant matrix U (Gershgorin),
  so the Cartesian distance of any atom outside them is still beyond the cutoff:
  0.5 for a hexagonal cell. The distances there include the multiplication by U, in
  the same order as upstream.
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
  test checks only atoms in nearby cells (in any cell shape, as for the lattice step).
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

## Exact labelling in parallel, and with less memory

Exact labelling (`exact_percolation` in `percolation.f90`) runs on threads: each thread
labels its own block of z planes with union-find, the blocks are joined across the planes
between them, and every site is pointed at its root in parallel. The root of a set is its
smallest site, so a cluster's root is its first site in scan order, which is the order
clusters are numbered in, and the result does not depend on the number of threads. It
keeps one 4-byte parent per site, instead of the 8 bytes (labels and parents) of
`clusteranalysis_exact` and `span`, and no table sized by the number of clusters (the old
exact path stopped at 100,000 clusters in the limiting diameter's analysis).
`tests/test_percolation_mpi.f90` checks it against `clusteranalysis_exact` and `span` on
300 random lattices, on any number of threads and ranks.

## Wide pores and large structures in the lattice step

- **A widening search instead of every atom.** A cube whose nearest atom or surface could
  lie beyond the cutoff (a pore wider than about twice the cutoff) is retried with the
  atoms within twice the radius, then four times, and so on, until its answer is settled
  or the list holds every atom. Upstream, and the fork before, checked every atom at
  once. The answer is the same: the atoms added are beyond the cutoff, so the
  Lennard-Jones sum and the overlap test do not change, and they are farther than the
  nearest atom and surface found. The wider list is built once per cell-list cell.
- **Candidate atoms gathered from the cells.** Each cell-list cell's candidate atoms are
  gathered from the neighbouring cells and sorted, instead of marking them and scanning
  every atom in the structure. The cost now follows the atoms nearby, not the size of the
  structure (it was the number of cells times the number of atoms).

## Splitting the grid over MPI ranks

`make -C src mpi` builds `poreblazer-mpi.exe` with `mpif90` (`-DUSE_MPI`; module
`parallel`, which is a set of single-rank stubs in the serial build). Run it under mpirun
or srun, with exact labelling, which only it can split:

```
OMP_NUM_THREADS=4 mpirun -n 8 ./poreblazer-mpi.exe < input.dat
```

Each rank holds a slab of whole z planes of the grid; every rank reads the whole
structure, as the atoms are few beside the grid cubes. Each step:

- **Lattice:** each rank computes its own planes; no communication.
- **Percolation (helium, nitrogen, limiting diameter):** each rank labels its slab, then
  sends the roots of its bottom plane to the rank below (the halo, periodic). The pairs of
  roots that meet across a boundary, a few per cluster that crosses one, are shared by
  every rank and joined with union-find, so a cluster's global root is its first site in
  the whole grid. Clusters crossing ranks add up their sizes and occupied planes over the
  ranks; clusters inside one slab need nothing from the others.
- **Helium volume:** the Boltzmann factors are computed on every rank in parallel; the sum
  is passed from rank to rank in plane order, so it is added in upstream's order.
- **Surface area:** every rank draws every random number (the stream stays in step), and
  tests only the trial points in its own planes; the counts are summed.
- **Pore size distribution:** each sample's cube is found by the rank holding it, and
  every rank finds the largest sphere among its own cubes containing the point; the
  sample's answer is the largest over the ranks.
- **Limiting diameter:** each rank sorts its own radii; the bisection needs only the
  smallest and largest radius, how many lie below the probe, and the radii next to that
  position, which the ranks find together.
- **Output:** rank 0 writes the log and every file; the other ranks send it their planes
  of `nitrogen_network.grd` and their parts of `nitrogen_network.xyz`, in order. Only rank
  0 reads stdin, and it passes the three lines on.

On any number of ranks the output is identical to the serial build's: `tests/run_tests.py
--fork src/poreblazer-mpi.exe --mpi 1 2 3 4` checks the log (less the line naming the
ranks) and every file against `tests/reference_exact.json`. The log says `MPI ranks: N`
when there are more than one. Memory per rank is the fixed part plus the slab's share of
the grid. Global cube counts are 8-byte, so a grid may hold more than 2**31 cubes; one rank
may hold up to 2**31.

`POREBLAZER_TIMING=1` makes rank 0 write each step's wall time to stderr (not the log).

## Where the grid falls

Poreblazer starts its grid at the lowest atom. `POREBLAZER_GRID_OFFSET="fx fy fz"` moves the
structure against the grid by those fractions of a cubelet along the cell's sides (the log
says `Grid offset (cubelets): ...`), so the spread of a result over offsets shows how much
of it is where the grid happens to fall. A whole number of cubelets changes no result in a
cell that is a whole number of cubelets on every side (`tests/run_tests.py --grid-offset`).
Only the first side is fitted to the cubelet, as in 3.0.5: on the others, up to one
cubelet's width of the cell is left out of the grid, and an offset moves that strip.

Measured on the three Ambuild cells in `tests/cases`, 16 random offsets at each cubelet
size, against the mean of 4 offsets at 0.1 A (exact labelling, November 2026):

| Result | Spread over offsets (sd) at 0.2 A | Error at 0.2 A, one offset | Averaging offsets |
| --- | --- | --- | --- |
| Pore limiting diameter | 0.01-0.03 A | -0.02 to -0.09 A | does not help: the error is a bias, about -0.6 to -0.7 A per A of cubelet in the two open cells |
| Maximum pore diameter | 0.02-0.03 A | -0.05 to -0.08 A | does not help, as above |
| Geometric volume | up to 0.004% of the volume | up to 0.003% | removes it (mean within 0.0003%) |
| Helium volume | under 0.00002% | none | converged already |
| Surface area | none up to 0.35 A | none | does not help at 0.7 A, where it is 2% low: a bias |
| Percolation dimensions | none at 0.2-0.7 A | see below | |

In the dense 30 A cell the network percolates in 3 dimensions at 2 of the 4 offsets at
0.1 A, and in 1 dimension at the other two and at every one of the 64 runs from 0.2 to
0.7 A: a channel whose narrowest point is close to the probe's size is found or missed
by where the grid falls. A dimension count from one grid is uncertain for such cells.
