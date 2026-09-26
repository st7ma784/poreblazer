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
- The Makefile builds with `-fopenmp`. `OMP_NUM_THREADS` sets the number of threads.

Results do not depend on the number of threads, and match upstream's serial build.
Ambuild's profile of Poreblazer and the measurements behind these changes are in its
`docs/benchmarks.md`.
