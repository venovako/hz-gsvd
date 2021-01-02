# hz-gsvd
The Hari–Zimmermann variant of the Falk–Langemeyer algorithm for the generalized SVD.

This software is a supplementary material for the paper
doi:[10.1016/j.parco.2015.06.004](https://doi.org/10.1016/j.parco.2015.06.004 "Blocking and parallelization of the Hari–Zimmermann variant of the Falk–Langemeyer algorithm for the generalized SVD").

Only the block-oriented algorithm variants are included, up to the shared-memory thread-parallel (OpenMP) level.
Distributed-memory (MPI) level is still under cleanup and will be added eventually.

This work has been supported in part by Croatian Science Foundation under the project IP-2014-09-3670 ([MFBDA](https://web.math.pmf.unizg.hr/mfbda/)).
