# hz-gsvd
The Hari–Zimmermann variant of the Falk–Langemeyer algorithm for the generalized SVD.

This software is a supplementary material for the paper
doi:[10.1016/j.parco.2015.06.004](https://doi.org/10.1016/j.parco.2015.06.004 "Blocking and parallelization of the Hari–Zimmermann variant of the Falk–Langemeyer algorithm for the generalized SVD").

Only the block-oriented algorithm variants are included, up to the shared-memory thread-parallel (OpenMP) level.
Distributed-memory (MPI) level is still under cleanup and will be added eventually.

## HDF5

For now, the HDF5 library cannot be built with the ``-fdefault-integer-8`` (GFortran) or ``-i8`` (Intel Fortran) flag.
Please use a pre-2.0.0 HDF5, e.g., from the tag [final_autotools_develop](https://github.com/HDFGroup/hdf5/releases/tag/final_autotools_develop) instead.
An example for configuring it on macOS (Intel):
```bash
./autogen.sh
CC=gcc-15 FC=gfortran-15 FCFLAGS="-fdefault-integer-8" CXX=g++-15 ./configure --prefix=/opt/hdf5 --enable-build-mode=production --enable-dependency-tracking --enable-fortran --enable-cxx --enable-test-express=0 --disable-doxygen-doc --enable-pic --disable-parallel --without-zlib --without-szlib
```

This work has been supported in part by Croatian Science Foundation under the project IP-2014-09-3670 ([MFBDA](https://web.math.pmf.unizg.hr/mfbda/)).
