# Highest Cusped Waves for the Burgers-Hilbert equation

![Figure showing the approximation of u with error
bounds](figures/BH-u.png)

This repository contains the code for the computer assisted parts of
the proofs for the paper [Highest Cusped Waves for the Burgers-Hilbert
equation](https://doi.org/10.1007/s00205-023-01904-6). If you are
interested in seeing the results of the computations without running
anything yourself you can look at the html-file in the `proofs`
directory. You will likely have to download the file and open it
locally in your browser.

## Reproducing the proof

The html-file mentioned above is generated from a
[Pluto](https://github.com/fonsp/Pluto.jl) notebook. You can reproduce
the proof by running the notebook yourself as described below.

The proofs were run with Julia version 1.8.5 but should likely work
with later versions as well. This repository contains the same
`Manifest.toml` file as was used when running the proofs, this allows
us to install exactly the same versions of the packages. To do this
start by downloading this repository, enter the directory and start
Julia. Now run the code

``` julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

This will likely take some time the first time you run it since it has
to download all the packages. You can see if it seems to work by
running `Pkg.test()`. This should give you some output related to all
the installed packages and test summary hopefully indicating that all
tests passed.

To run the notebooks do

``` julia
using Pluto
Pluto.run()
```

which should open a Pluto window in your browser. Now you can open the
notebook inside the `proofs` directory through this and it should run
the proof.

It is also possible to run the proof directly from the terminal with

``` shell
cd proofs
julia -t auto bh.jl
```

# Notes about the implementation
The code uses [Arblib.jl](https://github.com/kalmarek/Arblib.jl),
which is an interface to [Arb](https://www.arblib.org/), for rigorous
numerics. Some of the methods used for computing bounds are from
[ArbExtras.jl](https://github.com/Joel-Dahne/ArbExtras.jl), most
notably `ArbExtras.enclose_maximum` which is used for computing most
of the bounds.

The code includes implementations of some special functions, these are
found in `src/special-functions/`. Most notably are the
implementations of the Clausen functions, these are found in
`src/special-functions/clausenc.jl` and
`src/special-functions/clausens.jl`. The other functions are mostly
wrappers for Arb implementations, sometimes with slight modifications.

The approximation u₀ is represented by the type `BHAnsatz`, the
default arguments give the approximation used in the paper, so it can
be computed with`u0 = BHAnsatz{Arb}()`. Most of the code for
evaluating `u0` and computing bounds of the required constants are in
`src/BurgersHilbert/`. For the construction it first computes an
approximation for the Fractional KdV equation, the code for this is
found in `src/FractionalKdV/`.

The methods for handling removable singularities are found in
`src/series.jl`.
