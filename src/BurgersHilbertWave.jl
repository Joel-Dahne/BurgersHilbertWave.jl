module BurgersHilbertWave

using ArbExtras
using Arblib
using NLsolve
using OffsetArrays
using OrderedCollections
using SpecialFunctions

include("fmpz.jl")
include("arb.jl")
include("series.jl")
include("special-functions/special-functions.jl")
include("special-functions/polylog.jl")
include("special-functions/clausenc.jl")
include("special-functions/clausens.jl")
include("types.jl")
include("evaluation.jl")

include("bounds.jl")
include("estimates.jl")
include("T0.jl")
include("proof.jl")

include("FractionalKdV/FractionalKdVAnsatz.jl")
include("FractionalKdV/determination.jl")
include("FractionalKdV/evaluation.jl")

include("BurgersHilbert/BHAnsatz.jl")
include("BurgersHilbert/determination.jl")
include("BurgersHilbert/evaluation.jl")
include("BurgersHilbert/n0.jl")
include("BurgersHilbert/delta0.jl")
include("BurgersHilbert/T0.jl")
include("BurgersHilbert/T01.jl")
include("BurgersHilbert/T02.jl")
include("BurgersHilbert/D0.jl")

end # module
