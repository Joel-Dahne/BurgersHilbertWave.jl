using Arblib, BurgersHilbertWave, OffsetArrays, SpecialFunctions, Test

setprecision(Arb, 128) do
    @testset "BurgersHilbertWave" begin
        @testset "Special functions" begin
            include("arb.jl")
            include("clausenc.jl")
            include("clausens.jl")
        end

        @testset "BurgersHilbert" begin
            include("BurgersHilbert/evaluation.jl")
            include("BurgersHilbert/T0.jl")
        end
    end
end
