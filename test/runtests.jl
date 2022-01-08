using preCOPSE
using Test

const pCO2₀ = 285e-6

@testset "pre-industrial" begin
    @test integrate(1e6, mac₀) ≈ pCO2₀
    @test integrate(1e6, whak₀) ≈ pCO2₀
end

t = LinRange(1, 1e6, 100)
p = initparams(A₀=1.28e20)
@test all(integrate(t, mac₀, p) .> integrate(t, whak₀, p))