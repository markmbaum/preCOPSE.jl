using preCOPSE
using Test

const pCO2₀ = 285e-6

@testset "pre-industrial" begin
    @test integrate(1e6, mac₀) ≈ pCO2₀
    @test integrate(1e6, whak₀) ≈ pCO2₀
end
