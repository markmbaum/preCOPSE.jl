using Roots
using Optim
using preCOPSE
import GEOCLIM: mac, whak
using PyPlot

##

# present day atmosphere/ocean CO2 [mole]
const A₀ = 3.193e18
# present day ocean phosphate [mole]
const P₀ = 6e15
# target pCO2 value [bar]
const pCO2₀ = 285e-6
#target weathering rate [mole/year]
const W₀ = 7.5e12
# weathering constants
const Tₑ = 11.1
const T₀ = 288

## simplified functions for temperature and runoff against pCO2

𝒻T(pCO2) = 288 + 4.5*log2(pCO2/285e-6)

𝒻q(pCO2) = max(0.2*(1/year)*(1 + 0.03*(𝒻T(pCO2) - 288)), 0)

## calibrate h for correc pCO2 concentration

#need to scale the value of h to an appropriate magnitude for root finding
h = 1e20*find_zero(x->𝒻pCO2(A₀, x*1e20) - pCO2₀, 2)
println("h = $h → pCO2 = $(𝒻pCO2(A₀,h))")

## calibrate MAC for balanced pre-industrial weathering by tuning Λ

function macₚ(pCO2; kw...)
    mac(𝒻q(pCO2), 𝒻T(pCO2), pCO2, Tₑ, T₀, pCO2₀; kw...)*0.3*𝐒ₑ*year
end

Λ = find_zero(x->macₚ(pCO2₀, Λ=x) - W₀, 0.005)
println("Λ = $Λ → macₚ = $(macₚ(pCO2₀, Λ=Λ))")

## calibrate WHAK for balanced pre-industrial weathering by tuning k

function whakₚ(pCO2, k)
    whak(𝒻q(pCO2), 𝒻T(pCO2), pCO2, k, Tₑ, T₀, pCO2₀)*0.3*𝐒ₑ*year
end

k = find_zero(x->whakₚ(pCO2₀, x) - W₀, 0.25)
println("k = $k → whacₚ = $(whakₚ(pCO2₀, k))")

## optimize for MAC steady state

function macloss(k₇, k₈, Λ)
    params = initparams(k₇=k₇, k₈=k₈)
    sum(precopse(P₀, A₀, pCO2->macₚ(pCO2; Λ=Λ), params).^2)
end

res = optimize(
    x -> macloss(x[1]*1e12, x[2]*1e10, Λ),
    [6.3, 3.3]
)

k₇, k₈ = res.minimizer
k₇ *= 1e12
k₈ *= 1e10
println("MAC steady case: k₇ = $k₇, k₈ = $k₈")
params = initparams(k₇=k₇, k₈=k₈)
t, pCO2 = integrate(1e6, pCO2->macₚ(pCO2, Λ=Λ), params)
plot(t, pCO2)

## optimize for WHAK steady state

function whakloss(k₇, k₈, k)
    params = initparams(k₇=k₇, k₈=k₈)
    sum(precopse(P₀, A₀, pCO2->whakₚ(pCO2, k), params).^2)
end

res = optimize(
    x -> whakloss(x[1]*1e12, x[2]*1e10, k),
    [6.3, 3.3]
)

k₇, k₈ = res.minimizer
k₇ *= 1e12
k₈ *= 1e10
println("WHAK steady case: k₇ = $k₇, k₈ = $k₈")
params = initparams(k₇=k₇, k₈=k₈)
t, pCO2 = integrate(1e6, pCO2->whakₚ(pCO2, k), params)
plot(t, pCO2)
