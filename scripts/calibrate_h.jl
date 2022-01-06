using Roots
using preCOPSE: 𝒻pCO2

##

# present day atmosphere/ocean CO2 [mole]
A₀ = 3.193e18
# target pCO2 value [bar]
pCO2 = 285e-6
#initial guess for h
h₀ = 2e20

#need to scale the value of h to an appropriate magnitude for root finding
h = 1e20*find_zero(x->𝒻pCO2(A₀, x*1e20) - pCO2, h₀/1e20)

println("at h = $h, pCO2 = $pCO2")