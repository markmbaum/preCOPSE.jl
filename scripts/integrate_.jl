using preCOPSE
using PyPlot
using OrdinaryDiffEq

pygui(true)

##

#the Earth's mean radius [m]
const 𝐑ₑ = 6.371e6

#the Earth's surface area [m^2]
const 𝐒ₑ = 4π*𝐑ₑ^2

#seconds in a year
const year = 31536000.0

## define weathering as a function of pCO2 alone

function weathering(pCO2)
    #temperature
    T = 278 + 4.5*log2(pCO2/285e-6)
    #runoff
    q = max(0.2*(1/year)*(1 + 0.03*(T - 288)), 0)
    #weathering rate per unit area
    w = mac(q, T, pCO2, 11.1, 288, 285e-6)
    #scale to 30 % land fraction
    w*0.3*𝐒ₑ*year
end

## parameter values

param = initparams()

## integrate

p = (weathering, param)
u0 = [3.1e15, 1e20]
tspan = (0.0, 1e7)
prob = ODEProblem(ℱ!, u0, tspan, p)
sol = solve(prob, RadauIIA5())