using preCOPSE
using PyPlot
using OrdinaryDiffEq

pygui(true)

## define weathering as a function of pCO2 alone

function weathering(pCO2)
    #temperature
    T = 288 + 4.5*log2(pCO2/285e-6)
    #runoff
    q = max(0.2*(1/year)*(1 + 0.03*(T - 288)), 0)
    #weathering rate per unit area
    w = mac(q, T, pCO2, 11.1, 288, 285e-6,Λ=0.005455)#,L=1.5587652489,Λ=0.0084)
    #scale to 30 % land fraction
    w*0.3*𝐒ₑ*year
end

## parameter values

#param = initparams(h=2.326e20,k1=7.5e12,k7=7.5e12,V=7.5e12)
#param = initparams(h=2.326e20,k1=1.05e13,k7=1/05e13,V=7.5e12,k3=0)
param = initparams(h=2.326e20,k8=3.3e10,k3=0,k7=4.5e12,V=7.4992382839015205e12)#,k1=0)#,k1=0,k7=0,V=7.5e12)

## integrate

p = (weathering, param)
u0 = [3.1e15, 3.193e18] #1e20
tspan = (0.0, 1e7)
prob = ODEProblem(ℱ!, u0, tspan, p)
sol = solve(prob, RadauIIA5())

figure()
plot(sol.t,𝒻pCO2.(sol[2,:],param.h)*1e6)

figure()
plot(sol.t,sol[1,:])