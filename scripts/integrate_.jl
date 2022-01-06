using preCOPSE
using PyPlot
using OrdinaryDiffEq

pygui(true)

## define weathering as a function of pCO2 alone

function weathering_mac(pCO2)
    #temperature
    T = 288 + 4.5*log2(pCO2/285e-6)
    #runoff
    q = max(0.2*(1/year)*(1 + 0.03*(T - 288)), 0)
    #weathering rate per unit area
    w = mac(q, T, pCO2, 11.1, 288, 285e-6,Λ=0.005455)#,L=1.5587652489,Λ=0.0084)
    #scale to 30 % land fraction
    w*0.3*𝐒ₑ*year
end

wref = W/(0.3*𝐒ₑ*year)
function weathering_whak(pCO2)
    #temperature
    T = 288 + 4.5*log2(pCO2/285e-6)
    #runoff
    q = max(0.2*(1/year)*(1 + 0.03*(T - 288)), 0)
    #weathering rate per unit area
    w = whak(q, T, pCO2, 0.24506705893859612, 11.1, 288, 285e-6)#,L=1.5587652489,Λ=0.0084)
    #scale to 30 % land fraction
    w*0.3*𝐒ₑ*year
end

## parameter values

#param = initparams(h=2.326e20,k1=7.5e12,k7=7.5e12,V=7.5e12)
#param = initparams(h=2.326e20,k1=1.05e13,k7=1/05e13,V=7.5e12,k3=0)
param = initparams(h=2.326e20,k8=3.3e10,k3=0,k7=6.375e12,V=7.4992382839015205e12)#,k1=0)#,k1=0,k7=0,V=7.5e12)

## integrate
p = (weathering_mac, param)
u0 = [3.7e15, 3.193e18] #1e20
u0 = [3.7e15, 1e20] #1e20
u0 = [6e15, 1.25e20] #1e20
tspan = (0.0, 1e6)
prob = ODEProblem(ℱ!, u0, tspan, p)
sol1 = solve(prob, RadauIIA5())

## integrate
p = (weathering_whak, param)
u0 = [3.7e15, 3.193e18] #1e20
u0 = [3.7e15, 1e20] #1e20
u0 = [6e15, 1.25e20] #1e20
tspan = (0.0, 1e5)
prob = ODEProblem(ℱ!, u0, tspan, p)
sol2 = solve(prob, RadauIIA5())

figure()
plot(sol2.t,𝒻pCO2.(sol2[2,:],param.h)*1e6)

figure()
plot(sol1.t,𝒻pCO2.(sol1[2,:],param.h)*1e6)
plot(sol2.t,𝒻pCO2.(sol2[2,:],param.h)*1e6)


# figure()
# plot(sol.t,sol[1,:])