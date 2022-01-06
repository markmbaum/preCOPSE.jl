using preCOPSE
using OrdinaryDiffEq

#------------------------------------------------------------------------------
# immutable physical constants

#the Earth's mean radius [m]
const 𝐑ₑ = 6.371e6

#the Earth's surface area [m^2]
const 𝐒ₑ = 4π*𝐑ₑ^2

#surface gravity [m/s^2]
const 𝐠 = 9.8

#seconds in a year
const 𝐲𝐫 = 31536000.0

#stefan-boltzmann constant [m^-2 K^-4]
const 𝛔 = 5.67e-8 

# Runoff linearization 
function q(T,T₀,ϵ,Γ,Pref)
    q = max(Γ*Pref*(1 + ϵ*(T-T₀)),0)
end

#------------------------------------------------------------------------------
function mac(pCO2, 
    Tₑ=11.1, 
    T₀=288,
    pCO2₀=285e-6,
    n=0.316,
    Λ=1.4e-3,
    L=1.0,
    ϕ=0.1,
    ρ=12728.0,
    k₀=8.7e-6,
    𝐀=1e2,
    X=0.36,
    tₛ=1e5,
    m=0.27,
    μ=exp(2),
    β=0.2,
    ϵ=0.03,
    Γ=0.2,
    Pref=1/𝐲𝐫,
    λ=4.5)
T = 278 + λ*log2(pCO2/pCO2₀)
#println(T)
r = q(T,T₀,ϵ,Γ,Pref)
#println(r*𝐲𝐫)
#defined for convenience
α = L*ϕ*ρ*𝐀*X*μ
#equilibrium concentration
Ceq = 1e3*Λ*(pCO2^n) #conversion from mol/liter to mol/m3
#temperature dependence
a = exp((T - T₀)/Tₑ)
#pCO2 dependence
b = (pCO2/pCO2₀)^β
#denominator
d = 1/(k₀*a*b) + m*𝐀*tₛ + α/(r*𝐲𝐫*Ceq)
#weathering per unit area 
(α/d)/𝐲𝐫
end

k1    = 9e-15    # mol C yr⁻¹ Total organic carbon burial
k2    = 1.5e10   # mol P yr⁻¹ Ca associated phosphorus burial
k3    = 6e9      # mol P yr⁻¹ Fe associated phosphorus burial
k4    = 0.86     # -          Initial oxic fraction
k5    = 6.65e12  # mol C yr⁻¹ Silicate weathering
k6    = 1.335e13 # mol C yr⁻¹ Carbonate weathering
k7    = 7.75e12  # mol C yr⁻¹ Oxidative weathering
k8    = 5.7e10   # mol P yr⁻¹ Reactive phosphorus weathering
k9    = 1.25e12  # mol C yr⁻¹ Organic carbon degassing
k10   = 6.65e12  # mol C yr⁻¹ Carbonate carbon degassing
CPsea = 250      # mol:mol    C:P burial ratio
A₀    = 3.193e18 # mol        Present day atmosphere/ocean CO2
P₀    = 3.1e15   # mol        Present day ocean phosphate
O₀    = 3.7e19   # mol        Present day atmosphere/ocean O2
h     = 2.32746  # mol        Partitioning value for pCO2
W₀    = 7.5e12   # mol/yr     Past value of 

p     = [mac,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,CPsea,A₀,P₀,O₀,h,W₀]
u0    = [3.1e15;1e20]
tspan = (0.0,1000*𝐲𝐫)
prob  = ODEProblem(precopse!,u0,tspan,p)
sol   = solve(prob)