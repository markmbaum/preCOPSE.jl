using preCOPSE

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

