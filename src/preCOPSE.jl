module preCOPSE

using UnPack
using OrdinaryDiffEq
using GEOCLIM: mac, whak
export mac, whak

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

#molar mass of CO2 [kg/mole]
const 𝛍 = 0.044

#stefan-boltzmann constant [m^-2 K^-4]
const 𝛔 = 5.67e-8

#------------------------------------------------------------------------------
# initializing parameter values for integration

export initparams

function initparams(;
    W0::Real=7.5e12,
    h::Real=2.32746,
    k1::Real=4.5e12,
    k2::Real=1.5e10,
    k3::Real=6e9,
    k7::Real=7.75e12,
    k8::Real=5.7e10,
    CPsea::Real=250,
    P0::Real=3.1e15,
    O::Real=1.76e18,
    O0::Real=3.7e19,
    V::Real=7.9e12)
    #construct a named tuple containing all the parameters
    (
        W0    = Float64(W0),
        h     = Float64(h),
        k1    = Float64(k1),
        k2    = Float64(k2),
        k3    = Float64(k3),
        k7    = Float64(k7),
        k8    = Float64(k8),
        CPsea = Float64(CPsea),
        P0    = Float64(P0),
        O     = Float64(O),
        O0    = Float64(O0),
        V     = Float64(V)
    )
end

#------------------------------------------------------------------------------
# the system of ODES

export 𝒻ϕ, 𝒻pCO2, 𝒻mocb, 𝒻fepb, ℱ!

#fraction of carbon in the atmososphere [-]
# A - total ocean-atmosphere CO2 [mole]
𝒻ϕ(A, h) = 0.78*A/(A + h)

#partial pressure of CO2 [bar]
# A - total ocean-atmososphere CO2 [mole]
𝒻pCO2(A, h) = (𝒻ϕ(A,h)*A*𝛍*𝐠/𝐒ₑ)/1e5

𝒻phosw(k₈, W, W₀) = k₈*(W/W₀ + 5/12)

𝒻mopb(mocb, CPsea) = mocb/CPsea

𝒻capb(k₂, mocb, mocb₀) = k₂*mocb/mocb₀

# Iron-sorbed phosphate burial 
#𝒻fepb(P, P₀, O, O₀, k₁, k₃) = k₃*(O/O₀*P₀/P)*𝒻mocb(k₃,P,P₀)
𝒻fepb(k₃, O, O₀, P, P₀, mocb, mocb₀) = k₃*(O/O₀)*(P₀/P)*(mocb/mocb₀)

# Marine organic carbon burial 
𝒻mocb(k₁, P, P₀) = k₁*(P/P₀)^2

𝒻oxidw(k₇, Wglacial=1.0) = k₇*Wglacial

function ℱ!(du, u, param, t)::Nothing
    #unpack phosphate and total carbon from input vector
    P, A = u
    #the weathering function is the first element of param
    𝒻W = param[1]
    #numeric parameters are in a named tuple in the second element
    @unpack h, k8, W0, k1, P0, CPsea, k2, k3, O, O0, k7, V = param[2]
    #atmospheric carbon concentration [bar]
    pCO2 = 𝒻pCO2(A, h)
    #weathering rate [mol/yr]
    W = 𝒻W(pCO2)
    #derived quantities
    phosw = 𝒻phosw(k8, W, W0)
    mocb = 𝒻mocb(k1, P, P0)
    mocb₀ = k1
    mopb = 𝒻mopb(mocb, CPsea)
    capb = 𝒻capb(k2, mocb, mocb₀)
    fepb = 𝒻fepb(k3, O, O0, P, P0, mocb, mocb₀)
    oxidw = 𝒻oxidw(k7)
    #evaluate dP/dt
    du[1] = phosw - mopb - capb - fepb
    #evaluate dA/dt
    du[2] = -mocb - W + oxidw + V
    #empty return value
    nothing
end

end
