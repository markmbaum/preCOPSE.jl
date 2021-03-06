module preCOPSE

using UnPack
using OrdinaryDiffEq
using GEOCLIM: mac, whak
export mac, whak

#------------------------------------------------------------------------------
# immutable physical constants

export ğâ, ğâ, ğ , year

#the Earth's mean radius [m]
const ğâ = 6.371e6

#the Earth's surface area [m^2]
const ğâ = 4Ï*ğâ^2

#surface gravity [m/s^2]
const ğ  = 9.8

#seconds in a year
const year = 31536000.0

#molar mass of CO2 [kg/mole]
const ğ = 0.044

#------------------------------------------------------------------------------
# initializing parameter values for integration

export initparams

#==============================================================================
Handy function for creating named tuple of parameters with default values
which can be overridden by keywords. Parameters will all have the same type,
as specified by the only non-keyword arg.

Parameters:
 Pâ - initial ocean phosphate reservoir [mole]
 Aâ - initial ocean-atmosphere CO2 reservoir [mole]
 Wâ - reference weathering rate [mole/yr]
 h - parameter fit for partitioning CO2 between ocean & atmosphere [mole]
 kâ - total organic carbon burial rate [mole/yr]
 kâ - Ca associated phosphorus burial [mole/yr]
 kâ - Fe associated phosphorus burial [mole/yr]
 kâ - oxidative weathering [mole/yr]
 kâ - reactive phosphorus weathering [mole/yr]
 CPsea - C:P burial ratio
 O - ocean-atmosphere oxygen reservoir [mole/yr]
 Oâ - reference (present-day) ocean-atmosphere oxygen reservoir [mole/yr]
 V - volcanic CO2 outgassing [mole/yr]
==============================================================================#
function initparams(ğ¯::Type=Float64;
                    Pâ::Real=6e15,
                    Aâ::Real=3.193e18,
                    Wâ::Real=7.5e12,
                    h::Real=2.3269250670587494e20,
                    kâ::Real=4.5e12,
                    kâ::Real=1.5e10,
                    kâ::Real=6e9,
                    kâ::Real=4.5e12,
                    kâ::Real=2.3495580286e10,
                    CPsea::Real=250,
                    O::Real=1.76e18,
                    Oâ::Real=3.7e19,
                    V::Real=7.5e12)
    (
        Pâ    = convert(ğ¯, Pâ),
        Aâ    = convert(ğ¯, Aâ),
        Wâ    = convert(ğ¯, Wâ),
        h     = convert(ğ¯, h),
        kâ    = convert(ğ¯, kâ),
        kâ    = convert(ğ¯, kâ),
        kâ    = convert(ğ¯, kâ),
        kâ    = convert(ğ¯, kâ),
        kâ    = convert(ğ¯, kâ),
        CPsea = convert(ğ¯, CPsea),
        O     = convert(ğ¯, O),
        Oâ    = convert(ğ¯, Oâ),
        V     = convert(ğ¯, V)
    )
end

#------------------------------------------------------------------------------
# the system of ODES

export ğ»Ï, ğ»pCO2, ğ»mocb, ğ»fepb, precopse, precopse!

#fraction of carbon in the atmososphere [-]
# A - total ocean-atmosphere CO2 [mole]
ğ»Ï(A, h) = 0.78*A/(A + h)

#partial pressure of CO2 [bar]
# A - total ocean-atmososphere CO2 [mole]
ğ»pCO2(A, h) = (ğ»Ï(A,h)*A*ğ*ğ /ğâ)/1e5

ğ»phosw(kâ, W, Wâ) = kâ*(W/Wâ + 5/12)

ğ»mopb(mocb, CPsea) = mocb/CPsea

ğ»capb(kâ, mocb, mocbâ) = kâ*mocb/mocbâ

# Iron-sorbed phosphate burial 
#ğ»fepb(P, Pâ, O, Oâ, kâ, kâ) = kâ*(O/Oâ*Pâ/P)*ğ»mocb(kâ,P,Pâ)
ğ»fepb(kâ, O, Oâ, P, Pâ, mocb, mocbâ) = kâ*(O/Oâ)*(Pâ/P)*(mocb/mocbâ)

# Marine organic carbon burial 
ğ»mocb(kâ, P, Pâ) = kâ*(P/Pâ)^2

ğ»oxidw(kâ, Wglacial=1.0) = kâ*Wglacial

function precopse(P, A, ğ»W::F, params) where {F}
    #unpack numeric parameters
    @unpack Pâ, Wâ, h, kâ, kâ, kâ, kâ, kâ, CPsea, O, Oâ, V = params
    #atmospheric carbon concentration [bar]
    pCO2 = ğ»pCO2(A, h)
    #weathering rate [mol/yr]
    W = ğ»W(pCO2)
    #derived quantities
    phosw = ğ»phosw(kâ, W, Wâ)
    mocb = ğ»mocb(kâ, P, Pâ)
    mocbâ = kâ
    mopb = ğ»mopb(mocb, CPsea)
    capb = ğ»capb(kâ, mocb, mocbâ)
    fepb = ğ»fepb(kâ, O, Oâ, P, Pâ, mocb, mocbâ)
    oxidw = ğ»oxidw(kâ)
    #evaluate dP/dt
    dP = phosw - mopb - capb - fepb
    #evaluate dA/dt
    dA = -mocb - W + oxidw + V
    return dP, dA
end

function precopse!(d, P, A, ğ»W::F, params)::Nothing where {F}
    d[1], d[2] = precopse(P, A, ğ»W, params)
    nothing
end

function precopse!(du, u, p, t)::Nothing
    @inbounds precopse!(du, u[1], u[2], p[1], p[2])
    nothing
end

#------------------------------------------------------------------------------
# integration of the model

export integrate

function integrate(tspan::Tuple,
                   ğ»W::F,
                   params::NamedTuple=initparams();
                   kw...) where {F}
    #initial conditions
    uâ = Float64[params[:Pâ], params[:Aâ]]
    #bundle function with numeric parameters
    p = (ğ»W, params)
    #problem definition
    prob = ODEProblem(precopse!, uâ, map(Float64, tspan), p)
    #run the solver
    solve(prob, Rodas4P(); kw...)
end

function integrate(t::Real,
                   ğ»W::F,
                   params::NamedTuple=initparams();
                   kw...) where {F}
    #run the solver
    sol = integrate((0.0, Float64(t)), ğ»W, params; kw...)
    #return only the CO2 concentration and time
    ğ»pCO2(sol[2,end], params[:h])
end

function integrate(t::AbstractVector,
                   ğ»W::F,
                   params::NamedTuple=initparams();
                   kw...) where {F}
    #run the solver
    sol = integrate((0.0, maximum(t)), ğ»W, params; kw...)
    #dense output for carbon reservoir
    A = sol(t, idxs=2)
    #return only the CO2 concentration and time
    ğ»pCO2.(A, params[:h])
end

#------------------------------------------------------------------------------
# versions of MAC and WHAK that are at rest with default parameters

export ğ»T, ğ»q, macâ, whakâ

const Tâ = 11.1
const Tâ = 288
const pCO2â = 285e-6

ğ»T(pCO2) = 288 + 4.5*log2(pCO2/285e-6)

ğ»q(pCO2) = max(0.2*(1/year)*(1 + 0.03*(ğ»T(pCO2) - 288)), 0)

function macâ(pCO2) 
    0.3*ğâ*year*mac(ğ»q(pCO2), ğ»T(pCO2), pCO2, Tâ, Tâ, pCO2â, Î=0.005455935160109789)
end

function whakâ(pCO2)
    0.3*ğâ*year*whak(ğ»q(pCO2), ğ»T(pCO2), pCO2, 0.24506705893859618, Tâ, Tâ, pCO2â)
end

end