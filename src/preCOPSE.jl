module preCOPSE

using UnPack
using OrdinaryDiffEq
using GEOCLIM: mac, whak
export mac, whak

#------------------------------------------------------------------------------
# immutable physical constants

export 𝐑ₑ, 𝐒ₑ, 𝐠, year

#the Earth's mean radius [m]
const 𝐑ₑ = 6.371e6

#the Earth's surface area [m^2]
const 𝐒ₑ = 4π*𝐑ₑ^2

#surface gravity [m/s^2]
const 𝐠 = 9.8

#seconds in a year
const year = 31536000.0

#molar mass of CO2 [kg/mole]
const 𝛍 = 0.044

#------------------------------------------------------------------------------
# initializing parameter values for integration

export initparams

#===
Handy function for creating named tuple of parameters with default values
which can be overridden by keywords. Parameters will all have the same type,
as specified by the only non-keyword arg.
===#
function initparams(𝒯::Type=Float64;
                    P₀::Real=6e15,
                    A₀::Real=3.193e18,
                    W₀::Real=7.5e12,
                    h::Real=2.326925e20,
                    k₁::Real=4.5e12,
                    k₂::Real=1.5e10,
                    k₃::Real=6e9,
                    k₇::Real=4.5e12,
                    k₈::Real=2.349558e10,
                    CPsea::Real=250,
                    O::Real=1.76e18,
                    O₀::Real=3.7e19,
                    V::Real=7.5e12)
    (
        P₀    = convert(𝒯, P₀),
        A₀    = convert(𝒯, A₀),
        W₀    = convert(𝒯, W₀),
        h     = convert(𝒯, h),
        k₁    = convert(𝒯, k₁),
        k₂    = convert(𝒯, k₂),
        k₃    = convert(𝒯, k₃),
        k₇    = convert(𝒯, k₇),
        k₈    = convert(𝒯, k₈),
        CPsea = convert(𝒯, CPsea),
        O     = convert(𝒯, O),
        O₀    = convert(𝒯, O₀),
        V     = convert(𝒯, V)
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

function precopse(P, A, 𝒻W::F, params) where {F}
    #unpack numeric parameters
    @unpack P₀, W₀, h, k₁, k₂, k₃, k₇, k₈, CPsea, O, O₀, V = params
    #atmospheric carbon concentration [bar]
    pCO2 = 𝒻pCO2(A, h)
    #weathering rate [mol/yr]
    W = 𝒻W(pCO2)
    #derived quantities
    phosw = 𝒻phosw(k₈, W, W₀)
    mocb = 𝒻mocb(k₁, P, P₀)
    mocb₀ = k₁
    mopb = 𝒻mopb(mocb, CPsea)
    capb = 𝒻capb(k₂, mocb, mocb₀)
    fepb = 𝒻fepb(k₃, O, O₀, P, P₀, mocb, mocb₀)
    oxidw = 𝒻oxidw(k₇)
    #evaluate dP/dt
    dP = phosw - mopb - capb - fepb
    #evaluate dA/dt
    dA = -mocb - W + oxidw + V
    return dP, dA
end

function precopse!(d, P, A, 𝒻W::F, params)::Nothing where {F}
    d[1], d[2] = precopse(P, A, 𝒻W, params)
    nothing
end

function precopse!(du, u, p, t)::Nothing
    @inbounds precopse!(du, u[1], u[2], p[1], p[2])
    nothing
end

#------------------------------------------------------------------------------
# integration of the model

export integrate

function integrate(t, 𝒻W::F, params::NamedTuple=initparams()) where {F}
    #initial conditions
    u₀ = Float64[params[:P₀], params[:A₀]]
    #time span
    tspan = (0.0, Float64(t))
    #bundle function with numeric parameters
    p = (𝒻W, params)
    #problem definition
    prob = ODEProblem(precopse!, u₀, tspan, p)
    #run the solver
    sol = solve(prob, Rodas4P())
    #return only the CO2 concentration and time
    sol.t, 𝒻pCO2.(sol[2,:], params[:h])
end

end
