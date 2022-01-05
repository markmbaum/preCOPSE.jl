module preCOPSE

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

#------------------------------------------------------------------------------
# functions supporting the ODEs

#fraction of carbon in the atmososphere [-]
# A - total ocean-atmosphere CO2 [mole]
𝒻ϕ(A, h) = 0.78*A/(A + h)

#partial pressure of CO2 [bar]
# A - total ocean-atmososphere CO2 [mole]
𝒻pCO2(A, h) = (𝒻ϕ(A,h)*A*0.04401*𝐠/𝐒ₑ)/1e5

# Marine organic carbon burial 
𝒻mocb(k1,P,P₀) = k1*(P/P₀)^2

# Iron-sorbed phosphate burial 
𝒻fepb(P,P₀,O,O₀,k1,k3) = k3*(O/O₀*P₀/P)*𝒻mocb(k1,P,P₀)


function precopse!(du, u, p, t)::Nothing
    #unpack phosphate and total carbon
    P, A = u
    #unpack weathering function W(pCO2)
    𝒲    = p[1]
    #unpack all the constants
    k1    = p[2]  # mol C yr⁻¹ Total organic carbon burial
    k2    = p[3]  # mol P yr⁻¹ Ca associated phosphorus burial
    k3    = p[4]  # mol P yr⁻¹ Fe associated phosphorus burial
    k4    = p[5]  # -          Initial oxic fraction
    k5    = p[6]  # mol C yr⁻¹ Silicate weathering
    k6    = p[7]  # mol C yr⁻¹ Carbonate weathering
    k7    = p[8]  # mol C yr⁻¹ Oxidative weathering
    k8    = p[9]  # mol P yr⁻¹ Reactive phosphorus weathering
    k9    = p[10] # mol C yr⁻¹ Organic carbon degassing
    k10   = p[11] # mol C yr⁻¹ Carbonate carbon degassing
    CPsea = p[12] # mol:mol    C:P burial ratio
    A₀    = p[13] # mol        Present day atmosphere/ocean CO2
    P₀    = p[14] # mol        Present day ocean phosphate
    O₀    = p[15] # mol        Present day atmospheric oxgen
    h     = p[16] # mol        Partitioning value for pCO2
    W₀    = p[17] # mol C yr⁻¹ Reference modern weathering
    #atmospheric carbon concentration 
    pCO2 = 𝒻pCO2(A, h)
    #evaluate dP/dt
    du[1] = k8*(𝒲(pCO2)/𝒲(pCO2₀)*(7/12) + (5/12)) - mocb(P,P₀)/CPsea - k2*mocb(P,P₀)/k1 - 𝒻fepb(P,P₀,O,O₀,k1,k3)
    #evaluate dA/dt
    du[2] = k9*D + k10*D - 𝒲(pCO2) - mocb(P)
    nothing
end

end
