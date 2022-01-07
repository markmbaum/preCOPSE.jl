using preCOPSE
using Roots
using PyPlot
using OrdinaryDiffEq: ContinuousCallback

pygui(true)

const millscap = 0.3*𝐒ₑ*0.125

## find the pCO2 where whak hits the cap

tol = big(1e-19)
const pCO2cap = Float64(
    find_zero(
        pCO2 -> whak₀(pCO2) - millscap,
        (big(0.0001), big(0.1)),
        strict=true,
        xatol=tol,
        xrtol=tol,
        atol=tol,
        rtol=tol
    )
)

## use capped whak for the Mills function

mills(pCO2) = pCO2 > pCO2cap ? millscap : whak₀(pCO2)

## parameters and time samples

t = exp10.(LinRange(-2, 8, 10000))
params = initparams(A₀=1.28e20)

##

figure()

semilogx(t, integrate(t, whak₀, params), label="WHAK")
semilogx(t, integrate(t, mac₀, params), label="MAC")

condition(u, t, integrator) = 𝒻pCO2(u[2], params[:h]) - pCO2cap
affect!(integrator) = nothing
cb = ContinuousCallback(condition, affect!, save_positions=(false,false))
semilogx(t, integrate(t, mills, params, callback=cb), label="Mills")

xlabel("Time [yr]")
ylabel("pCO2 [bar]")
legend()