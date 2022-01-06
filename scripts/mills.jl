using preCOPSE
using PyPlot

pygui(true)

##

mills(pCO2) = min(whak₀(pCO2), 0.3*𝐒ₑ*0.125)

##

t = exp10.(LinRange(0, 8, 1000))
params = initparams(A₀=1.28e20)

##

figure()
semilogy(t, integrate(t, whak₀, params), label="WHAK")
semilogy(t, integrate(t, mac₀, params), label="MAC")
#semilogx(t, integrate(t, mills, params), label="Mills")

xlabel("Time [yr]")
ylabel("pCO2 [bar]")
legend()