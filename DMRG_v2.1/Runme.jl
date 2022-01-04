include("Hamiltonian-maker.jl")


Ns = 20
θ = 0.3
H = XVBSmake(Ns,cos(θ),0.25*sin(θ))
psi = makePsi0(Ns)

dmrg(psi, H, maxm = 10, cutoff = 1E-3, sweeps = 100)
print("coarse run completed")
dmrg(psi, H, maxm = 50, cutoff = 1E-8, sweeps = 20, cvgE = false)
print("intermediate run completed")
dmrg(psi, H, maxm = 100, cutoff = 1E-12, sweeps = 10, cvgE = false)
print("fine run completed")
variance = expect(psi,H,H)- (expect(psi, H))^2
print("variance is, ",variance)
