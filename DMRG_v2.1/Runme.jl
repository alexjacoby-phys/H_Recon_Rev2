include("Hamiltonian-maker.jl")


Ns = 20
H = XVBSmake(Ns,0.1,0.1)
psi = makePsi0(Ns)


dmrg(psi, H, maxm = 40, cutoff = 1E-9, noise = 1E-9, allSvNbond = true)
move!
