include("Hamiltonian-maker.jl")


Ns = 20
H = XVBSmake(Ns,0.1,0.)
psi = makePsi0(Ns)

for i in 1:80
    @time dmrg(psi, H, maxm = 10, cutoff = 1E-4, noise = 1E-5)
end
print("coarse run completed")
for i in 1:20
    @time dmrg(psi, H, maxm = 40, cutoff = 1E-8)
end
print("intermediate run completed")
for i in 1:10
    @time dmrg(psi, H, maxm = 80, cutoff = 1E-15)
end
print("fine run completed")
variance = expect(psi,H,H)- (expect(psi, H))^2
print("variance is, ",variance)

psi[10][]
move!(psi, 5)
svd(psi)
