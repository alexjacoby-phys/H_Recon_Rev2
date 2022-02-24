#########################################################################
#
#  Density Matrix Renormalization Group (and other methods) in julia (DMRjulia)
#                              v0.1
#
#########################################################################
# Made by Thomas E. Baker (2018)
# See accompanying license with this program
# This code is native to the julia programming language (v1.1.0) or (v1.5)
#
import Pkg
Pkg.add("LinearAlgebra")
Pkg.add("DelimitedFiles"):q1
using LinearAlgebra
using DelimitedFiles
path = "../"
include(join([path,"DMRjulia.jl"]))
using .DMRjulia

θ = 0.65
Ns = 10
MXM = 40 #6^Ns
psi = makePsi0(Ns)
H = XVBSmake(Int(Ns),#=Float64(1/2)=#cos(θ),#=Float64(1/12)=#0.25*sin(θ))
#@time params = dmrg(psi,H,maxm=MXM,sweeps=250,cutoff=1E-9)
nD = 100 # number of elements of sqrt(rho) to save
dummy_D = zeros(nD)
@time params = dmrg(psi,H,maxm=MXM,sweeps=100,cutoff=1E-6, storeD=dummy_D,allSvNbond=true)

bilin = XVBSmake(Int(Ns),1.,0.)
biquad = XVBSmake(Int(Ns),0.,1.)

bldis = expect(psi,bilin)
bqdis = expect(psi,biquad)
A = expect(psi,bilin,bilin)- bldis^2
B = expect(psi,biquad,biquad)-bqdis^2
C = expect(psi,bilin,biquad)-bldis*bqdis
savearray = [A,B,C]
filename = string("corrmat-theta=",θ,".txt")
touch(filename)
open(filename, "w") do io
    writedlm(filename, savearray)
end
thet = []
fluctuation = []
for i in 0:20
    K = θ-0.1+0.01*i
    X = XVBSmake(Int(Ns),cos(K),0.25*sin(K))
    Y = expect(psi,X,X)-expect(psi,X)^2
    append!(thet,K)
    append!(fluctuation,Y)
end


filename = string("Fluctuation-theta=",θ,".txt")
touch(filename)
savearray = [Vector{Float64}(thet),Vector{Float64}(fluctuation)]
open(filename,"w") do io
    writedlm(filename,savearray)
end
