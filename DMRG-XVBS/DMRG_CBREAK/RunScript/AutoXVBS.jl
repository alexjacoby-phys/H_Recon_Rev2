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
Pkg.add("DelimitedFiles")

using LinearAlgebra
using DelimitedFiles

path = "../"
include(join([path,"DMRjulia.jl"]))
using .DMRjulia
function SavePsi(psi,filename::String)
#must have DelimitedFiles package in environment
    mkdir(filename)
    A = string(filename,"/size.txt")
    B = string(filename,"/T.txt")
    touch(A)
    touch(B)
    Ns = length(psi)
    sizes = []
    Ts = []
    for i in 1:Ns
        push!(sizes,psi[i].size)
        push!(Ts,psi[i].T)
    end
    open(A,"w") do io
        writedlm(A,sizes)
    end
    open(B,"w") do io
        writedlm(B,Ts)
    end
end


function ReadPsi(filename)
    # requires both DelimitedFiles and MPutil.jl from DMRjulia to be loaded ahead of time
    folder = string(pwd(),"/",filename)
    sizes = open(string(filename,"/size.txt")) do file
         readdlm(file)
    end
    Ts = open(string(filename,"/T.txt")) do file
         readdlm(file)
    end
    Ns = length(sizes[:,3])
    psi = []
    for i in 1:Ns
        push!(psi, tens(Vector{Int64}(sizes[i,:]), Vector{Float64}(deleteat!(Ts[i,:],Ts[i,:] .== "" ))))
    end
    psi_typed = Vector{tens}(psi)
    return MPS(psi_typed)
end

fn  ="L=30-M=60at-theta0-10sweeps"
string(fn,"/T.txt")

function UpdatePsi(psi,filename::String)
    A = string(filename,"/size.txt")
    B = string(filename,"/T.txt")

end

θarray = [0]

# for θ in θarray
#     Ns = 5
#     MXM = 10 #6^Ns
#     psi = makePsi0(Ns)
#     H = XVBSmake2(Int(Ns),cos(θ),0.25*sin(θ),zeros(Float64,4,4))
#     #@time params = dmrg(psi,H,maxm=MXM,sweeps=250,cutoff=1E-9)
#     nD = MXM # number of elements of sqrt(rho) to save
#     dummy_D = zeros(nD)
#     @time params = dmrg(psi,H,maxm=MXM,sweeps=1,cutoff=1E-8, storeD=dummy_D,allSvNbond=true)
# end

θ = 0.
Ns = 30
MXM = 60 #6^Ns
psi = ReadPsi("L=30-M=60at-theta0-10sweeps")#makePsi0(Ns)
H = XVBSmake2(Int(Ns),cos(θ),0.25*sin(θ),zeros(Float64,4,4))
#@time params = dmrg(psi,H,maxm=MXM,sweeps=250,cutoff=1E-9)
nD = MXM # number of elements of sqrt(rho) to save
dummy_D = zeros(nD)
@time params = dmrg(psi,H,maxm=MXM,sweeps=10,cutoff=1E-9, storeD=dummy_D,allSvNbond=false#=true=#,cvgE = true, noise = 0., noise_incr = 0.)
expect(psi,H)


expect(ReadPsi("L=30-M=60at-theta0-20sweeps"),H)
expect(ReadPsi("L=30-M=60at-theta0-10sweeps"),H)

# Y1 = params.SvNvec
# Y2 = params.Dvec
# popfirst!(Y1)
# popfirst!(Y2)
# savearray = [θ,Vector(2:Ns),Y1,Y2]
# filename = string("DMRG_EE_ThetaIt=",i,"_Sites=",Ns,".txt")
# touch(filename)
# open(filename,"w") do io
#     writedlm(filename, savearray)
# end
