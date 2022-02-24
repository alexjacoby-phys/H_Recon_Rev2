import LinearAlgebra
using DMRJtensor


function xVBSOps()
  #SU(4) generators + Relevant Operators for VBS
  QS = 6
  N = LinearAlgebra.zeros(Float64, QS, QS)
  Gen = LinearAlgebra.zeros(Float64,4,4,QS,QS)
  #Cartan Subalgebra (Rnk 3+1-- not linearly independent, but doesn't matter)
  Gen[1,1,1,1] = Gen[1,1,2,2] = Gen[1,1,3,3] = 0.5
  Gen[1,1,4,4] = Gen[1,1,5,5] = Gen[1,1,6,6] = -0.5
  #####
  Gen[2,2,1,1] = Gen[2,2,4,4] = Gen[2,2,5,5] = 0.5
  Gen[2,2,2,2] = Gen[2,2,3,3] = Gen[2,2,6,6] = -0.5
  #####
  Gen[3,3,1,1] = Gen[3,3,3,3] = Gen[3,3,5,5] = 0.5
  Gen[3,3,2,2] = Gen[3,3,4,4] = Gen[3,3,6,6] = -0.5
  #####
  Gen[4,4,3,3] = Gen[4,4,5,5] = Gen[4,4,6,6] = 0.5
  Gen[4,4,1,1] = Gen[4,4,2,2] = Gen[4,4,4,4] = -0.5
  #####
  #The rest of the generators
  ##### first block in c code
  Gen[1,2,4,2] = Gen[2,1,2,4] = Gen[1,2,5,3] = Gen[2,1,3,5] = 1.0
  ##### third block in c code
  Gen[1,4,5,1] = Gen[4,1,1,5] = Gen[1,4,6,2] = Gen[4,1,2,6] = -1.0
  ##### fourth block in c code
  Gen[2,3,2,1] = Gen[3,2,1,2] = Gen[2,3,6,5] = Gen[3,2,5,6] = 1.0
  ##### sixth block in c code
  Gen[3,4,3,2] = Gen[4,3,2,3] = Gen[3,4,5,4] = Gen[4,3,4,5] = 1.0
  ##### second block in c code
  Gen[1,3,4,1] = Gen[3,1,1,4] = 1.0
  Gen[1,3,6,3] = Gen[3,1,3,6] = -1.0
  ##### fifth block in c code
  Gen[2,4,3,1] = Gen[4,2,1,3] = 1.0
  Gen[2,4,6,4] = Gen[4,2,4,6] = -1.0
 ##### Manual reshape of generators-- reshape function didn't work
  indx = [[1,1],[1,2],[1,3],[1,4],[2,1],[2,2],[2,3],[2,4],[3,1],[3,2],[3,3],[3,4],[4,1],[4,2],[4,3],[4,4]]
  Generators = fill(N,16)
  for i = 1:16
    Generators[i] = Gen[indx[i][1],indx[i][2],:,:]
  end
  ###### Bilinear is already covered
  ###### Biquadratic
  BiQuadLeft = []
  BiQuadRight = []
  for i in 1:16, j in 1:16
    if Generators[i]*Generators[j] !== LinearAlgebra.zeros(Float64,6,6)
      push!(BiQuadLeft,Generators[i]*Generators[j])
    end
    if LinearAlgebra.transpose(Generators[i])*LinearAlgebra.transpose(Generators[j]) !== LinearAlgebra.zeros(Float64,6,6)
      push!(BiQuadRight,LinearAlgebra.transpose(Generators[i])*LinearAlgebra.transpose(Generators[j]))
    end
  end
  Id = copy(N)+LinearAlgebra.I
  #removes zero elements to fix an annoying bug
  BiQuadLeftPop = copy(BiQuadLeft)
  BiQuadRightPop = copy(BiQuadRight)
  for i in 1:length(BiQuadLeftPop)
      if iszero(BiQuadLeftPop[1+length(BiQuadLeftPop)-i]) | iszero(BiQuadRightPop[1+length(BiQuadLeftPop)-i])
          popat!(BiQuadRight,1+length(BiQuadLeftPop)-i)
          popat!(BiQuadLeft,1+length(BiQuadLeftPop)-i)
      end
  end
  return Gen, Id, Generators, BiQuadLeft, BiQuadRight
end
Gen , Id, Generators, BiQuadLeft, BiQuadRight = xVBSOps()


function Z1(a::Float64,b::Float64)
  G = length(Generators)
  BQ = length(BiQuadLeft)
  FirstSite = LinearAlgebra.zeros(Float64,1,6,6,1+G+BQ)
  for i in 1:G
    FirstSite[1,:,:,i] = a*Generators[i]
  end
  for i in 1:BQ
    FirstSite[1,:,:,G+i]= b* BiQuadLeft[i]
  end
  FirstSite[1,:,:,1+G+BQ] = Id
  FirstSiteDone = tens(FirstSite)
  return FirstSiteDone
end

function Z5()
  G = length(Generators)
  BQ = length(BiQuadLeft)
  LastSite = LinearAlgebra.zeros(Float64,1+G+BQ,6,6,1)
  LastSite[1,:,:,1] = Id
  for i in 1:G
    LastSite[1+i,:,:,1] = LinearAlgebra.transpose(Generators[i])
  end
  for i in 1:BQ
    LastSite[1+G+i,:,:,1]=BiQuadRight[i]
  end
  LastSiteDone = tens(LastSite)
  return LastSiteDone
end

############################################### OPLIST for next to ends of chains
function Z2(a::Float64,b::Float64)
  G = length(Generators)
  BQ = length(BiQuadLeft)
  NextToFirstSite = LinearAlgebra.zeros(Float64,1+G+BQ,6,6,2+G+BQ)
  NextToFirstSite[1+G+BQ,:,:,2+G+BQ] = Id
  for i in 1:G
    NextToFirstSite[i,:,:,1]  = LinearAlgebra.transpose(Generators[i])
  end
  for i in 1:BQ
    NextToFirstSite[G+i,:,:,1]  = BiQuadRight[i]
  end
  for i in 1:G
      NextToFirstSite[1+G+BQ,:,:,1+i] = a*Generators[i]
  end
  for i in 1:BQ
    NextToFirstSite[1+G+BQ,:,:,i+G+1] = b*BiQuadLeft[i]
  end
  NextToFirstSiteDone = tens(NextToFirstSite)
  return NextToFirstSiteDone
end
function Z4(a::Float64,b::Float64)
  G = length(Generators)
  BQ = length(BiQuadLeft)
  NextToLastSite = LinearAlgebra.zeros(Float64,2+G+BQ,6,6,1+G+BQ)
  NextToLastSite[1,:,:,1] = Id
  for i in 1:G
    NextToLastSite[i+1,:,:,1] = transpose(Generators[i])
  end
  for i in 1:BQ
    NextToLastSite[i+G+1,:,:,1] = BiQuadRight[i]
  end
  for i in 1:G
    NextToLastSite[G+BQ+2,:,:,i+1] = a*Generators[i]
  end
  for i in 1:BQ
    NextToLastSite[G+BQ+2,:,:,i+G+1] = b*BiQuadLeft[i]
  end
  NextToLastSiteDone = tens(NextToLastSite)
  return NextToLastSiteDone
end
############################################### OPLIST for middle sites
function Z3(a::Float64,b::Float64)
  G = length(Generators)
  BQ = length(BiQuadLeft)
  MiddleSite = LinearAlgebra.zeros(Float64,G+BQ+2,6,6,G+BQ+2)
  MiddleSite[1,:,:,1] = Id
  MiddleSite[G+BQ+2,:,:,G+BQ+2] = Id
  for i in 1:G
    MiddleSite[i+1,:,:,1] = LinearAlgebra.transpose(Generators[i])
  end
  for i in 1:BQ
    MiddleSite[i+1+G,:,:,1] = BiQuadRight[i]
  end
  for i in 1:G
    MiddleSite[G+BQ+2,:,:,1+i] = a*Generators[i]
  end
  for i in 1:BQ
    MiddleSite[G+BQ+2,:,:,i+G+1] = b*BiQuadLeft[i]
  end
  MiddleSiteDone = tens(MiddleSite)
  return MiddleSiteDone
end
############################################### OPLIST for 2 site case
function Z6(a::Float64,b::Float64)
  G = length(Generators)
  BQ = length(BiQuadLeft)
  TwoSiteL = LinearAlgebra.zeros(Float64,1,6,6,BQ+G)
  for i in 1:G
    TwoSiteL[1,:,:,i] = a*Generators[i]
  end
  for i in 1:BQ
    TwoSiteL[1,:,:,i+G] = b*BiQuadLeft[i]
  end
  TwoSiteLDone = tens(TwoSiteL)
  return TwoSiteLDone
end
function Z7()
  G = length(Generators)
  BQ = length(BiQuadLeft)
  TwoSiteR = LinearAlgebra.zeros(Float64,G+BQ,6,6,1)
  for i in 1:G
    TwoSiteR[i,:,:,1] = LinearAlgebra.transpose(Generators[i])
  end
  for i in 1:BQ
    TwoSiteR[i+G,:,:,1] = BiQuadRight[i]
  end
  TwoSiteRDone = tens(TwoSiteR)
  return TwoSiteRDone
end
############################################### OPLIST for 3 site case
function Z8(a::Float64,b::Float64)
  G = length(Generators)
  BQ = length(BiQuadLeft)
  ThreeSiteMid = LinearAlgebra.zeros(Float64,G+BQ+1,6,6,G+BQ+1)
  for i in 1:G
    ThreeSiteMid[i,:,:,1]  = LinearAlgebra.transpose(Generators[i])
  end
  for i in 1:BQ
    ThreeSiteMid[i+G,:,:,1] = BiQuadRight[i]
  end
  for i in 1:G
    ThreeSiteMid[BQ+G+1,:,:,1+i] = a*Generators[i]
  end
  for i in 1:BQ
    ThreeSiteMid[BQ+G+1,:,:,1+G+i] = b*BiQuadLeft[i]
  end
  ThreeSiteMidDone = tens(ThreeSiteMid)
  return ThreeSiteMidDone
end

G = length(Generators)
BQ = length(BiQuadLeft)
FL = G+BQ
BLL2 = zeros(Float64,6,6,G)
BLR2 = zeros(Float64,G,6,6)
BQL2 = zeros(Float64,6,6,BQ)
BQR2 = zeros(Float64,BQ,6,6)
for i in 1:G
  BLL2[:,:,i] = Generators[i]
  BLR2[i,:,:] = Generators[i]
end
for i in 1:BQ
  BQL2[:,:,i] = BiQuadLeft[i]
  BQR2[i,:,:] = BiQuadRight[i]
end












############################################### Hamiltonian Maker-- Stitches all the ops together
function XVBSmake(Ns::Integer, a::Float64, b::Float64)
   localopsarray = Vector{tens{Float64}}(undef, Ns)
   if Ns == 1
     return print("This is an error message! You need more than one site friend!")
   elseif Ns >= 5
     localopsarray[1] = Z1(a,b)
     localopsarray[Ns] = Z5()
     ###############################################
     localopsarray[2] = Z2(a,b)
     localopsarray[Ns-1] = Z4(a,b)
     ###############################################
     for i in 3:Ns-2
       localopsarray[i] = Z3(a,b)
     end
   elseif Ns == 4
     localopsarray[1] = Z1(a,b)
     localopsarray[2] = Z2(a,b)
     localopsarray[3] = Z4(a,b)
     localopsarray[4] = Z5()
   elseif Ns == 2
     localopsarray[1] = Z6(a,b)
     localopsarray[2] = Z7()
   elseif Ns == 3
     localopsarray[1] = Z1(a,b)
     localopsarray[2] = Z8(a,b)
     localopsarray[3] = Z5()
   end
   Hamiltonian = MPO(localopsarray)
  return Hamiltonian
end

function makePsi0(Ns::Integer)
  hereQS = convert(Int64,6)
  QS = cld(hereQS,2)
  initTensor = [LinearAlgebra.zeros(1,hereQS,1) for i=1:Ns]
  for i = 1:Ns
     initTensor[i][1,i%2 == 1 ? 1 : 2,1] = 1.0
  end
  psi = MPS(initTensor,1)
  return psi
end
