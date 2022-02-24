import LinearAlgebra
using DMRJtensor

# returns Sx,Sy,Sz for given spin s, in a basis of diagonalized Sz
function SU2_spinOps(s::Float64=0.5) 
    # make sure spin is integer/half-integer
    if (mod(2*s,1) != 0)
        throw(DomainError(s,"spin must be integer or half integer!"))
    end
    
    n = Int64(2*s+1) # size of matrices
    Sx = LinearAlgebra.zeros(Complex, n, n) 
    Sy = LinearAlgebra.zeros(Complex, n, n) 
    Sz = LinearAlgebra.zeros(Complex, n, n)
    for d = 1:n # iterate over matric indices
        m = s - (d-1) # "spin" index for state, e.g. |s,m >
        am = sqrt(s*(s+1) - m*(m-1)) # S_- operator
        Sz[d,d] = m # Sz is diagonal
        if (d < n)
            Sx[d,d+1] = am # Sx and Sy both have S_- as element
            Sy[d,d+1] = am
        end
    end
    Sx = 0.5*(Sx + transpose(Sx)) # make rest of Sx
    Sy = complex(0,1)*0.5*(Sy - transpose(Sy)) # make rest of Sy

    return Sx,Sy,Sz
    
end

function makePsi0(Ns::Integer)
  # todo
end
