import LinearAlgebra
using DMRJtensor

# makes an initial MPS for spin s and size Ns
function makePsi0(s::Float64,Ns::Integer)

    n = Int64(2*s+1) # size of matrices
    initTensor = [zeros(ComplexF64,1,n,1) for i=1:Ns]
    for i = 1:Ns 
        #initTensor[i][1,1,1] = 1 # spin-up start
        for j = 1:n   # random start?
            initTensor[i][1,j,1] = rand(1)[1] + 1im*rand(1)[1]
        end
    end
    psi = MPS(ComplexF64,initTensor,1)
    return psi
end

# iteratively constructs the MPO for a given "power" of spin terms, e.g. power 3 -> \sum_{ijk} S_i S_j S_k
function S_pow_iterative(H, spin_mag::Float64, curr_op, start_index::Integer, npowers::Integer, curr_power::Integer, J)
    # npowers is the number of terms in the spin sum, e.g. 3 in the above example
    # curr_power is where we are in the nested sum, e.g. S_i -> 1, S_j -> 2, and S_k -> 3 in above example
    # start_index keeps track of where we are in the larger MPO (which is a direct product of different "npowers" calls to this iterative function)
    # curr_op just keeps track of the preceding terms in the spin product, e.g. for S_k (curr_power = 3) curr_op would be S_i S_j.
    
    Sx,Sy,Sz,Sp,Sm,O,Id = spinOps(spinmag)
    Sop_arr = [Sx,Sy,Sz] # spin operators
    nSpin = Int(2*spin_mag + 1) # dimension of spin space
    nterms = size(H,1) # total size of the MPO's onsite H
    
    if (curr_power == npowers) # if at bottom of iterative tree
                
        for dim = 1:3 # loop over (x,y,z)
            
            idx_start = start_index + nSpin*(dim-1) # starting point along first column / last row
            lr_start = nterms - nSpin + 1 # index location of last row

            op_h = curr_op*Sop_arr[dim]; # generate final operator
            
            H[idx_start:(idx_start+nSpin-1),1:nSpin] = -J*op_h # first column (has J!)
                
            H[lr_start:(lr_start+nSpin-1),idx_start:(idx_start+nSpin-1)] = op_h # last column (no J!)
                
        end
    else # otherwise, keep iterating
        for dim = 1:3 # loop over (x,y,z)
            op_h = curr_op*Sop_arr[dim] # add term to operator product
            stride_idx = 3^(npowers-curr_power) # figure out how wide the next block will be
            # iterate!
            S_pow_iterative(H, spin_mag, op_h, start_index + nSpin*stride_idx*(dim-1), npowers, curr_power+1, J) 
        end
    end
        
end

# construct the Heisenberg chain MPO
function H_SU2(s::Float64, J_arr)
    # J_arr is an array of prefactors -> H = \sum_p J_p (S*S)^p
             
    Sx,Sy,Sz,Sp,Sm,O,Id = spinOps(spinmag)
    n = Int(2*s + 1) # dimension of spin space
    max_p = n-1 # max allowed power of (S*S)^p terms
    nterms = n*2 # keep track of Id and O terms in MPO
    for p = 1:max_p # find size over all allowed powers
        nterms = nterms + n*(3^p)
    end
    H_n = zeros(ComplexF64, nterms, nterms) # initialize
                     
    # set first column and last row identities
    H_n[1:n,1:n] = Id
    H_n[(nterms-n+1):nterms,(nterms-n+1):nterms] = Id
                        
    # we start in H just after the first [n x n] row block (or column block, if on last row)
    start_index = 1 + n;
    
    for p = 1:max_p # add terms associated with each power, (S*S)^p
        S_pow_iterative(H_n, s, Id, start_index, p, 1, J_arr[p]) # iterative construction of (S*S)^p
        start_index = start_index + n*(3^p) # keep track of the H index for each (S*S)^p group!
    end
    return H_n 
end
