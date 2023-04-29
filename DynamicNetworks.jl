# function to obtain reproductive value of each individual, i.e. pi_{i}^{beta} (see Equation 4 in the main text)
# input -----
# edge_seq: edge list 
# transition: network transition matrix 
# output -----
# an N-by-L matrix, where the element in the row i and column beta is pi_{i}^{beta}
function pi_DB(edge_seq::Array{Float64,2}, transition::Array{Float64,2})
    # network size
    N = Int64(maximum(edge_seq[:,1:2]))
    # number of networks
    L = Int64(maximum(edge_seq[:,4]))
    P = spzeros(Float64, L*N, N)
    for i = 1:L
        edge_list_i = findall(isequal(i), edge_seq[:,4])
        M_i = sparse(edge_seq[edge_list_i,1], edge_seq[edge_list_i,2], edge_seq[edge_list_i,3], N, N)
        deg = sum(M_i, dims = 1)
        P[(i-1)*N+1:i*N,:] = M_i./deg
    end
    Q = sparse(transition[:,1], transition[:,2], transition[:,3], L, L)
    MatA = spzeros(Float64, L*N, L*N)
    for i = 1:L
        for j = 1:L
            MatA[(i-1)*N+1:i*N, (j-1)*N+1:j*N] = Diagonal(Q[i,j]*ones(Float64,N))
        end
    end
    for i = 1:L
        Pbeta = P[(i-1)*N+1:i*N,:]/N
        for j = 1:L
            MatA[(i-1)*N+1:i*N, (j-1)*N+1:j*N] = MatA[(i-1)*N+1:i*N, (j-1)*N+1:j*N] + Pbeta*Q[i,j]
        end
    end
    for i = 1:L
        Pbeta = sum(P[(i-1)*N+1:i*N,:]/N, dims = 1)
        for j = 1:L
            MatA[(i-1)*N+1:i*N, (j-1)*N+1:j*N] = MatA[(i-1)*N+1:i*N, (j-1)*N+1:j*N] - Diagonal(Pbeta[1,:])*Q[i,j]
        end
    end
    MatA = MatA - Diagonal(ones(Float64,L*N))
    MatB = zeros(Float64, L*N)
    seq1 = [i*N for i = 1:L]
    seq2 = setdiff(1:L*N, seq1)
    MatB = sum(MatA[:,seq1], dims = 2)
    MatB_reduced = -MatB[seq2]
    for i = 1:L
        MatA[:,(i-1)*N+1:i*N-1] = MatA[:,(i-1)*N+1:i*N-1].- MatA[:,i*N]
    end
    MatA_reduced = MatA[seq2,seq2]
    pi_solution_reduced = zeros(Float64, L*(N-1))
    if sum(Diagonal(Q)) == L
        for i = 1:L
            pi_solution_reduced[(i-1)*(N-1)+1:i*(N-1)] = MatA_reduced[(i-1)*(N-1)+1:i*(N-1),(i-1)*(N-1)+1:i*(N-1)]\MatB_reduced[(i-1)*(N-1)+1:i*(N-1)]
        end
    else
        pi_solution_reduced = idrs(MatA_reduced,MatB_reduced)
    end
    pi_solution = zeros(Float64,L*N)
    for i = 1:L
        pi_solution[(i-1)*N+1:i*N-1] = pi_solution_reduced[(i-1)*(N-1)+1:i*(N-1)]
        pi_solution[i*N] = 1-sum(pi_solution_reduced[(i-1)*(N-1)+1:i*(N-1)])
    end
    return reshape(pi_solution,:,L)
end

# function to obtain eta_ij (see Equation 5 in the main text)
# input -----
# edge_seq: edge list 
# transition: network transition matrix 
# output -----
# an NL-by-N matrix, where the element in the row (beta-1)*N+i and column j is eta_{ij}^{beta}
function eta_DB(edge_seq::Array{Float64,2}, transition::Array{Float64,2})
    N = Int64(maximum(edge_seq[:,1:2]))
    L = Int64(maximum(edge_seq[:,4]))
    P = spzeros(Float64, L*N, N)
    for i = 1:L
        edge_list_i = findall(isequal(i), edge_seq[:,4])
        M_i = sparse(edge_seq[edge_list_i,1], edge_seq[edge_list_i,2], edge_seq[edge_list_i,3], N, N)
        deg = sum(M_i, dims = 1)
        P[(i-1)*N+1:i*N,:] = M_i./deg
    end
    Q = zeros(Float64, L, L) + sparse(transition[:,1], transition[:,2], transition[:,3], L, L)
    vec_solution = ones(Float64,L)/L
    if sum(Diagonal(Q)) != L
        MatA = transpose(Q)-Diagonal(ones(Float64,L))
        MatB = -MatA[:,L]
        MatB_reduced = MatB[1:L-1]
        MatA = MatA .- MatA[:,L]
        MatA_reduced = MatA[1:L-1,1:L-1]
        vec_solution_reduced = MatA_reduced\MatB_reduced
        vec_solution[1:L-1] = vec_solution_reduced
        vec_solution[L] = 1-sum(vec_solution_reduced)
    end
    # \tau_{ij}
    MatA = spzeros(Float64, L*N^2, L*N^2)
    for i = 1:L
        for j = 1:L
            Pgamma = transpose(P[(j-1)*N+1:j*N,:]*Q[j,i]/N)
            inds = findall(!iszero, Pgamma)
            a = getindex.(inds, 1)
            b = getindex.(inds, 2)
            vec = [i for i = 1:N]
            vec = reshape(vec,1,N)
            X1 = N*(a.-1)*ones(Int64,1,N)+ones(Int64,length(a))*vec
            Y1 = N*(b.-1)*ones(Int64,1,N)+ones(Int64,length(b))*vec
            W1 = Pgamma[inds]*ones(Int64,1,N)
            X11 = reshape(X1,:,1)
            Y11 = reshape(Y1,:,1)
            W11 = reshape(W1,:,1)
            MatA[(i-1)*N^2+1:i*N^2,(j-1)*N^2+1:j*N^2] = MatA[(i-1)*N^2+1:i*N^2,(j-1)*N^2+1:j*N^2] + sparse(X11[:,1],Y11[:,1],W11[:,1],N^2,N^2)
            vec = [(i-1)*N for i = 1:N]
            vec = reshape(vec,1,N)
            X2 = a*ones(Int64,1,N)+ones(Int64,size(a,1))*vec
            Y2 = b*ones(Int64,1,N)+ones(Int64,size(b,1))*vec
            W2 = Pgamma[inds]*ones(Int64,1,N)
            X22 = reshape(X2,:,1)
            Y22 = reshape(Y2,:,1)
            W22 = reshape(W2,:,1)
            MatA[(i-1)*N^2+1:i*N^2,(j-1)*N^2+1:j*N^2] = MatA[(i-1)*N^2+1:i*N^2,(j-1)*N^2+1:j*N^2] + sparse(X22[:,1],Y22[:,1],W22[:,1],N^2,N^2)
        end
    end
    for i = 1:L
        for j = 1:L
            MatA[(i-1)*N^2+1:i*N^2,(j-1)*N^2+1:j*N^2] = MatA[(i-1)*N^2+1:i*N^2,(j-1)*N^2+1:j*N^2] + Diagonal(ones(Float64,N^2))*Q[j,i]
            Pgamma = sum(P[(j-1)*N+1:j*N,:]*Q[j,i]/N,dims=1)
            Pgamma2 = ones(Float64, N)*Pgamma + transpose(ones(Float64, N)*Pgamma)
            Pgamma3 = reshape(Pgamma2,:,1)
            MatA[(i-1)*N^2+1:i*N^2,(j-1)*N^2+1:j*N^2] = MatA[(i-1)*N^2+1:i*N^2,(j-1)*N^2+1:j*N^2] - Diagonal(Pgamma3[:,1])
        end
    end
    MatA = MatA - Diagonal(ones(Float64,L*N^2))
    vec = [(i-1)*N+i+(j-1)*N^2 for i = 1:N, j = 1:L]
    vec2 = reshape(vec,:,1)
    vec3 = setdiff(1:L*N^2, vec2)
    MatA_reduced = MatA[vec3,vec3]
    MatB = reshape(ones(Float64, N^2)*transpose(vec_solution),:,1)
    MatB = 2*MatB/N
    MatB_reduced = MatB[vec3]
    eta_solution_reduced = idrs(MatA_reduced,MatB_reduced)
    eta_solution = spzeros(Float64, L*N^2)
    eta_solution[vec3] = eta_solution_reduced
    eta_solution = reshape(eta_solution, N, :)
    return transpose(eta_solution)
end

# function to obtain u0, v0, u2, v2, and the critical benefit-to-cost ratio (see Equations 10 and 11 in the main text)
# input -----
# edge_seq: edge list 
# transition: network transition matrix 
# output -----
# critical benefit-to-cost ratio (b/c)^*
function bc_DB(edge_seq::Array{Float64,2}, transition::Array{Float64,2})
    N = Int64(maximum(edge_seq[:,1:2]))
    L = Int64(maximum(edge_seq[:,4]))
    P = spzeros(Float64, L*N, N)
    M = spzeros(Float64, L*N, N)
    for i = 1:L
        edge_list_i = findall(isequal(i), edge_seq[:,4])
        M_i = sparse(edge_seq[edge_list_i,1], edge_seq[edge_list_i,2], edge_seq[edge_list_i,3], N, N)
        M[(i-1)*N+1:i*N, :] = M_i
        deg = sum(M_i, dims = 1)
        P[(i-1)*N+1:i*N,:] = M_i./deg
    end
    Q = zeros(Float64, L, L) + sparse(transition[:,1], transition[:,2], transition[:,3], L, L)
    pi_solution = pi_DB(edge_seq,transition)
    eta_solution = eta_DB(edge_seq,transition)
    u0 = 0.0
    for i = 1:L
        P_temp1 = transpose(P[(i-1)*N+1:i*N,:])
        P_temp2 = (sum(transpose(Q[i,:]).*pi_solution,dims=2).*P_temp1)*(transpose(M[(i-1)*N+1:i*N,:]).*eta_solution[(i-1)*N+1:i*N,:])/N
        u0 = u0 + sum(P_temp2)
    end
    v0 = 0.0
    u2 = 0.0
    for i = 1:L
        P_temp1 = transpose(P[(i-1)*N+1:i*N,:])
        P_temp2 = (sum(transpose(Q[i,:]).*pi_solution,dims=2).*P_temp1).*(P_temp1*transpose(M[(i-1)*N+1:i*N,:])*transpose(eta_solution[(i-1)*N+1:i*N,:]))/N
        u2 = u2 + sum(P_temp2)
    end
    v2 = 0.0
    for i = 1:L
        P_temp1 = transpose(P[(i-1)*N+1:i*N,:])
        P_temp2 = (sum(transpose(Q[i,:]).*pi_solution,dims=2).*P_temp1).*(P_temp1*(sum(M[(i-1)*N+1:i*N,:],dims=2).*transpose(eta_solution[(i-1)*N+1:i*N,:])))/N
        v2 = v2 + sum(P_temp2)
    end
    return (v0-v2)/(u0-u2)
end


# main function 
# calculate the benefit-to-cost ratio (b/c)^* for the two-clique network in Figure 2a
################################################################################
#---
using Printf
using Statistics
using LinearAlgebra
using DelimitedFiles
using TickTock
using IterativeSolvers
using SparseArrays
# rng = MersenneTwister(1234)

# input the edge list of a two-clique network with size N = 40 and a = 0.5
# elements in each row: [the source node id, the target node id, the edge weight, the network id]
edge_seq = readdlm("edge_list1", '\t', Int64)

# set the transition pattern
# elements in each row: [the current network id, the next network id, the transition probability]
N = 40
transition = [1 1 1-1/N;
              1 2 1/N;
              2 1 1/N;
              2 2 1-1/N]

# calculate the critical benefit-to-cost ratio
bc = bc_DB(1.0*edge_seq, transition)
println("For the two-clique network with N = 40 and a = 0.5, the critical benefit-to-cost ratio (b/c)* is ", bc)

# input the edge list of a two-clique network with size N = 40 and a = 0.7
edge_seq = readdlm("edge_list2", '\t', Int64)
bc = bc_DB(1.0*edge_seq, transition)
println("For the two-clique network with N = 40 and a = 0.7, the critical benefit-to-cost ratio (b/c)* is ", bc)
