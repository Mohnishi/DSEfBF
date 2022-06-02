"""
Eigen estimation 
"""

#------------------Packages-------------------#
using Base.Iterators: partition
using ElasticArrays
# data buffer
struct DataBuffersE{T<:AbstractFloat}
    statebuf::ElasticArray{T}      
    rewbuf::ElasticArray{T} 
    sumbuf1::AbstractArray{Complex}
    sumbuf2::AbstractArray{Complex}
    lastoutput::AbstractArray{Complex}
    counter::AbstractArray{Int}
    function DataBuffersE{T}(dim) where T<:AbstractFloat
        new(ElasticArray{T}(undef,dim, 0),
            ElasticArray{T}(undef,1, 0),
            zeros(Complex, dim, dim),
            zeros(Complex, dim, dim),
            zeros(Complex, dim, dim),
            zeros(Int, 2)
           )
    end
end

# Alg structure
struct LinAlg{T,E}
    T_s::Int
    T_e::Int
    dim::Int
    q::Int
    κ::T
    Δ::T
    N::Int
    eps::T
    env::E
    actionmat::Matrix{T}
    buffers::DataBuffersE    
    function LinAlg{T}(
                 T_s,
                 T_e,
                 dim,
                 q,
                 κ,
                 Δ,
                 N,
                 eps,
                 env;
                 prebuffer = nothing,
                ) where T<:AbstractFloat

        # check errors
        1 <=  T_s || throw(ArgumentError("T_s must be >= 1"))
        1 <=  T_e || throw(ArgumentError("T_e must be >= 1"))
        1 <=  N || throw(ArgumentError("N must be >= 1"))
        0 < eps || throw(ArgumentError("eps must be > 0"))
        1 <= q || throw(ArgumentError("q must be >= 1"))
        1 <= dim || throw(ArgumentError("dim must be >= 1"))
        1 <= κ || throw(ArgumentError("κ must be >= 1"))
        0 < Δ <= 1 || throw(ArgumentError("Δ must be > 0 and <= 1"))
        actionmat = randn(T, dim, dim)
        for i in 1:dim
            actionmat[:,1] .= actionmat[:,1]./norm(actionmat[:,1])
        end
        new{
            T,
            typeof(env),
           }(
             T_s,
             T_e,
             dim,
             q,
             κ,
             Δ,
             N,
             eps,
             env,
             actionmat,
             prebuffer == nothing ? DataBuffersE{Float32}(dim) : prebuffer,
            )
    end
end

function ConstructA(rew::AbstractFloat, buffers::DataBuffersE, m::Int, q::Int)
    c1 = buffers.counter[1]
    c2 = buffers.counter[2]
    α = c1%(2*d); β = mod(c1-1, d)
    if 1 <= α <= d 
        buffers.sumbuf1[m,end-β] = buffers.sumbuf1[m,end-β] +
                 exp(im * 2 * π  * c2^2 / (4*q)) * rew
    else
        buffers.sumbuf2[m,end-β] = buffers.sumbuf2[m,end-β] +
                 exp(im * 2 * π  * c2^2 / (4*q)) * rew
    end
end

function LowrankSVD(buffers::DataBuffersE, eps, N)
    buffers.sumbuf1 .= buffers.sumbuf1 ./ N
    buffers.sumbuf2 .= buffers.sumbuf2 ./ N
    A1_svd = svd(buffers.sumbuf1)
    s = A1_svd.S
    r = length(s[s .>= eps])
    if r == 0
        return zeros(Complex, size(buffers.sumbuf1))
    end
    U = A1_svd.U[:, 1:r]
    Σ = Diagonal(A1_svd.S[1:r])
    Vt = A1_svd.Vt[1:r, :]
    return buffers.sumbuf2 * (Vt' * inv(Σ) * U')
end
#Iterator 
function Base.iterate(estimation::LinAlg{DT}, i = 0) where {DT}
    T_s, T_e, dim, q, κ, Δ, N, eps, env, actionmat, buffers = 
    estimation.T_s, estimation.T_e, estimation.dim, estimation.q, 
    estimation.κ, estimation.Δ, estimation.N, estimation.eps, 
    estimation.env, estimation.actionmat, estimation.buffers


    if (i == 0) || (mod(i, 1000000) == 0); @info "Iterations:" i; end;
    step!(env)
    if i == T_e + 1
        mat = LowrankSVD(buffers, eps, N)
        buffers.lastoutput .= mat
    elseif (i >= T_s) && (i <= T_e)
        
        buffers.counter[1] = buffers.counter[1] + 1
       
        m0 = mod(i - T_s, 2 * d^2) + 1
        m = Int(ceil(m0/(2*d)))
        
        rew = getreward(actionmat[:,m], env)
        # save some trajectories
        if buffers.counter[2] < 10
            append!(buffers.statebuf, env.state)
            append!(buffers.rewbuf, rew)
        end
        
        ConstructA(rew, buffers, m, q)
        
        if buffers.counter[1] == 2 * d^2; buffers.counter[1] = 0; 
            buffers.counter[2] = buffers.counter[2] + 1; end
       
    end
    #------------------updating state---------------#
    result = (
              iter = i,
             )

    return result, i + 1
end
