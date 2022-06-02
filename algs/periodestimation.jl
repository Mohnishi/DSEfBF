"""
Period estimation 
"""

#------------------Packages-------------------#
using Base.Iterators: partition
using ElasticArrays
# data buffer for cellular automata
struct DataBuffers{T<:AbstractFloat}
    statebuf::ElasticArray{Int}      
    obsbuf::ElasticArray{T} 
    rewbuf::ElasticArray{T} 
    function DataBuffers{T}(height, width) where T<:AbstractFloat
        new(ElasticArray{Int}(undef,height,width, 0),
            ElasticArray{T}(undef,height,width, 0),
            ElasticArray{T}(undef,1, 0),
           )
    end
end
# data buffer for other types
struct DataBuffers2{T<:AbstractFloat}
    statebuf::ElasticArray{T}      
    obsbuf::ElasticArray{T} 
    rewbuf::ElasticArray{T} 
    function DataBuffers2{T}(d) where T<:AbstractFloat
        new(ElasticArray{T}(undef,d, 0),
            ElasticArray{T}(undef,d, 0),
            ElasticArray{T}(undef,1, 0),
           )
    end
end

# Alg structure
struct PeriodAlg{T,E}
    T_p::Int
    dim::Int
    eps::T
    env::E
    buffers::Union{DataBuffers, DataBuffers2}
    function PeriodAlg{T}(
                 T_p::Int,
                 dim::Int,
                 eps::T,
                 env;
                 prebuffer = nothing,
                ) where T<:AbstractFloat

        # check errors
        0 < T_p  || throw(ArgumentError("T_p must be > 0"))
        0 < eps || throw(ArgumentError("eps must be > 0"))
        0 < dim || throw(ArgumentError("dim must be > 0"))
        new{
            T,
            typeof(env)
           }(
             T_p,
             dim,
             eps,
             env,
             prebuffer == nothing ? DataBuffers{Float32}(height,width) : prebuffer,
            )
    end
end

function Rval(rewbuf::ElasticArray, ℓ::Int,α::Int,β::Int,s::Int,m::Int,T_p::Int)
    sumval = 0
    for j=0:Int(floor(T_p/β))-1
        sumval = sumval + rewbuf[T_p*(m-1)+1+β*j+s]*exp(im * 2 * π * α * j / ℓ)
    end
    return abs(sumval / (floor(T_p/β)))
end

function findperiod(rewbuf::ElasticArray, epsilon::AbstractFloat, Lmax::Int, dim::Int)
    L1=ElasticArray{Int}(undef,1, 0)
    append!(L1, 1)
    for m = 1:dim
        β = 1
        val = 1
        while true
            for s in 0:β-1, ℓ in 1:Int(floor(Lmax/β)), α in 1:ℓ-1
                R = Rval(rewbuf,ℓ,α,β,s,m,T_p)
                if R > epsilon
                    val = Int(ℓ/(gcd(α,ℓ)))
                    break
                end
            end
            β = β * val
            if (val == 1) || (floor(Lmax/β) < 2)
                break
            end
            val = 1
        end
        append!(L1, β)
    end
    return lcm(L1)
end
#Iterator 
function Base.iterate(estimation::PeriodAlg{DT}, i = 1) where {DT}
    T_p, dim, eps, env, buffers = 
    estimation.T_p, estimation.dim, estimation.eps, estimation.env, estimation.buffers
  

    if (i == 1) || (mod(i, 100000) == 0); @info "Iterations:" i; end;

    if (i<30) # save some trajectories
        append!(buffers.statebuf, env.state)
        append!(buffers.obsbuf, env.obs)
    end
    step!(env)
    curr_dim = min(Int(ceil(i / T_p)), dim)
    est = 0
    rew = getreward(curr_dim, env)
    append!(buffers.rewbuf, rew)
    if i == T_p * dim 
        est = findperiod(buffers.rewbuf, eps, Lmax, dim)
    end
    #------------------updating state---------------#
    result = (
              iter = i,
              rew = rew,
              est = est,
             )

    return result, i + 1
end
