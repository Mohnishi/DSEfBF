"""
Eigen estimation 
"""
#------------------Packages-------------------#
using LinearAlgebra, Random, Statistics, UniversalLogger
using LyceumBase.Tools, JLSO

Random.seed!(seednum)

include("../algs/eigenestimation.jl")
include("../utils/linear.jl")

#------------------Functions-------------------#
function SampleComputation(d::Int, q::Int, B::AbstractFloat, R::AbstractFloat, 
    delta::AbstractFloat, Δ::AbstractFloat, κ::AbstractFloat, C::AbstractFloat)
    #n = max( ( (4*d-3)/Δ )^2 + ( (5*d+4)^2 * log(3*d*(κ+1)) /Δ), 16q^2)
    n = max( (-(d+1) * log(Δ) + log(B*κ^2) + d + 6)/Δ + d, 16q^2)
    return Int(ceil(n*C))
end
function EpsComputation(d, R, delta, N)
    return (sqrt(4*d^2*R^2*log(4*d^2/delta))+1)/sqrt(N)
end

#------------------constants-------------------#
#ConstN = 5
d = 5
delta = 0.2
R = 0.3
q = 24
κ = 6.
B = 1.
Δ = 0.1
N = SampleComputation(d, q, B, R, delta, Δ, κ, ConstN)
T_s = d + N + 1
T_e = d + N + 2 * N * d^2
@info T_e
#T_p = 10000
eps = EpsComputation(d, R, delta, N)
M=[0 0 0 1 0;0 1 0 0 0;1 0 0 0 0;0 0 1 0 0;0 0 0 0 0.7]
env = LinearDyn{Float64}(d, M, R, noise = Uniform(-R, R))
#------------------Algorithm-------------------#
lin_alg = LinAlg{Float64}(
        T_s,
        T_e,   
        d,
        q,
        κ,
        Δ,
        N,
        eps,
        env
    )


#------------------Main-------------------#
function linear_alg(alg::LinAlg; NITER=1000)
    # save data to the following file
    exper = Experiment("log/eigen/eigen.jlso", overwrite = false)

    lg = ULogger()
    for (i, state) in enumerate(alg)
        if i >= NITER
            exper[:buffers] = alg.buffers
            break
        end
        # save all the log data
        #push!(lg, :algstate, filter_nt(state, exclude = ()))
    end
    exper, lg
end

exper, lg = linear_alg(lin_alg; NITER=T_e+2);

exper[:logs] = get(lg)
finish!(exper); # flushes everything to disk
    