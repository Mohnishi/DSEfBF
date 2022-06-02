"""
Period estimation LifeGame
"""
#------------------Packages-------------------#
using LinearAlgebra, Random, Statistics, UniversalLogger
using LyceumBase.Tools, JLSO

Random.seed!(seednum)

include("../algs/periodestimation.jl")
include("../utils/ca.jl")

#------------------Functions-------------------#
function SampleComputation(d::Int, Lmax::Int, B::AbstractFloat, μ::AbstractFloat,
    R::AbstractFloat, rho::AbstractFloat, delta::AbstractFloat, eps::AbstractFloat)
    r = 1.01*(μ/eps)
    r < 1  || throw(ArgumentError("r error"))
    A = R^2*log(4*Lmax^2*d*log(Lmax)/delta)
    T_p = ceil((48 * d * Lmax^2 * A) / (rho^2*(1-r)^2) + 
        (108*B*sqrt(d)*Lmax^3) / (rho*(1-r)) )
    return Int(T_p)
end
function EpsComputation(rho, d, L_max)
    return rho/(6*sqrt(d)*Lmax)
end

#------------------constants-------------------#
height = 12
width = 12
d = 5#(height) * (width)
rho =0.98
delta = 0.2
Lmax = 10
R = 0.3
μ = 0.
B = sqrt(d)
eps = EpsComputation(rho, d, Lmax)
T_p = SampleComputation(d, Lmax, B, μ, R, rho, delta, eps)
#T_p = 10000
env = LifeGame{Float64}(width, height, R)
#------------------Algorithm-------------------#
period_alg = PeriodAlg{Float64}(
        T_p,
        d,
        eps,
        env
    )


#------------------Main-------------------#
function ca_alg(alg::PeriodAlg; NITER=1000)
    # save data to the following file
    exper = Experiment("log/period/period.jlso", overwrite = false)

    lg = ULogger()
    for (i, state) in enumerate(alg)
        if i >= NITER
            exper[:buffers] = alg.buffers
            break
        end
        # save all the log data
        push!(lg, :algstate, filter_nt(state, exclude = ()))
    end
    exper, lg
end

exper, lg = ca_alg(period_alg; NITER=T_p*d+1);

exper[:logs] = get(lg)
finish!(exper); # flushes everything to disk
    