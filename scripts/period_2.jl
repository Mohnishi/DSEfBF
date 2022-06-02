"""
Period estimation for mu periodic
"""
#------------------Packages-------------------#
using LinearAlgebra, Random, Statistics, UniversalLogger
using LyceumBase.Tools, JLSO

Random.seed!(seednum)

include("../algs/periodestimation.jl")
include("../utils/mup.jl")

#------------------Functions-------------------#
function SampleComputation(d::Int, Lmax::Int, B::AbstractFloat, μ::AbstractFloat,
    R::AbstractFloat, rho::AbstractFloat, delta::AbstractFloat, eps::AbstractFloat)
    r = 1.01*(μ/eps)
    r < 1  || throw(ArgumentError("r error"))
    A = R^2*log(4*Lmax^2*d*log(Lmax)/delta)
    T_p = ceil((72 * d * Lmax^2 * A) / (rho^2*(1-r)^2) + 
        (108*B*sqrt(d)*Lmax^3) / (rho*(1-r)) )
    return Int(T_p)
end
function EpsComputation(rho, d, L_max)
    return rho/(6*sqrt(d)*Lmax)
end

#------------------constants-------------------#
muL = 5
ir = pi * 1
d = 2
rho =0.3
delta = 0.2
Lmax = 8
R = 0.3
μ = 0.001
B = 2.
eps = EpsComputation(rho, d, Lmax)
T_p = SampleComputation(d, Lmax, B, μ, R, rho, delta, eps)
#T_p = 10000
env = MuNearlyPeriod{Float64}(μ, ir, muL, R, noise = Uniform(-R, R))
#------------------Algorithm-------------------#
period_alg2 = PeriodAlg{Float64}(
        T_p,
        d,
        eps,
        env,
        prebuffer = DataBuffers2{Float32}(d)
    )


#------------------Main-------------------#
function mup_alg(alg::PeriodAlg; NITER=1000)
    # save data to the following file
    exper = Experiment("log/muperiod/muperiod.jlso", overwrite = false)

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

exper, lg = mup_alg(period_alg2; NITER=T_p*d+1);

exper[:logs] = get(lg)
finish!(exper); # flushes everything to disk
    