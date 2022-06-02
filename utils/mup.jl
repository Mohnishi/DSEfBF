"""
mu nearly periodic
"""

#------------------Packages-------------------#
using LinearAlgebra, Random, Statistics, Distributions


struct MuNearlyPeriod{N, T} 
    μ::T
    ir::T
    muL::Int
    noise::N
    state::Vector{T}
    obs::Vector{T}
    init::Vector{T}
    time::Vector{Int}

    function MuNearlyPeriod{T}(μ, ir, muL, R; noise=nothing, init=nothing) where T <: AbstractFloat
        state = zeros(T, 2)
        obs = zeros(T, 2)
        time = zeros(Int, 1)
        if init == nothing; init = [1.0 - μ/2; 0.] end
        if noise == nothing; noise = Normal(0, R) end
        env = new{typeof(noise), T}(μ, ir, muL, noise, state, obs,
                                init, time)
        reset!(env)
    end
end


function getstate!(state, env::MuNearlyPeriod) 
    state .= env.state
end
function setstate!(env::MuNearlyPeriod, state) 
    env.state .= state
end

function getobs!(obs, env::MuNearlyPeriod)
    state = zeros(size(env.state))
    getstate!(state, env) 
    obs .= state .+ rand(env.noise, size(state))
end

function getreward(index, env::MuNearlyPeriod)
    rew = env.obs[index]
    rew
end

function reset!(env::MuNearlyPeriod)
    env.state .= copy(env.init)
    getobs!(env.obs, env)
    env.time[1] = 0
    env
end
function randreset!(env::MuNearlyPeriod)
    env.state[1] .= 1.0 - env.μ/2
    env.state[2] .= 2*pi*rand()
    getobs!(env.obs, env)
    env.time[1] = 0
    env
end


function step!(env::MuNearlyPeriod)
    state, muL, μ, ir = env.state, env.muL, env.μ, env.ir
    currstate = copy(state)
    currstate[2] = currstate[2] + 2*pi/muL 
    currstate[2] = mod(currstate[2], 2*pi)
    currstate[1] = μ * (ir * (currstate[1]-1)/ μ - ceil(ir * (currstate[1]-1)/ μ)) + 1
    env.state .= currstate
    getobs!(env.obs, env)
    env.time[1] += timestep(env)
    env
end



timestep(env::MuNearlyPeriod) = 1
Base.time(env::MuNearlyPeriod) = env.time[1]