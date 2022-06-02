"""
linear dyn
"""

#------------------Packages-------------------#
using LinearAlgebra, Random, Statistics, Distributions


struct LinearDyn{N, T} 
    dim::Int
    M::Matrix{T}
    noise::N
    state::Vector{T}
    obs::Vector{T}
    init::Vector{T}
    time::Vector{Int}

    function LinearDyn{T}(dim, M, R; noise=nothing, init=nothing) where T <: AbstractFloat
        state = zeros(T, dim)
        obs = zeros(T, dim)
        time = zeros(Int, 1)
        if init == nothing; init = lininit(dim) end
        if noise == nothing; noise = Normal(0, R) end
        env = new{typeof(noise), T}(dim, M, noise, state, obs,
                                init, time)
        reset!(env)
    end
end


function getstate!(state, env::LinearDyn) 
    state .= env.state
end
function setstate!(env::LinearDyn, state) 
    env.state .= state
end

function getobs!(obs, env::LinearDyn)
    state = zeros(size(env.state))
    getstate!(state, env) 
    obs .= state
end

function getreward(action::AbstractVector, env::LinearDyn)
    rew = env.state' * action + rand(env.noise)
    rew
end

function reset!(env::LinearDyn)
    env.state .= copy(env.init)
    getobs!(env.obs, env)
    env.time[1] = 0
    env
end
function randreset!(env::LinearDyn)
    env.state .= 0
    env.state .= randn(size(env.state))
    env.state ./ norm(env.state)
    getobs!(env.obs, env)
    env.time[1] = 0
    env
end


function step!(env::LinearDyn)
    state, M = env.state, env.M
    currstate = copy(state)
    env.state .= M * currstate
    getobs!(env.obs, env)
    env.time[1] += timestep(env)
    env
end

function lininit(dim) 
    state = randn(dim)
    state .= state ./ norm(state)
    return state
end

timestep(env::LinearDyn) = 1
Base.time(env::LinearDyn) = env.time[1]