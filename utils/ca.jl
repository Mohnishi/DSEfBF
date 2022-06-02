"""
cellular automata
"""

#------------------Packages-------------------#
using LinearAlgebra, Random, Statistics, Distributions


struct LifeGame{N, T} 
    width::Int
    height::Int
    noise::N
    state::Matrix{Int}
    obs::Matrix{T}
    init::Matrix{Int}
    time::Vector{Int}

    function LifeGame{T}(width, height, R; noise=nothing, init=nothing) where T <: AbstractFloat
        state = zeros(Int, height, width)
        obs = zeros(T, height, width)
        time = zeros(Int, 1)
        if init == nothing; init = period8init(width, height) end
        if noise == nothing; noise = Normal(0, R) end
        env = new{typeof(noise), T}(width, height, noise, state, obs,
                                init, time)
        reset!(env)
    end
end


function getstate!(state, env::LifeGame) 
    state .= env.state
end
function setstate!(env::LifeGame, state) 
    env.state .= state
end

function getobs!(obs, env::LifeGame)
    state = zeros(Int, size(env.state))
    getstate!(state, env) 
    obs .= state .+ rand(env.noise, size(state))
end

function getreward(index, env::LifeGame)
    rew = env.obs[index+4, 6]
    rew
end

function reset!(env::LifeGame)
    env.state .= copy(env.init)
    getobs!(env.obs, env)
    env.time[1] = 0
    env
end
function randreset!(env::LifeGame)
    env.state .= 0
    env.state[2:end-1,2:end-1] .= rand(0:1, env.height-2, env.width-2)
    getobs!(env.obs, env)
    env.time[1] = 0
    env
end


function step!(env::LifeGame)
    state, width, height = env.state, env.width, env.height
    currstate = copy(state)
    for i = 2:height-1
        for j = 2:width-1
            val = state[i,j]
            neighbor = sum(state[i-1:i+1, j-1:j+1]) - state[i,j]
            if (val == 1) && ((neighbor == 2) || (neighbor == 3))
                currstate[i,j] = 1
            elseif (val == 1) && ((neighbor == 1) || (neighbor > 3))
                currstate[i,j] = 0
            elseif (val == 0) && (neighbor == 3)
                currstate[i,j] = 1
            else
                currstate[i,j] = 0
            end
        end
    end
    env.state .= currstate
    getobs!(env.obs, env)
    env.time[1] += timestep(env)
    env
end

function period8init(width, height)
    state = zeros(Int, height, width)
    state[3,8] = state[4,7] = state[4,9] = state[5,6] = state[5,10] = 1
    state[6,5] = state[6,9] = state[7,4] = state[7,8] = state[8,3] = 1
    state[8,7] = state[9,4] = state[9,6] = state[10,5] = 1
    return state
end

timestep(env::LifeGame) = 1
Base.time(env::LifeGame) = env.time[1]