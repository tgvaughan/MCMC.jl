"""
Module for phylogenetics-related MCMC.
"""
module MCMC

import Base.show, Base.log, Base.close
import StatsBase.sample, StatsBase.WeightVec

type UnimplementedMethodException <: Exception end

abstract TargetDistribution
function getLogDensity(d::TargetDistribution) end
getDeps(d::TargetDistribution) = throw(UnimplementedMethodException())

abstract Operator
getDeps(op::Operator) = throw(UnimplementedMethodException())
function propose(op::Operator) end

type State{T}
    name
    value::T
    storedValue::T
end

State{T}(name, value::T) = State{T}(name, value, value)

function store(state::State)
    state.storedValue = state.value
end

function restore(state::State)
    state.value, state.storedValue = state.storedValue, state.value
end

getLogNames(state::State) = [state.name]
getLogValues(state::State) = [state.value]

include("Floats.jl")

abstract Logger
function init(logger::Logger) end
function log(logger::Logger, iter::Int) end
function close(logger::Logger) end
summary(logger::Logger) = print("This logger does not provide a summary.")
trace(logger::Logger) = print("This logger does not provide a trace.")

include("Loggers.jl")


include("Trees.jl")

"""
Run MCMC chain.
"""
function run{O<:Operator,L<:Logger}(d::TargetDistribution,
    weightedOperators::Array{Tuple{O,Float64},1}, loggers::Array{L,1}, nIters::Int)

    # Initialize loggers
    for logger in loggers
        init(logger)
    end

    oldLogDensity = getLogDensity(d)
    weightVec = WeightVec(map(x->x[2], weightedOperators))
    opVec = map(x->x[1], weightedOperators)

    for iter in 1:nIters

       # Propose new state
       op = sample(opVec, weightVec)
       store(getDeps(op))

       HF = propose(op)

       if HF > -Inf
           # Evaluate target density of new state
           newLogDensity = getLogDensity(d)

           # Evaluate acceptance probability
           alpha = exp(newLogDensity - oldLogDensity + HF)

           # Determine fate of proposal
           if alpha > 1 || rand()<alpha
               oldLogDensity = newLogDensity
           else
               restore(getDeps(op))
           end
       else
           restore(getDeps(op))
       end

       # Log state if necessary
       for logger in loggers
           log(logger, iter)
       end
    end

    # Finalize loggers
    for logger in loggers
        close(logger)
    end
end

function run{O<:Operator,L<:Logger}(d::TargetDistribution,
            operators::Array{O,1}, loggers::Array{L,1}, nIters::Int)
    run(d, [op => 1.0 for op = operators], loggers, nIters)
end


function store{S<:State}(stateArray::Array{S,1})
    for state in stateArray
        store(state)
    end
end

function restore{S<:State}(stateArray::Array{S,1})
    for state in stateArray
        restore(state)
    end
end

end
