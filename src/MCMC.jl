"""
Module for phylogenetics-related MCMC.
"""
module MCMC

import Base.show, Base.log, Base.close
import StatsBase.sample, StatsBase.WeightVec

type UnimplementedMethodException <: Exception end

# Abstract probability distribution type

abstract TargetDistribution
function getLogDensity(d::TargetDistribution) end
getDeps(d::TargetDistribution) = throw(UnimplementedMethodException())

# Abstract proposal operator type

abstract Operator
getDeps(op::Operator) = throw(UnimplementedMethodException())
function propose(op::Operator) end

# Parameterized state type with default store/restore

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

# Abstract logger type

abstract Logger
function init(logger::Logger) end
function log(logger::Logger, iter::Int) end
function close(logger::Logger) end
summary(logger::Logger) = print("This logger does not provide a summary.")
trace(logger::Logger) = print("This logger does not provide a trace.")

include("Floats.jl")
include("Loggers.jl")

include("phylogenetics/Phylogenetics.jl")

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
    
    opAcceptFreq = zeros(Int64, length(opVec))
    opChoiceFreq = zeros(Int64, length(opVec))

    for iter in 1:nIters

       # Choose operator
       opIdx = sample(1:length(opVec), weightVec)
       op = opVec[opIdx]
       opChoiceFreq[opIdx] += 1
       store(getDeps(op))

       # Propose new state
       HF = propose(op)

       if HF > -Inf
           # Evaluate target density of new state
           newLogDensity = getLogDensity(d)

           # Evaluate acceptance probability
           alpha = exp(newLogDensity - oldLogDensity + HF)

           # Determine fate of proposal
           if alpha > 1 || rand()<alpha
               oldLogDensity = newLogDensity
               opAcceptFreq[opIdx] += 1
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

    # Print operator summary
    println("\nOperator usage statistics")
    println("-------------------------")
    println("Operator Type\tSelected\tAccepted\tFraction")
    for i in 1:length(opVec)
        println(typeof(opVec[i]), "\t", opChoiceFreq[i], "\t", opAcceptFreq[i], "\t", opAcceptFreq[i]/opChoiceFreq[i])
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
