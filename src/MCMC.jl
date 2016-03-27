"""
Module for phylogenetics-related MCMC.
"""
module MCMC

import Base.show, Base.log, Base.close

export State, run

abstract TargetDistribution
function getLogDensity(d::TargetDistribution) end

abstract Operator
function propose(op::Operator) end

type State{T}
    name
    value::T
    storedValue::T

    State(name, value) = new(name, value, value)
end

State{T}(name, value::T) = State{T}(name, value)

function store(state::State)
    state.storedValue = state.value
end

function restore(state::State)
    state.value = state.storedValue
end

getLogName(state::State) = state.name
getLogValue(state::State) = state.value

include("Floats.jl")

abstract Logger
function init(logger::Logger) end
function log(logger::Logger, iter::Integer) end
function close(logger::Logger) end
summary(logger::Logger) = print("This logger does not provide a summary.")
trace(logger::Logger) = print("This logger does not provide a trace.")

include("Loggers.jl")


include("Trees.jl")

"""
Run MCMC chain.
"""
function run{S<:State,O<:Operator,L<:Logger}(d::TargetDistribution,
    initialStateArray::Array{S,1}, operators::Array{O,1}, loggers::Array{L,1},
    nIters::Integer)

    # Initialize loggers
    for logger in loggers
        init(logger)
    end

    oldLogDensity = getLogDensity(d)
    stateArray = initialStateArray

    for iter in 1:nIters

       # Propose new state
       HF = propose(rand(operators))

       if HF > -Inf
           # Evaluate target density of new state
           newLogDensity = getLogDensity(d)

           # Evaluate acceptance probability
           alpha = exp(newLogDensity - oldLogDensity + HF)

           # Determine fate of proposal
           if alpha > 1 || rand()<alpha
               oldLogDensity = newLogDensity
               accept(stateArray)
           else
               reject(stateArray)
           end
       else
           reject(stateArray)
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

function accept{S<:State}(stateArray::Array{S,1})
    for state in stateArray
        store(state)
    end
end

function reject{S<:State}(stateArray::Array{S,1})
    for state in stateArray
        restore(state)
    end
end


# Testing
function main()
    x = State{Float64}("x", 1.0)
    #op = ScaleOperator(1.2, x)
    op = UniformOperator(0.5, x)
    d = GaussianDistribution(10, 0.1, x)
    fileLogger = FlatTextLogger("samples.log", [x], 100)
    screenLogger = ScreenLogger([x], 100000)

    run(d, [x], [op], [fileLogger, screenLogger], 1000000)
end

end
