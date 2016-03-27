"""
Module for phylogenetics-related MCMC.
"""
module MCMC

import Base.show, Base.log, Base.close

#export TargetDistribution, getLogDensity,
#    Operator, propose,
#    State, store, restore, getLogValue, getLogName,
#    Logger, log,
#    run


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

getLogValue(state::State) = state.value
getLogName(state::State) = state.name

abstract Logger
function log(logger::Logger, iter::Integer) end
function close(logger::Logger) end


function run{S<:State,O<:Operator,L<:Logger}(d::TargetDistribution,
    initialStateArray::Array{S,1}, operators::Array{O,1}, loggers::Array{L,1},
    nIters::Integer)

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


type ScaleOperator <: Operator
    scaleFactor::Float64
    param::State{Float64}
end

function propose(op::ScaleOperator)
    minF = min(op.scaleFactor, 1/op.scaleFactor)

    f = minF + rand()*(1/minF - minF)

    op.param.value *= f

    return -log(f)
end

type UniformOperator <: Operator
    windowSize::Float64
    param::State{Float64}
end

function propose(op::UniformOperator)
    op.param.value += (rand()-0.5)*op.windowSize

    return 0
end



type GaussianDistribution <: TargetDistribution
    mean::Float64
    variance::Float64
    x::State{Float64}
end

function getLogDensity(d::GaussianDistribution)
    return -(d.x.value - d.mean)^2/(2*d.variance)
end



type FlatTextLogger <: Logger
    outStream::IOStream
    states::Array{State,1}
    samplePeriod::Integer
end

function FlatTextLogger{T}(fileName::ASCIIString, states::Array{State{T},1}, samplePeriod::Integer)
    outStream = open(fileName, "w")

    print(outStream, "Sample")

    for state in states
        print(outStream, string("\t", getLogName(state)))
    end

    print(outStream, "\n")

    FlatTextLogger(outStream, states, samplePeriod)
end

function log(logger::FlatTextLogger, iter::Integer)

    if iter % logger.samplePeriod != 0
        return
    end

    print(logger.outStream, iter)
    for state in logger.states
        print(logger.outStream, string("\t", getLogValue(state)))
    end

    print(logger.outStream, "\n")
end

function close(logger::FlatTextLogger)
    close(logger.outStream)
end


function main()
    x = State{Float64}("x", 1.0)
    #op = ScaleOperator(1.2, x)
    op = UniformOperator(0.5, x)
    d = GaussianDistribution(10, 0.1, x)
    logger = FlatTextLogger("samples.log", [x], 10)

    run(d, [x], [op], [logger], 100000)
end

end
