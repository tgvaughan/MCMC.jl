# Basic loggers

# Flat text file logger

type FlatTextLogger{S<:State} <: Logger
    outStream::IOStream
    states::Array{S,1}
    samplePeriod::Integer
end

function FlatTextLogger{S<:State}(fileName::ASCIIString, states::Array{S,1}, samplePeriod::Integer)
    outStream = open(fileName, "w")
    FlatTextLogger(outStream, states, samplePeriod)
end

function init(logger::FlatTextLogger)
    print(logger.outStream, "Sample")

    for state in logger.states
        print(logger.outStream, string("\t", getLogName(state)))
    end

    print(logger.outStream, "\n")
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


# Screen logger

type ScreenLogger{S<:State} <: Logger
    states::Array{S,1}
    samplePeriod::Integer
    startTime::Float64

    ScreenLogger(states, samplePeriod) = new(states, samplePeriod, 0.0)
end

ScreenLogger{S<:State}(states::Array{S,1}, samplePeriod) = ScreenLogger{S}(states, samplePeriod)

function init(logger::ScreenLogger)
    print("Sample")

    for state in logger.states
        print(string("\t", getLogName(state)))
    end

    logger.startTime = time()

    print("\n")
end

function log(logger::ScreenLogger, iter::Integer)

    if iter % logger.samplePeriod != 0
        return
    end

    print(iter)
    for state in logger.states
        print(string("\t", getLogValue(state)))
    end

    speed = (time() - logger.startTime)/iter*1e6
    print("\t($speed seconds/MSample)")

    print("\n")
end


