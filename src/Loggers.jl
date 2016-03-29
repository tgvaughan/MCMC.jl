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

getFlatTextLogName(state) = getLogName(state)
getFlatTextLogValue(state) = getLogValue(state)

function init(logger::FlatTextLogger)
    print(logger.outStream, "Sample")

    for state in logger.states
        print(logger.outStream, string("\t", getFlatTextLogName(state)))
    end

    println(logger.outStream)

    flush(logger.outStream)
end

function log(logger::FlatTextLogger, iter::Integer)

    if iter % logger.samplePeriod != 0
        return
    end

    print(logger.outStream, iter)
    for state in logger.states
        print(logger.outStream, string("\t", getFlatTextLogValue(state)))
    end

    println(logger.outStream)

    flush(logger.outStream)
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

getScreenLogName(state) = getLogName(state)
getScreenLogValue(state) = getLogValue(state)

function init(logger::ScreenLogger)
    print("Sample")

    for state in logger.states
        print(string("\t", getScreenLogName(state)))
    end

    logger.startTime = time()

    println()

    flush(STDOUT)
end

function log(logger::ScreenLogger, iter::Integer)

    if iter % logger.samplePeriod != 0
        return
    end

    print(iter)
    for state in logger.states
        print(string("\t", getScreenLogValue(state)))
    end

    speed = (time() - logger.startTime)/iter*1e6
    print("\t($(round(speed,3)) seconds/MSample)")

    println()

    flush(STDOUT)
end
