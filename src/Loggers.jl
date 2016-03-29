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

getFlatTextLogNames(state) = getLogNames(state)
getFlatTextLogValues(state) = getLogValues(state)

function init(logger::FlatTextLogger)
    print(logger.outStream, "Sample")

    for state in logger.states
        for name in getFlatTextLogNames(state)
            print(logger.outStream, string("\t", name))
        end
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
        for value in getFlatTextLogValues(state)
            print(logger.outStream, string("\t", value))
        end
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

getScreenLogNames(state) = getLogNames(state)
getScreenLogValues(state) = getLogValues(state)

function init(logger::ScreenLogger)
    print("Sample")

    for state in logger.states
        for name in getScreenLogNames(state)
            print(string("\t", name))
        end
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
        for value in getScreenLogValues(state)
            print(string("\t", value))
        end
    end

    speed = (time() - logger.startTime)/iter*1e6
    print("\t($(round(speed,3)) seconds/MSample)")

    println()

    flush(STDOUT)
end
