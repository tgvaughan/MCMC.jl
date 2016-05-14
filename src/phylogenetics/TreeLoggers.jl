## Summary statistics

function getTreeLength(tree::TimeTree)
    l = 0.0
    for node in getNodes(tree)
        if !isRoot(node)
            l += node.parent.age - node.age
        end
    end

    return l
end

## Screen and flat text logging

getScreenLogNames{T<:TimeTree}(state::State{T}) =
    [string(state.name, "_height"), string(state.name, "_length")]
getScreenLogValues{T<:TimeTree}(state::State{T}) = 
    [state.value.root.age, getTreeLength(state.value)]

getFlatTextLogNames{T<:TimeTree}(state::State{T}) = getScreenLogNames(state)
getFlatTextLogValues{T<:TimeTree}(state::State{T}) = getScreenLogValues(state)

## Tree logging

type TreeLogger{T<:TimeTree} <: Logger
    outStream::IOStream
    treeState::State{T}
    samplePeriod::Int
end

function TreeLogger{T<:TimeTree}(fileName::ASCIIString, treeState::State{T}, samplePeriod::Int)
    outStream = open(fileName, "w")
    TreeLogger(outStream, treeState, samplePeriod)
end

function init(logger::TreeLogger)
    println(logger.outStream, "#nexus")
    println(logger.outStream, "begin trees;")
end

function log(logger::TreeLogger, iter::Int)
    if iter % logger.samplePeriod != 0
        return
    end

    println(logger.outStream, "tree tree_$iter = $(getNewick(logger.treeState.value))")

    flush(logger.outStream)
end

function close(logger::TreeLogger)
    println(logger.outStream, "end;")
    close(logger.outStream)
end

