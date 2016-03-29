using TimeTrees

# State initialization

State{T<:TimeTree}(name, tree::T) = State(name, tree, getCopy(tree))

# Tree state storing (only topology and edge lengths are stored)

function store{T<:TimeTree}(state::State{T})
    tree = state.value
    storedTree = state.storedValue

     for (i, node) in enumerate(tree.nodes)
        storedTree.nodes[i].age = node.age

        if isRoot(node)
            storedTree.nodes[i].parent = storedTree.nodes[i]
            storedTree.root = storedTree.nodes[i]
        else
            storedTree.nodes[i].parent = storedTree.nodes[node.parent.number]
        end

        for (ci, child) in enumerate(node.children)
            storedTree.nodes[i].children[ci] = storedTree.nodes[child.number]
        end

    end
        
end

function restore{T<:TimeTree}(state::State{T})
    state.value, state.storedValue = state.storedValue, state.value
end

# Simulate coalescent tree
function CoalescentTree(leafAges::Dict{ASCIIString,Float64}, popSize::Float64)

    # Create leaf nodes
    unprocessedLeaves = Array{Node,1}()
    for (taxon, height) in leafAges
        leaf = Node()
        leaf.label = taxon
        leaf.age = height
        push!(unprocessedLeaves, leaf)
    end

    # Sort leaves in order of increasing age 
    sort!(unprocessedLeaves, by=leaf -> leaf.age)

    # Create empty list of active nodes
    activeNodes = Array{Node,1}()

    t = 0.0
    while length(activeNodes)>1 || !isempty(unprocessedLeaves)

        k = length(activeNodes)
        coalRate = k*(k-1)/popSize

        # Sample time of next event
        if coalRate>0.0
            t += randexp()/coalRate
        else
            t = Inf
        end

        if !isempty(unprocessedLeaves) && t > unprocessedLeaves[1].age
            t = unprocessedLeaves[1].age
            push!(activeNodes, unprocessedLeaves[1])
            splice!(unprocessedLeaves, 1)
            continue
        end

        # Choose nodes to coalesce
        left = rand(activeNodes)
        right = rand(activeNodes)
        while left == right
            right = rand(activeNodes)
        end

        parent = Node()
        parent.age = t
        addChild!(parent, left)
        addChild!(parent, right)
        splice!(activeNodes, findfirst(n->n==left, activeNodes))
        splice!(activeNodes, findfirst(n->n==right, activeNodes))
        push!(activeNodes, parent)
    end

    return TimeTree(first(activeNodes))
end

# Distributions

type CoalescentDistribution <: TargetDistribution
    popSize::Float64
    treeState::State{TimeTree}
end
getDeps(d::CoalescentDistribution) = [d.treeState]

function getLogDensity(d::CoalescentDistribution)
    tree = d.treeState.value
    logP = 0.0

    nodes = getNodes(tree)[:]
    sort!(nodes, by=n->n.age)

    k = 1
    t = 0.0
    for i in 2:length(nodes)
        node = nodes[i]

        dt = node.age - t
        if dt>0.0
            logP += -dt*0.5*k*(k-1)/d.popSize
        end
        t = node.age

        if isLeaf(node)
            k += 1
        else
            logP += log(1/d.popSize)
            k -= 1
        end
    end

    return logP
end



# Operators

type TreeScaler{T<:TimeTree} <: Operator
    scaleFactor::Float64
    treeState::State{T}
end
getDeps(op::TreeScaler) = [op.treeState]

function propose(op::TreeScaler)
    fmin = min(op.scaleFactor, 1/op.scaleFactor)
    f = fmin + (1/fmin - fmin)*rand()

    tree = op.treeState.value

    for node in getInternalNodes(tree)
        node.age *= f;
    end

    for leaf in getLeaves(tree)
        if leaf.age > leaf.parent.age
            return -Inf
        end
    end

    return (length(getInternalNodes(tree))-2)*log(f)
end

type TreeUniform{T<:TimeTree} <: Operator
    treeState::State{T}
end
getDeps(op::TreeUniform) = [op.treeState]

function propose(op::TreeUniform)
    tree = op.treeState.value

    node = rand(getInternalNodes(tree))
    while isRoot(node)
        node = rand(getInternalNodes(tree))
    end

    maxAge = node.parent.age
    minAge = -Inf
    for child in node.children
        minAge = max(minAge, child.age)
    end

    node.age = minAge + (maxAge - minAge)*rand()

    return 0.0
end

type TreeWilsonBalding{T<:TimeTree} <: Operator
    alpha::Float64
    treeState::State{T}
end
getDeps(op::TreeWilsonBalding) = [op.treeState]

function propose(op::TreeWilsonBalding)
    tree = op.treeState.value

    if (length(getNodes(tree))<3)
        return -Inf
    end

    nodeI = rand(getNodes(tree))
    nodeJ = rand(getNodes(tree))
    while nodeJ == nodeI
        nodeJ = rand(getNodes(tree))
    end

    if !isRoot(nodeJ) && nodeJ.parent.age < nodeI.age || !isRoot(nodeI) 
    end

    # TODO
end

# Loggers

getScreenLogName{T<:TimeTree}(state::State{T}) = string(state.name, "_height")
getScreenLogValue{T<:TimeTree}(state::State{T}) = state.value.root.age

getFlatTextLogName{T<:TimeTree}(state::State{T}) = string(state.name, "_height")
getFlatTextLogValue{T<:TimeTree}(state::State{T}) = state.value.root.age

type TreeLogger{T<:TimeTree} <: Logger
    outStream::IOStream
    treeState::State{T}
    samplePeriod::Integer
end

function TreeLogger{T<:TimeTree}(fileName::ASCIIString, treeState::State{T}, samplePeriod::Integer)
    outStream = open(fileName, "w")
    TreeLogger(outStream, treeState, samplePeriod)
end

function init(logger::TreeLogger)
    println(logger.outStream, "#nexus")
    println(logger.outStream, "begin trees;")
end

function log(logger::TreeLogger, iter::Integer)
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

# Testing
function testCoalescent()
    t = State("tree", CoalescentTree([string(i) => 0.0 for i = 1:10], 1.0))
    ops = [TreeScaler(0.8, t)]

    d = CoalescentDistribution(1.0, t)

    loggers = [ScreenLogger([t], 100000),
                FlatTextLogger("samples.log", [t], 100),
                TreeLogger("trees.log", t, 1000)]

    run(d, ops, loggers, 1000000)
end

