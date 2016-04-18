# Distributions and operators relevant to phylogenetic trees

using TimeTrees

# State representation

State{T<:TimeTree}(name, tree::T) = State(name, tree, getCopy(tree))

# Tree state storing (only topology and edge lengths are stored)

function store{T<:TimeTree}(state::State{T})
    tree = state.value
    storedTree = state.storedValue

    for i in 1:length(tree.nodes)
        node = tree.nodes[i]
        storedTree.nodes[i].age = node.age

        if isRoot(node)
            storedTree.nodes[i].parent = storedTree.nodes[i]
            storedTree.root = storedTree.nodes[i]
        else
            storedTree.nodes[i].parent = storedTree.nodes[node.parent.number]
        end

        for ci in 1:length(node.children)
            storedTree.nodes[i].children[ci] = storedTree.nodes[node.children[ci].number]
        end

    end
        
end

function restore{T<:TimeTree}(state::State{T})
    state.value, state.storedValue = state.storedValue, state.value
end

# State initialization

## Coalescent simulation

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
        coalRate = 0.5*k*(k-1)/popSize

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

## Coalescent distribution

type CoalescentDistribution <: TargetDistribution
    popSize::Float64
    treeState::State{TimeTree}
end
getDeps(d::CoalescentDistribution) = [d.treeState]

function getLogDensity(d::CoalescentDistribution)
    tree = d.treeState.value
    logP = 0.0

    nodes = getNodes(tree)
    ages = zeros(length(nodes))
    for i in 1:length(nodes)
        ages[i] = nodes[i].age
    end

    nodeIndices = sortperm(ages)

    k = 1
    t = 0.0
    for i in nodeIndices[2:length(nodes)]
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

## Operator utility functions

function getSibling(node::Node)
    node.parent.children[1] == node ? node.parent.children[2] : node.parent.children[1]
end

function replaceChild(parent::Node, oldChild::Node, newChild::Node)
    newChild.parent = parent
    for i in 1:length(parent.children)
        if parent.children[i] == oldChild
            parent.children[i] = newChild
            break
        end
    end
end

## Tree Scalers

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

type TreeRootScaler{T<:TimeTree} <: Operator
    scaleFactor::Float64
    treeState::State{T}
end
getDeps(op::TreeRootScaler) = [op.treeState]

function propose(op::TreeRootScaler)
    fmin = min(op.scaleFactor, 1/op.scaleFactor)
    f = fmin + (1/fmin - fmin)*rand()

    tree = op.treeState.value

    tree.root.age *= f;
    for child in tree.root.children
        if tree.root.age < child.age
            return -Inf
        end
    end

    return -log(f)
end

## Uniform node height shifter

type TreeUniform{T<:TimeTree} <: Operator
    treeState::State{T}
end
getDeps(op::TreeUniform) = [op.treeState]

function propose(op::TreeUniform)
    tree = op.treeState.value

    if getLeafCount(tree)<3
        return -Inf
    end

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

## Wilson-Balding move

type TreeWilsonBalding{T<:TimeTree} <: Operator
    alpha::Float64
    treeState::State{T}
end
getDeps(op::TreeWilsonBalding) = [op.treeState]

function invalidSrcNode(srcNode::Node)
    if isRoot(srcNode)
        return true
    end

    sister = getSibling(srcNode)

    if isRoot(srcNode.parent)
        return isLeaf(sister) || srcNode.age >= sister.age
    end

    return false
end

function invalidDestNode(srcNode::Node, destNode::Node)
    if destNode == srcNode || destNode == srcNode.parent || destNode.parent == srcNode.parent
        return true
    end

    if !isRoot(destNode) && (destNode.parent.age <= srcNode.age)
        return true
    end

    return false
end

function propose(op::TreeWilsonBalding)
    tree = op.treeState.value

    if getLeafCount(tree)<3
        return -Inf
    end

    nodes = getNodes(tree)

    srcNode = rand(nodes)
    while invalidSrcNode(srcNode)
        srcNode = rand(nodes)
    end

    destNode = rand(nodes)
    while invalidDestNode(srcNode, destNode)
        destNode = rand(nodes)
    end

    if isRoot(destNode)
        # Forward root move

        logHR = 0.0

        srcNodeS = getSibling(srcNode)
        srcNodeG = srcNode.parent.parent

        logHR += -log(srcNodeG.age - max(srcNode.age, srcNodeS.age))

        replaceChild(srcNodeG, srcNode.parent, srcNodeS)
        replaceChild(srcNode.parent, srcNodeS, destNode)

        tree.root = srcNode.parent
        srcNode.parent.parent = srcNode.parent

        srcNode.parent.age = destNode.age + randexp()*op.alpha

        logHR -= -(srcNode.parent.age-destNode.age)/op.alpha + log(1/op.alpha)

        return logHR
    end

    if isRoot(srcNode.parent)
        # Reverse root move

        logHR = 0.0

        srcNodeS = getSibling(srcNode)

        logHR += -(srcNode.parent.age - srcNodeS.age)/op.alpha + log(1/op.alpha)

        destNodeP = destNode.parent

        replaceChild(destNodeP, destNode, srcNode.parent)
        replaceChild(srcNode.parent, srcNodeS, destNode)

        tree.root = srcNodeS
        srcNodeS.parent = srcNodeS

        minAge = max(srcNode.age, destNode.age)
        maxAge = destNodeP.age
        srcNode.parent.age = minAge + rand()*(maxAge - minAge)

        logHR -= -log(maxAge-minAge)

        return logHR
    end

    # Non-root move

    logHR = 0.0

    srcNodeS = getSibling(srcNode)
    srcNodeG = srcNode.parent.parent

    logHR += -log(srcNodeG.age-max(srcNodeS.age,srcNode.age))

    destNodeP = destNode.parent

    replaceChild(srcNodeG, srcNode.parent, srcNodeS)
    replaceChild(destNodeP, destNode, srcNode.parent)
    replaceChild(srcNode.parent, srcNodeS, destNode)

    minAge = max(destNode.age, srcNode.age)
    maxAge = destNodeP.age
    srcNode.parent.age = minAge + rand()*(maxAge - minAge)

    logHR -= -log(maxAge-minAge)

    return logHR

end

## Subtree exchange move

type SubtreeExchange{T<:TimeTree} <: Operator
    alpha::Float64
    isWide::Bool
    treeState::State{T}
end
getDeps(op::SubtreeExchange) = [op.treeState]

SubtreeExchangeNarrow{T<:TimeTree}(alpha::Float64, treeState::State{T}) =
    SubtreeExchange{T}(alpha, false, treeState)

SubtreeExchangeWide{T<:TimeTree}(alpha::Float64, treeState::State{T}) =
    SubtreeExchange{T}(alpha, true, treeState)

function propose(op::SubtreeExchange)
    tree = op.treeState.value
    nodes = getNodes(tree)

    if !op.isWide
        srcNode = rand(nodes)
        while isRoot(srcNode) || isRoot(srcNode.parent)
            srcNode = rand(nodes)
        end

        destNode = getSibling(srcNode.parent)
    else
        srcNode = rand(nodes)
        while isRoot(srcNode)
            srcNode = rand(nodes)
        end

        destNode = rand(nodes)
        while destNode == srcNode ||
                isRoot(destNode) ||
                destNode.parent == srcNode.parent ||
                destNode == srcNode.parent ||
                srcNode == destNode.parent

            destNode = rand(nodes)
        end

    end

    if (srcNode.parent.age < destNode.age) || (destNode.parent.age < srcNode.age)
        return -Inf
    end

    srcNodeP = srcNode.parent
    destNodeP = destNode.parent
    replaceChild(srcNodeP, srcNode, destNode)
    replaceChild(destNodeP, destNode, srcNode)

    return 0.0
end

## Subtree slide move

type SubtreeSlide{T<:TimeTree} <: Operator
    size::Float64
    treeState::State{T}
end
getDeps(op::SubtreeSlide) = [op.treeState]

function propose(op::SubtreeSlide)
    tree = op.treeState.value
    nodes = getNodes(tree)

    node = rand(nodes)
    while isRoot(node)
        node = rand(nodes)
    end

    delta = randn()*op.size
    
end

# Loggers

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

# Testing

function testCoalescent(;sim=true)
    taxa = [string(i) => 0.0 for i = 1:1000]

    t = State("tree", CoalescentTree(taxa, 1.0))
    ops = [(TreeScaler(0.5, t), 1.0),
        (TreeRootScaler(0.2, t), 1.0),
        (TreeUniform(t), 1.0),
        (TreeWilsonBalding(0.1, t), 1.0),
        (SubtreeExchangeNarrow(0.1, t), 1.0),
        (SubtreeExchangeWide(0.1, t), 1.0)]

    d = CoalescentDistribution(1.0, t)

    loggers = [ScreenLogger([t], 10000),
                FlatTextLogger("samples.log", [t], 100),
                TreeLogger("trees.log", t, 1000)]

    run(d, ops, loggers, 1000000)

    # Simulation for comparison
    if sim
        print("\nSimulating coalescent trees for comparison")

        outStream = open("sims.log", "w")
        println(outStream, "Sample\ttree_height\ttree_length")
        nSims = 1000
        for i in 1:nSims
            if i % (nSims/10) == 0
                print(".")
            end

            tree = CoalescentTree(taxa, 1.0)
            println(outStream, "$(i-1)\t$(tree.root.age)\t$(getTreeLength(tree))")
            flush(outStream)
        end
        close(outStream)
        println("done")
    end
end

