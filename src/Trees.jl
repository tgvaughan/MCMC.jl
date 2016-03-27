using TimeTrees

export CoalescentTree, CoalescentDistribution, TreeScaler, TreeUniform

# State initialization

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

function getLogDensity(d::CoalescentDistribution)
    tree = d.treeState.value
    logP = 0.0

    nodes = getNodes(tree)
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

type TreeScaler <: Operator
    scaleFactor::Float64
    treeState::State{TimeTree}
end

function propose(op::TreeScaler)
    fmin = min(op.scaleFactor, 1/op.scaleFactor)
    f = fmin + (1/fmin - fmin)*rand()

    tree = op.treeState.value

    for node in getInternalNodes(tree)
        node.age *= f;

        for child in node.children
            if isLeaf(child) && child.age > node.age
                return -Inf
            end
        end
    end

    return (length(getInternalNodes(tree))-2)*log(f)
end

type TreeUniform <: Operator
    treeState::State{TimeTree}
end

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


# Loggers

getScreenLogName{T<:TimeTree}(state::State{T}) = string(state.name, "_height")
getScreenLogValue{T<:TimeTree}(state::State{T}) = state.value.root.age

getFlatTextLogName{T<:TimeTree}(state::State{T}) = string(state.name, "_height")
getFlatTextLogValue{T<:TimeTree}(state::State{T}) = state.value.root.age
