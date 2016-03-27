using TimeTrees

export CoalescentTree

# Simulate coalescent tree
function CoalescentTree(leafAges::Dict{ASCIIString,Float64}, popSize::Float64)

    # Create leaf nodes
    unprocessedLeaves = Array{Node,1}()
    for (taxon, height) in leafAges
        leaf = Node()
        leaf.label = taxon
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
