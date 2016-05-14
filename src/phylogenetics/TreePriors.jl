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

