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


# Testing

function testCoalescent(;sim=true)
    taxa = [string(i) => 0.0 for i = 1:1000]

    t = State("tree", CoalescentTree(taxa, 1.0))
    ops = [(TreeScaler(0.5, t), 1.0),
        (TreeRootScaler(0.2, t), 1.0),
        (TreeUniform(t), 1.0),
        (TreeWilsonBalding(0.1, t), 1.0),
        (SubtreeSlide(0.1, t), 1.0),
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
        nSims = 10000
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
