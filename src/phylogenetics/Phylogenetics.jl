# Distributions and operators relevant to phylogenetic trees

using TimeTrees

include("TreeStates.jl")
include("TreePriors.jl")
include("TreeOperators.jl")
include("TreeLoggers.jl")

# Operators

# Loggers

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

