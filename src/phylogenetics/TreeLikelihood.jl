# Sequence alignment I/O

abstract SequenceDataType
getStateCount(dt::SequenceDataType) = throw(UnimplementedMethodException())
stateFromChar(dt::SequenceDataType, c::Char) = throw(InimplementedMethodException())

type DNA <: SequenceDataType end
getStateCount(dt::DNA) = 4
stateFromChar(dt::DNA, c::Char) = Dict( 'A' => 1,
                                        'C' => 2,
                                        'G' => 3,
                                        'T' => 4,
                                        '-' => -1,
                                        '?' => -1)[c]

type Alignment
    datatype::SequenceDataType

    taxa::Array{ASCIIString,1}
    patterns::Array{Int,2}
    weights::Array{Int,1}

    variablePatterns::Array{Int,1}
    constantPatterns::Array{Int,1}

    ambigiousPatterns::Set{Int}
end

getPatternCount(a::Alignment) = length(a.weights)
getSiteCount(a::Alignment) = sum(a.weights)

function show(io::IO, a::Alignment)
    print(io, "Alignment of ", length(a.taxa), " sequences containing ",
                getSiteCount(a), " sites with ",
                getPatternCount(a), " unique site patterns.")
end


function Alignment(sequences::Dict{ASCIIString, ASCIIString}; datatype::SequenceDataType = DNA())
    # Map sequence strings to integer arrays
    mappedSeqs = [taxon => map((c) -> stateFromChar(datatype, c), collect(uppercase(seq)))
                    for (taxon,seq) in sequences]

    # Ensure that sequences have the same length:

    nSites = nothing
    for (taxon, seq) in mappedSeqs
        if nSites == nothing
            nSites = length(seq)
        else
            if nSites != length(seq)
                error("Sequences do not have the same length.  Are they aligned?")
            end
        end
    end

    # Construct ordered taxon list

    taxa = collect(keys(mappedSeqs))
    nTaxa = length(taxa)

    # Collect site patterns

    patterns = Vector{Vector{Int}}(nSites)
    for s in 1:nSites
        pattern = Vector{Int}(nTaxa)
        for t in 1:nTaxa
            pattern[t] = mappedSeqs[taxa[t]][s]
        end
        patterns[s] = pattern
    end


    # Find unique patterns and counts

    sortedPat = sort(patterns, lt=(a,b) -> begin
        for site in 1:length(a)
            if a[site] < b[site]
                return true
            elseif a[site] > b[site]
                return false
            end
        end
    return false

    end)
    uniquePat = Array{Array{Int,1},1}()
    patternWeights = Array{Int,1}()

    push!(uniquePat, sortedPat[1])
    push!(patternWeights, 1)

    function patternEq(a::Array{Int,1}, b::Array{Int,1})
        for site in 1:length(a)
            if a[site] != b[site]
                return false
            end
        end

        return true
    end

    for s in 2:nSites
        if patternEq(sortedPat[s], sortedPat[s-1])
            patternWeights[length(patternWeights)] += 1
        else
            push!(uniquePat, sortedPat[s])
            push!(patternWeights, 1)
        end
    end

    patternMatrix = Array{Int,2}(length(uniquePat), nTaxa)
    for p in 1:length(uniquePat)
        for t in 1:nTaxa
            patternMatrix[p,t] = uniquePat[p][t]
        end
    end

    # Record indices of constant and variable patterns

    variablePatterns = Array{Int,1}()
    constantPatterns = Array{Int,1}()
    ambiguousPatterns = Set{Int}()
    for p in 1:size(patternMatrix)[1]
        isConstant = true
        isAmbiguous = false
        for t in 1:nTaxa
            if t>1 && isConstant && (patternMatrix[p,t] != patternMatrix[p,t-1])
                isConstant = false
            end

            if !isAmbiguous && patternMatrix[p,t]<0
                isAmbiguous = true
            end
        end

        if isConstant
            push!(constantPatterns, p)
        else
            push!(variablePatterns, p)
        end

        if isAmbiguous
            push!(ambiguousPatterns, p)
        end
    end

    # Assemble and return alignment

    return Alignment(datatype, taxa, patternMatrix, patternWeights,
        variablePatterns, constantPatterns, ambiguousPatterns)
end


function readFasta(filename; datatype::SequenceDataType = DNA())
    f = open(filename)

    sequences = Dict{ASCIIString, ASCIIString}()

    seq = nothing
    header = nothing
    for line in eachline(f)
        line = chomp(line)

        if startswith(line, ">")
            if seq != nothing
                sequences[header] = seq
            end

            header = line[2:length(line)]
            seq = ""
        else
            seq = string(seq, line)
        end
    end

    if seq != nothing
        sequences[header] = seq
    end

    close(f)

    return Alignment(sequences, datatype=datatype)
end

# Likelihood

## Substitution models

abstract SubstitutionModel{D}
updateTransitionMatrix(model::SubstitutionModel, matrix::Array{Float64,2}, dist::Float64) =
    throw(UnimplementedMethodException())

type JukesCantor <: SubstitutionModel{DNA} end
function updateTransitionMatrix(jc::JukesCantor, matrix::Array{Float64,2}, dist::Float64)
    pSame = (1 + 3*exp(-4/3*dist))/4
    pDiff = 1 - pSame

    fill!(matrix, pDiff)
    for i in 1:4
        matrix[i,i] = pSame
    end
end

type TreeLikelihood <: TargetDistribution
    alignment::Alignment
    clockRate::Float64
    substModel::SubstitutionModel
    treeState::State{TimeTree}
end
getDeps(d::TreeLikelihood) = [d.treeState]

function getLogDensity(d::TreeLikelihood)

    tree = d.treeState.value
    alignment = d.alignment
    substModel = d.substModel

    nodes = getNodes(tree)
    nNodes = length(nodes)

    nStates = getStateCount(alignment.datatype)

    partials = Array{Vector{Float64}}(nNodes,
                                getPatternCount(alignment))
    for n in 1:nNodes
        for pIdx in 1:getPatternCount(alignment)
            partials[n,pIdx] = zeros(nStates)
        end
    end
    partialsDirty = Vector{Bool}(nNodes)
    fill!(partialsDirty, true)

    transProbs = Vector{Matrix{Float64}}(nNodes)
    for n in 1:(nNodes)
        transProbs[n] = zeros(nStates,nStates)
    end
    transProbsDirty = Vector{Bool}(nNodes)
    fill!(transProbsDirty, true)
    
    function computePartials(nodeNr)

        if !partialsDirty[nodeNr]
            return 
        end

        node = nodes[nodeNr]

        left = node.children[1]
        leftNr = left.number

        right = node.children[2]
        rightNr = right.number

        if transProbsDirty[leftNr]
            updateTransitionMatrix(substModel, transProbs[leftNr], (node.age-left.age)/d.clockRate)
            transProbsDirty[leftNr] = false
        end

        if transProbsDirty[rightNr]
            updateTransitionMatrix(substModel, transProbs[rightNr], (node.age-right.age)/d.clockRate)
            transProbsDirty[rightNr] = false
        end

        if isLeaf(left) && isLeaf(right)
            for pIdx in alignment.variablePatterns
                partials[nodeNr, pIdx] =
                    transProbs[leftNr][:,alignment.patterns[leftNr, pIdx]] .*
                    transProbs[rightNr][:,alignment.patterns[pIdx,rightNr]]
            end
        elseif isLeaf(left)
            computePartials(rightNr)

            for pIdx in alignment.variablePatterns
                partials[pIdx, nodeNr] =
                    transProbs[leftNr][:,alignment.patterns[leftNr, pIdx]] .*
                    (transProbs[rightNr]*partials[pIdx, rightNr])
            end

        elseif isLeaf(right)
            computePartials(leftNr)

            for pIdx in alignment.variablePatterns
                partials[pIdx, nodeNr] =
                    (transProbs[leftNr]*partials[leftNr, pIdx]) .*
                    transProbs[rightNr][:,alignment.patterns[pIdx,rightNr]]
            end
        else
            computePartials(leftNr)
            computePartials(rightNr)

            for pIdx in alignment.variablePatterns
                partials[pIdx, nodeNr] =
                    (transProbs[leftNr]*partials[leftNr, pIdx]) .*
                    (transProbs[rightNr]*partials[pIdx, rightNr])
            end
        end
    end

    computePartials(tree.root.number)

    print(partials)
end

# Testing

function testLikelihood()
    t = State("tree", TimeTree("((A:1,B:1):0.5,C:1.5):0.0;"))

    of = open("align.fna", "w")
    println(of, ">A")
    println(of, "GTCA")
    println(of, ">B")
    println(of, "GCCA")
    println(of, ">C")
    println(of, "GTCA")
    close(of)

    msa = readFasta("align.fna")

    likelihood = TreeLikelihood(msa, 1.0, JukesCantor(), t)

    getLogDensity(likelihood)
end
