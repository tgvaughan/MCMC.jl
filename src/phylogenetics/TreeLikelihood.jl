# Sequence alignment I/O

abstract SequenceDataType
getNStates(dt::SequenceDataType) = throw(UnimplementedMethodException())
stateFromChar(dt::SequenceDataType, c::Char) = throw(InimplementedMethodException())

type DNA <: SequenceDataType end
getNStates(dt::DNA) = 4
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

    patterns = Array{Array{Int,1},1}(nSites)
    for s in 1:nSites
        pattern = Array{Int,1}(nTaxa)
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

    push!(uniquePat, patterns[1])
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
    pDiff = 1 - pStay

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

    partials = Array{Float64,4}(length(alignment.variablePatterns),
                                nNodes,
                                nStates)
    partialsDirty = Array{Bool,1}(nNodes)
    fill!(partialsDirty, true)

    transProbs = Array{Float64,3}(nNodes-1, nStates, nStates)
    transProbsDirty = Array{Bool,1}(nNodes-1)
    fill!(transProbsDirty, true)
    
    function computePartials(nodeNr, pIdx::Int, characterState::Int)

        if !partialsDirty[nodeNr]
            return partials[pIdx, nodeNr, characterState]
        end

        node = nodes[nodeNr]
        left = node.children[1]
        leftNr = left.number
        rightNr = right.number

        if transProbsDirty[leftNr]
            updateTransitionMatrix(substModel, sub(transProbs, leftNr,:,:), (node.age-left.age)/clockRate)
        end

        if transProbsDirty[rightNr]
            updateTransitionMatrix(substModel, sub(transProbs, rightNr,:,:), (node.age-right.age)/clockRate)
        end


        if isLeaf(node.children[1]) && isLeaf(node.children[2])
                
        elseif isLeaf(node.children[1])

        elseif isLeaf(node.children[2])

        else

        end

        if isLeaf(node)
            return characterState == alignment.patterns[pIdx, node.number] ? 1.0 : 0.0
        else
            for c in 1:nStates
                for cp in 1:nStates
                    for child in node.children
                    end
                end
            end
        end
    end
    for p in alignment.variablePatterns

    end
end
