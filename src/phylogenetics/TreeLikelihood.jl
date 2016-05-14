type Alignment
    taxa::Array{ASCIIString,1}
    patterns::Array{Int,2}
    weights::Array{Int,1}
end

function show(io::IO, a::Alignment)
    print(io, "Alignment of ", length(a.taxa), " sequences containing ",
                sum(a.weights), " sites with ",
                length(a.weights), " unique site patterns.")
end


function patternLT(a::Array{Int,1}, b::Array{Int,1})
    if length(a) != length(b)
        error("Patterns have different lengths.")
    end

    for site in 1:length(a)
        if a[site] < b[site]
            return true
        elseif a[site] > b[site]
            return false
        end
    end

    return false
end

function patternEq(a::Array{Int,1}, b::Array{Int,1})
    if length(a) != length(b)
        error("Patterns have different lengths.")
    end

    for site in 1:length(a)
        if a[site] != b[site]
            return false
        end
    end

    return true
end

dnaMapping = Dict( 'A' => 1,
            'C' => 2,
            'G' => 3,
            'T' => 4,
            '-' => -1,
            '?' => -1)

function readFasta(filename; mapping=dnaMapping)
    f = open(filename)

    sequences = Dict{ASCIIString, Array{Int,1}}()

    seq = nothing
    header = nothing
    for line in eachline(f)
        line = chomp(line)

        if startswith(line, ">")
            if seq != nothing
                sequences[header] = map((c) -> mapping[uppercase(c)], collect(seq))
            end

            header = line[2:length(line)]
            seq = ""
        else
            seq = string(seq, line)
        end
    end

    close(f)

    # Ensure that sequences have the same length:

    nSites = nothing
    for (taxon, seq) in sequences
        if nSites == nothing
            nSites = length(seq)
        else
            if nSites != length(seq)
                error("Sequences do not have the same length.  Are they aligned?")
            end
        end
    end

    # Construct ordered taxon list

    taxa = collect(keys(sequences))
    nTaxa = length(taxa)

    # Collect site patterns

    patterns = Array{Array{Int,1},1}(nSites)
    for s in 1:nSites
        pattern = Array{Int,1}(nTaxa)
        for t in 1:nTaxa
            pattern[t] = sequences[taxa[t]][s]
        end
        patterns[s] = pattern
    end

    # Find unique patterns and counts

    sortedPat = sort(patterns, lt=patternLT)
    uniquePat = Array{Array{Int,1},1}()
    patternWeights = Array{Int,1}()

    push!(uniquePat, patterns[1])
    push!(patternWeights, 1)

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

    # Assemble and return alignment

    return Alignment(taxa, patternMatrix, patternWeights)
end
