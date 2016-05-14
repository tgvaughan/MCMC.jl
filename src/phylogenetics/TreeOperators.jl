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

function getDescendents(node::Node, time::Float64)
    ancestors = [node]
    descendents = Array{Node,1}()

    done = false
    while !isempty(ancestors)
        ancestorsPrime = Array{Node,1}()
        for a in ancestors
            for c in a.children
                if c.age <=time
                    push!(descendents, c)
                else
                    push!(ancestorsPrime, c)
                end
            end
        end

        ancestors = ancestorsPrime
    end

    return descendents
end

function propose(op::SubtreeSlide)
    logHR = 0.0

    tree = op.treeState.value
    nodes = getNodes(tree)

    node = rand(nodes)
    while isRoot(node)
        node = rand(nodes)
    end
    parent = node.parent

    delta = randn()*op.size

    if delta > 0
        if !isRoot(parent) && parent.age + delta > parent.parent.age
            # Topology change

            sister = getSibling(node)
            grandParent = parent.parent

            destEdge = grandParent
            while !isRoot(destEdge) && destEdge.parent.age < parent.age + delta
                destEdge = destEdge.parent
            end

            # Identify number of descendents of destEdge at age of parent for HR
            logHR = -log(length(getDescendents(destEdge, parent.age)))

            if !isRoot(destEdge)
                destEdgeParent = destEdge.parent
                replaceChild(destEdgeParent, destEdge, parent)
            else
                tree.root = parent
                parent.parent = parent
            end
            replaceChild(parent, sister, destEdge)
            replaceChild(grandParent, parent, sister)
        end

    else
        if parent.age + delta < node.age
            return -Inf
        end

        sister = getSibling(node)
        if parent.age + delta < sister.age
            # Topology change

            # Identify potential destination edges
            potentials = getDescendents(sister, parent.age + delta)

            logHR = log(length(potentials))

            destEdge = rand(potentials)
            destEdgeParent = destEdge.parent

            if !isRoot(parent)
                grandParent = parent.parent
                replaceChild(grandParent, parent, sister)
            else
                tree.root = sister
                sister.parent = sister
            end
            replaceChild(destEdgeParent, destEdge, parent)
            replaceChild(parent, sister, destEdge)
        end
    end

    parent.age += delta

    return logHR
end


