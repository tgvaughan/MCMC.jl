# Distributions and operators relevant to real numbers

# Operators

type ScaleOperator <: Operator
    scaleFactor::Float64
    param::State{Float64}
end

function propose(op::ScaleOperator)
    minF = min(op.scaleFactor, 1/op.scaleFactor)

    f = minF + rand()*(1/minF - minF)

    op.param.value *= f

    return -log(f)
end

type UniformOperator <: Operator
    windowSize::Float64
    param::State{Float64}
end

function propose(op::UniformOperator)
    op.param.value += (rand()-0.5)*op.windowSize

    return 0
end


# Basic parametric distributions

type GaussianDistribution <: TargetDistribution
    mean::Float64
    variance::Float64
    x::State{Float64}
end

function getLogDensity(d::GaussianDistribution)
    return -(d.x.value - d.mean)^2/(2*d.variance)
end



