# Distributions and operators relevant to real numbers

# Operators

type ScaleOperator <: Operator
    scaleFactor::Float64
    param::State{Float64}
end
getDeps(op::ScaleOperator) = [op.param]

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
getDeps(op::UniformOperator) = [op.param]

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
getDeps(d::GaussianDistribution) = [d.x]

getLogDensity(d::GaussianDistribution) = -(d.x.value - d.mean)^2/(2*d.variance)

type ExponentialDistribution <: TargetDistribution
    mean::Float64
    x::State{Float64}
end
getDeps(d::ExponentialDistribution) = [d.x]

getLogDensity(d::ExponentialDistribution) = -d.x.value/d.mean - log(d.mean)


# Testing
function testGaussian()
    x = State("x", 1.0)
    op = UniformOperator(0.5, x)
    d = GaussianDistribution(10, 0.1, x)
    fileLogger = FlatTextLogger("samples.log", [x], 100)
    screenLogger = ScreenLogger([x], 100000)

    run(d, [op], [fileLogger, screenLogger], 1000000)
end


