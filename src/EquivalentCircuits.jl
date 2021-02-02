module EquivalentCircuits
    export circuitevolution
    export parameteroptimisation
    using Random, Combinatorics, FunctionWrappers, DelimitedFiles, Distributions, Optim, Plots
    import Base: isless, length

    include("Circuits.jl") 
    include("CircuitFunction.jl")
    include("EvolutionOperators.jl")
    include("ObjectiveFunction.jl")
    include("OptimizeParameters.jl")
    include("SimulateImpedance.jl")
    include("CircuitSimplification.jl")
    include("CircuitLibrary.jl")
    include("RedundancyTesting.jl")
    include("CircuitEvolution.jl")
end
