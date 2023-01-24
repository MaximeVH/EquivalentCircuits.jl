module EquivalentCircuits

    export circuit_evolution, circuit_search
    export parameteroptimisation, circuitfunction
    export loadpopulation
    export simulateimpedance_noiseless
    export Circuit, EquivalentCircuit
    using Random, Combinatorics, GeneralizedGenerated, DelimitedFiles, Distributions, Optim
    using BlackBoxOptim
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
