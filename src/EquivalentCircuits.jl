module EquivalentCircuits

    export circuit_evolution
    export parameteroptimisation
    export circuitfunction
    export loadpopulation
    export circuit_search
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
