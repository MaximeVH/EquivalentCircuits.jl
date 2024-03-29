module EquivalentCircuits

    export circuit_evolution, circuit_search, circuit_evolution_batch
    export parameteroptimisation, circuitfunction
    export loadpopulation, get_parameter_upper_bound
    export simulateimpedance_noiseless
    export Circuit, EquivalentCircuit
    using Random, Combinatorics, GeneralizedGenerated, DelimitedFiles, Distributions, Optim
    using BlackBoxOptim
    using Hwloc
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
    include("Miscellaneous.jl")
    include("CircuitEvolution.jl")
    # include("MAP_elites.jl")

end
