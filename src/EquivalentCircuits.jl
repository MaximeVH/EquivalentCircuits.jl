module EquivalentCircuits

using Random, Distributions, SyntaxTree, Optim, Plots, ExprRules


include("CircuitEvolution.jl")
include("CircuitFunction.jl")
include("CircuitObjective.jl")
include("EvolutionOperators.jl")
include("InitializeParameters.jl")
include("ObjectiveFunction.jl")
include("OptimizeParameters.jl")
include("RC_Circuit.jl")
include("SimulateImpedance.jl")

export InitializeCircuit, CircuitObjective, OptimizeParameters
export InitializePopulation, CircuitEvolution


end
