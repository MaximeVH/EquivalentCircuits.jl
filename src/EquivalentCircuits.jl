module EquivalentCircuits

using Random, Combinatorics, FunctionWrappers, Distributions, Optim, Plots
import Base: isless

include("Circuits.jl") 
include("CircuitFunction.jl")
include("EvolutionOperators.jl")
include("ObjectiveFunction.jl")
include("OptimizeParameters.jl")
include("SimulateImpedance.jl")
include("CircuitEvolution.jl")


export circuitevolution

end
