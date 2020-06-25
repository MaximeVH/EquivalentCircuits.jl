module EquivalentCircuits

using Random, Distributions, StatsBase, SyntaxTree, Optim, Plots, LinearAlgebra, ExprRules

#Optimal parameter fitting.

include("alt_expression_as_array.jl")
include("Circut_expression.jl")
include("Create_circuit_function.jl")
include("Get_parameter_dictionary.jl")
include("Get_parameter_elements.jl")
include("Get_parameter_values.jl")
include("new_objective.jl")
include("Optimal_circuit_parameters.jl")
include("Simulate_impedance_ex.jl")

export Optimal_circuit_parameters


end
