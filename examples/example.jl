using CSV
using DataFrames
using EquivalentCircuits

"""
Circuit identification is done with the `circuit_evolution` function. This function can take
several keyword arguments, allowing users to tune the circuit identification if necessary.
These can be found in the function documentation. The measurements and frequencies should
be provided either as a CSV file, or as arrays.
"""
# Example 1: Input is CSV file with as columns: Re(Z), Im(Z), frequencies measurements
circuit = circuit_evolution("example_measurements.csv")

# Example 2: Input is an array with as columns: Re(Z), Im(Z), frequencies measurements
data = CSV.read("example_measurements.csv", DataFrame)
Z = data[:, 1] .+ data[:, 2]im
freq = data[:, 3]
@time circuit = circuit_evolution(Z, freq)

"""
When the equivalent circuit to be used is known, its parameters can be fit using the
`parameteroptimisation` function.
"""
p = parameteroptimisation("[[C1,R2]-R3,P4]", Z, freq)

"""
If you want to generate multiple candidate circuits, you can use the `circuit_evolution_batch`
function. This function is parallelised on all available cores.
"""
circuits = circuit_evolution_batch(Z, freq; iters=6);
