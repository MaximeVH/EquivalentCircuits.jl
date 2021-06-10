using EquivalentCircuits, CSV

"""
Circuitidentification is done with the circuitevolution function. This function can take several keyword arguments,
allowing users to tune the circuitidentification if necessary. These can be found in the function documentation.
The measurements and frequencies should be provided either as a CSV file, or as arrays.

"""
#Measurements input as CSV file with as columns: real impedances, imaginary impedances and frequencies.
circuitevolution("example_measurements.csv")

#Example of input as arrays.
data = readdlm(raw"C:\Users\Dell\Documents\EquivalentCircuits.jl\example_measurements.csv",',')
measurements = data[:,1] .+ data[:,2]im
frequencies = data[:,3]
circuitevolution(measurements,frequencies)

"""
When the equivalent circuit to be used is known, its parameters can be fit using the parameteroptimisation function.

"""

parameteroptimisation("[[C1,R2]-R3,P4]", measurements, frequencies)