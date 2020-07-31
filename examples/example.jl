#Material science "failed coating" circuit example.

using EquivalentCircuits

# Circuit notation R: resistor, C: capacitor, square brackets indicate a parallel subcircuit.

failed_coating = "R1-[C2,[C3,R4]-R5]"

#Parameters : an array in which the parameter values' order corresponds to the cicuit elements'
#oder of occurence in the circuit notation string.

failed_coating_parameters = [20,4E-9,4E-6,2500,3400]

#Range of experimental frequencies in which the data collection (or in this case simulation)
#takes place.

frequencies = [10.0^i for i in LinRange(-1, 5, 50)]

#Convert the circuit string notation into a callable function that can generate impedance
#values corresponding to the circuits configuration and input frequencies.

Failed_coating_function = CircuitFunction(failed_coating)

# Generate the simulated impedance measurements, in this case the default noise ratio value
# of 0.01 is used.

measurements = SimulateImpedance(Failed_coating_function,failed_coating_parameters,frequencies)

# Execute one generation of the genetic programming-based circuit configuration design algorithm.
New_population = CircuitEvolution(measurements,frequencies)

#Print fittest circuit's configuration after one generation:
println(New_population[1].string)
