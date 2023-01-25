
"""
simulateimpedance_noiseless(circuitfunc::ggfunc-function,parameters::Vector{Any},frequencies::Array{Float64,1}) 

Simulate the impedance spectrum of a given equivalent electrical circuit. 

The inputs are circuitfunc which is a function that calculates the impedance of the circuit for a given parameterset and frequency,
the parameters of the circuit, and the frequencies over which the impedance spectrum is to be simulated. 
The output is an array of simulated complex-valued impedance measurements.

"""
function simulateimpedance_noiseless(circuitfunc,parameters,frequencies) 
    return [circuitfunc(parameters,fr) for fr in frequencies] 
end


"""
simulateimpedance_noiseless(circuit::Circuit,frequencies::Array{Float64,1}) 

Simulate the impedance spectrum of a given equivalent electrical circuit, it's corresponding parameters,
and an array of frequencies. 

The inputs are a Circuit (circuit), which is a mutable structure containing the encoding of the circuit used by the evolutionary algorithm,
the parameters and the fitness value obtained during the evolutionary algorithm; and the frequencies over which the impedance spectrum is to be simulated. 
The output is an array of simulated complex-valued impedance measurements.

"""
function simulateimpedance_noiseless(circuit,frequencies) 
    circuitfunc = tree_to_function(get_circuit_tree(circuit))
    return [circuitfunc(circuit.parameters,fr) for fr in frequencies] 
end

function simulateimpedance(circuitfunc,parameters,frequencies,noise_ratio=0.01)
    impedances = [circuitfunc(parameters,fr) for fr in frequencies] 
    R = real(impedances)
    I = imag(impedances)
    Z = abs.(impedances)
    return [R[i]+rand(Normal(0,Z[i]*noise_ratio))+I[i]*im+rand(Normal(0,Z[i]*noise_ratio))im for i in 1:length(impedances)]
end

function simulateimpedance(circuit::Circuit,frequencies,noise_ratio=0.01)
    circuitfunc = tree_to_function(get_circuit_tree(circuit))
    impedances = [circuitfunc(circuit.parameters,fr) for fr in frequencies] 
    R = real(impedances)
    I = imag(impedances)
    Z = abs.(impedances)
    return [R[i]+rand(Normal(0,Z[i]*noise_ratio))+I[i]*im+rand(Normal(0,Z[i]*noise_ratio))im for i in 1:length(impedances)]
end