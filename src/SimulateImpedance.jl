function simulateimpedance_noiseless(circuitfunc,parameters,frequencies) 
    return [circuitfunc(parameters,fr) for fr in frequencies] 
end

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