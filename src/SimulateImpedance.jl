function simulateimpedance_noiseless(circuitfunc,parameters,frequencies) 
    return [circuitfunc(parameters,fr) for fr in frequencies] 
end

function simulateimpedance_noiseless(circuit,frequencies) 
    circuitfunc = tree_to_function(get_circuit_tree(circuit))
    return [circuitfunc(circuit.parameters,fr) for fr in frequencies] 
end

function simulateimpedance(circuitfunc,parameters,frequencies,noise_ratio=0.01)
    impedances = [circuitfunc(parameters,fr) for fr in frequencies] 
    reals = real(impedances)
    imags = imag(impedances)
    return [R+rand(Normal(0,abs(R*noise_ratio)))+I*im + rand(Normal(0,abs(-I*noise_ratio)))im for (R,I) in zip(reals,imags)]
end

function simulateimpedance(circuit::Circuit,frequencies,noise_ratio=0.01)
    circuitfunc = tree_to_function(get_circuit_tree(circuit))
    impedances = [circuitfunc(circuit.parameters,fr) for fr in frequencies] 
    reals = real(impedances)
    imags = imag(impedances)
    return [R+rand(Normal(0,abs(R*noise_ratio)))+I*im + rand(Normal(0,abs(-I*noise_ratio)))im for (R,I) in zip(reals,imags)]
end

# function nyquist(measurements)
#     scatter(real(measurements),-imag(measurements))
# end