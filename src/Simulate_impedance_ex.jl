function Simulate_impedance_ex(circuit_function,parameters,frequencies,noise_ratio)
    # Get precise impedance values
    # impedances = [Solve_Circuit(Circuit,parameters,fr) for fr in frequencies]
    impedances = [circuit_function(parameters,fr) for fr in frequencies]
    reals = real(impedances)
    imags = imag(impedances)
    imepdance_measurements = [] # Add a suitable type preallocation here to improve the performance.
    #Add gaussian (mean = zero) measurement noise
    for i in 1:length(impedances)
        push!(imepdance_measurements,reals[i]+rand(Normal(0,reals[i]*noise_ratio))+imags[i]im + rand(Normal(0,-imags[i]*noise_ratio))im)
    end
    return imepdance_measurements
end
