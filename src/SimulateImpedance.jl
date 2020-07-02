function SimulateImpedance(CircuitFunc,parameters,frequencies,noise_ratio=0.01)
    # Get precise impedance values
    impedances = [CircuitFunc(parameters,fr) for fr in frequencies]
    reals = real(impedances)
    imags = imag(impedances)
    imepdance_measurements = Array{Any}(undef, length(frequencies))# Add a suitable type preallocation here to improve the performance.
    #Add gaussian (mean = zero) measurement noise
    for i in 1:length(impedances)
        imepdance_measurements[i] = reals[i]+rand(Normal(0,reals[i]*noise_ratio))+imags[i]im + rand(Normal(0,-imags[i]*noise_ratio))im
    end
    return imepdance_measurements
end
