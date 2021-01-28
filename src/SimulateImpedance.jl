function simulateimpedance_noiseless(circuitfunc,parameters,frequencies) 
    return [circuitfunc(parameters,fr) for fr in frequencies] 
end

function simulateimpedance(circuitfunc,parameters,frequencies,noise_ratio=0.01)
    impedances = [circuitfunc(parameters,fr) for fr in frequencies] 
    reals = real(impedances)
    imags = imag(impedances)
    imepdance_measurements = Array{Any}(undef, length(frequencies))
    for i in 1:length(impedances)
        imepdance_measurements[i] = reals[i]+rand(Normal(0,abs(reals[i]*noise_ratio)))+imags[i]im + rand(Normal(0,abs(-imags[i]*noise_ratio)))im
    end
    return imepdance_measurements
end

function nyquist(measurements)
    scatter(real(measurements),imag(measurements))
end