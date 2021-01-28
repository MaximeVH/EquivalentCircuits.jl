function simulateimpedance_noiseless(circuitfunc,parameters,frequencies) 
    return [circuitfunc(parameters,fr) for fr in frequencies] 
end

function simulateimpedance(circuitfunc,parameters,frequencies,noise_ratio=0.01)
    impedances = [circuitfunc(parameters,fr) for fr in frequencies] 
    reals = real(impedances)
    imags = imag(impedances)
    return [R+rand(Normal(0,abs(R*noise_ratio)))+I*im + rand(Normal(0,abs(-I*noise_ratio)))im for (R,I) in zip(reals,imags)]
end

function nyquist(measurements)
    scatter(real(measurements),imag(measurements))
end