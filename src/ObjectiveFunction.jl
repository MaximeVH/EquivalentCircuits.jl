"""
objectivefunction(circuitfunc,measurements,frequencies)
Returns the Least squares objective function between the measurements and a circuit's simulated
impedance for a given set of circuit parameters at a specified range of experimental frequencies.
This function is used to optimize the circuit parameters.

"""
function objectivefunction(circuitfunc,measurements,frequencies)
    function objective(x)
        N = length(measurements)
        model_output = [circuitfunc(x,fr) for fr in frequencies] 
        rme,ime,rmo,imo = real(measurements),imag(measurements),real(model_output),imag(model_output)
        return sum((rme.-rmo).^2 .+ (ime.-imo).^2) 
    end
    return objective
end
