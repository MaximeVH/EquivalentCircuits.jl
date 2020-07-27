"""
    ObjectiveFunction(CircFunc,measurements,frequencies)
Returns the Least squares objective function between the measurements and a circuit's simulated
impedance for a given set of circuit parameters at a specified range of experimental frequencies.
This function is used to optimize the circuit parameters.

"""
function ObjectiveFunction(CircFunc,measurements,frequencies)
    function Objective(x)
        N = length(measurements)
        model_output = [Base.invokelatest(CircFunc,x,fr) for fr in frequencies]
        rme,ime,rmo,imo = real(measurements),imag(measurements),real(model_output),imag(model_output)
        R = [(rme[d]-rmo[d])^2 for d in 1:N]
        I = [(ime[d]-imo[d])^2 for d in 1:N]
        objective = sum(R) + sum(I)
        return objective
    end
    return Objective
end
