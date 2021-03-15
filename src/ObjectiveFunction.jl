"""
objectivefunction(circuitfunc,measurements,frequencies)
Returns the Least squares objective function between the measurements and a circuit's simulated
impedance for a given set of circuit parameters at a specified range of experimental frequencies.
This function is used to optimize the circuit parameters.

"""
function objectivefunction(circuitfunc,measurements,frequencies) 
    function objective(x)
        model_output = [circuitfunc(x,fr) for fr in frequencies] 
        return  mean((abs.(measurements - model_output).^2)./(abs.(measurements).^2 .+ abs.(model_output).^2))
    end
    return objective
end

function objectivefunction1(circuitfunc,measurements,frequencies) #fitness measure from Ramos et al.
    function objective(x)
        model_output = [circuitfunc(x,fr) for fr in frequencies] 
        return  mean((abs.(measurements - model_output).^2)./(abs.(model_output).^2))
    end
    return objective
end

#Some alternate objective functions.

function objectivefunction2(circuitfunc,measurements,frequencies)
    function objective(x)
        model_output = [circuitfunc(x,fr) for fr in frequencies] 
        rme,ime,rmo,imo = real(measurements),imag(measurements),real(model_output),imag(model_output)
        return sum((rme.-rmo).^2 .+ (ime.-imo).^2) 
    end
    return objective
end

function objectivefunction3(circuitfunc,measurements,frequencies) #fitness measure from Ramos et al.
    function objective(x)
        model_output = [circuitfunc(x,fr) for fr in frequencies] 
        return  mean((abs.(measurements - model_output).^2)./(abs.(measurements).^2))
    end
    return objective
end

function objectivefunction4(circuitfunc,measurements,frequencies) #fitness measure from Ramos et al.
    function objective(x)
        model_output = [circuitfunc(x,fr) for fr in frequencies] 
        return  mean((abs.(measurements - model_output).^2))
    end
    return objective
end

function objectivefunction5(circuitfunc,measurements,frequencies) #fitness measure from Lempka et al.
    function objective(x)
        model_output = [circuitfunc(x,fr) for fr in frequencies] 
        rme,ime,rmo,imo = real(measurements),imag(measurements),real(model_output),imag(model_output)
        return sum(((rme.-rmo).^2)./(rmo.^2 .+rme.^2) .+ ((ime.-imo).^2)./(imo.^2 .+rme.^2))
    end
    return objective
end

function objectivefunction6(circuitfunc,measurements,frequencies) #fitness measure from Lempka et al.
    function objective(x)
        model_output = [circuitfunc(x,fr) for fr in frequencies] 
        rme,ime,rmo,imo = real(measurements),imag(measurements),real(model_output),imag(model_output)
        return sum(((rme.-rmo).^2)./(rmo.^2) .+ ((ime.-imo).^2)./(imo.^2)) 
    end
    return objective
end 

function objectivefunction7(circuitfunc,measurements,frequencies)
    function objective(x)
        model_output = [circuitfunc(x,fr) for fr in frequencies] 
        return  mean((abs.(measurements - model_output).^2)./(abs.(measurements).^2 .+ abs.(model_output).^2))
    end
    return objective
end 