function new_objective(measurements,circuit_function,parameters,frequencies)
    # model_output = [circuit_function(parameters,fr) for fr in frequencies]
    model_output = [eval_circuit_function(circuit_function,parameters,fr) for fr in frequencies]
    objective = 0
    rme = real(measurements)
    ime = imag(measurements)
    rmo = real(model_output)
    imo = imag(model_output)
    for i in 1:length(measurements)
        objective =+ ((rme[i] - rmo[i])^2 + (ime[i] - imo[i])^2)
    end
    return objective
end

function eval_circuit_function(circuit_function,parameters,frequency)
    return Base.invokelatest(circuit_function,parameters,frequency)
    end
