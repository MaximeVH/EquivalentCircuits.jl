function ObjectiveFunction(CircFunc,measurements,frequencies)
    function Objective(x)
        objective = 0
        model_output = [Base.invokelatest(CircFunc,x,fr) for fr in frequencies]
        rme,ime,rmo,imo = real(measurements),imag(measurements),real(model_output),imag(model_output)
        for i in 1:length(measurements)
            objective =+ ((rme[i] - rmo[i])^2 + (ime[i] - imo[i])^2)
        end
        return objective
    end
    return Objective
end
