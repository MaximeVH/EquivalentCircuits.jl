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
