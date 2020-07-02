function CircuitObjective(Circuit,measurements,frequencies)
    CircFunc = CircuitFunction(Circuit)
    ObjectFunc = ObjectiveFunction(CircFunc,measurements,frequencies)
    return ObjectFunc
end
