"""
    CircuitObjective(Circuit,measurements,frequencies)

Creates the Objective function of a circuit (see ObjectiveFunction) for a given set
of measurements and measurement frequencies, starting from the circuit string instead of the Circuit function.
"""
function CircuitObjective(Circuit,measurements,frequencies)
    CircFunc = CircuitFunction(Circuit)
    ObjectFunc = ObjectiveFunction(CircFunc,measurements,frequencies)
    return ObjectFunc
end
