function Optimal_circuit_parameters(circuit,init_params,frequencies,measurements)
    circuit_function = Create_Circuit_function(circuit)
    objective = Get_objective_function(impedance_measurements,circuit_function,frequencies)
    opt = optimize(objective,init_params,NelderMead())

    return opt.minimizer, opt.minimum
end

function Circuit_to_objective(Example_circuit,Measurements_2,Frequencies)
    circuit_function = Create_Circuit_function(Example_circuit)
    objective = Get_objective_function(Measurements_2,circuit_function,Frequencies)
    return objective
end
