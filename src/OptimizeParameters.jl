"""
    OptimizeParameters(Objective,Initial_parameters)
Returns the optimal parameters for a given circuit, provided its Objective function (obtained by the ObjectiveFunction function)
and a set of Initial_parameters (as obtained by the InitializeParameters function).
"""


function OptimizeParameters(Objective,Initial_parameters)
    if length(Initial_parameters) == 1
        return [1]
    else
    lower = zeros(length(Initial_parameters))
    upper = [10e3 for i in 1:length(Initial_parameters)]
    inner_optimizer = NelderMead()
    results = optimize(Objective, lower, upper, Initial_parameters, Fminbox(inner_optimizer),Optim.Options(time_limit = 20.0))
    return results.minimizer
end
end

function OptimizeParameters2(Objective,Initial_parameters,upper)
    lower = zeros(length(Initial_parameters))
    inner_optimizer = NelderMead()
    results = optimize(Objective, lower, upper, Initial_parameters, Fminbox(inner_optimizer), Optim.Options(time_limit = 20.0))
    return results.minimizer
end


function get_parameter_upper_bound(circuit,flat_params)
    param_length = length(flat_params)
    ranges = Dict('a'=>1000,'b'=>0.001,'c'=>1,'d'=>[100,1]) #to be expanded.
    elements = get_elements(circuit.circuit)
    lower = zeros(param_length)
    upper = Array{Float64}(undef, param_length)
    index_counter = 1
    for e in elements
        to_add = ranges[e]
        if e == 'd'
            upper[index_counter] = to_add[1]
            upper[index_counter+1] = to_add[2]
            index_counter +=2
        else
            upper[index_counter] = to_add
            index_counter +=1
        end
    end
    return upper
end


get_parameter_upper_bound(offspring[8],flatten_parameters)

get_parameters(offspring[8])


offspring[8].parameters[offspring[8].parameter_indices]

get_readable_circuit(offspring[8].circuit)