function get_parameter_upper_bound(tree)
    ranges = Dict('R'=>1.0e9,'C'=>0.01,'L'=>5,'P'=>[1.0e9,1],'+'=>0,'-'=>0) #Dict('R'=>5000,'C'=>0.001,'L'=>1,'P'=>[100,1],'+'=>0,'-'=>0) 
    return [ranges[node.Type] for node in tree]
end

function func_and_params_for_optim(tree) 
    circuit,circuit_parameters,param_inds = tree_to_circuit_with_inds(tree)
    circuitfunc = circuitfunction(circuit)
    upperbounds = get_parameter_upper_bound(tree)[param_inds]
    return circuitfunc, flatten(circuit_parameters) , flatten(upperbounds) , param_inds
end

function optimizeparameters(objective,initial_parameters,upper)
    lower = zeros(length(initial_parameters))
    inner_optimizer = NelderMead() 
    results = optimize(objective, lower, upper, initial_parameters, Fminbox(inner_optimizer), Optim.Options(time_limit = 20.0))
    return results.minimizer,results.minimum 
end

function deflatten_parameters(parameters,tree,param_inds)
    parametets = 
    correct_value_lengths = length.(get_tree_parameters(tree)[param_inds])
    correct_length = length(correct_value_lengths)
    deflattened_parameters = Array{Any}(undef,length(correct_value_lengths))
    flat_index_counter = 1
    for (e,v) in enumerate(correct_value_lengths)
        if v == 1
            deflattened_parameters[e] = parameters[flat_index_counter]
            flat_index_counter += 1
        else
            deflattened_parameters[e] = [parameters[flat_index_counter],parameters[flat_index_counter+1]]
            flat_index_counter += 2
        end
    end
    return deflattened_parameters
end