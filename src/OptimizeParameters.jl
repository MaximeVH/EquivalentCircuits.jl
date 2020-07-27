"""
    OptimizeParameters(Objective,Initial_parameters)
Returns the optimal parameters for a given circuit, provided its Objective function (obtained by the ObjectiveFunction function)
and a set of Initial_parameters (as obtained by the InitializeParameters function).
"""
function OptimizeParameters(Objective,Initial_parameters) #Avoid too much overhead for the user.
    lower = zeros(length(Initial_parameters))
    upper = [Inf for i in 1:length(Initial_parameters)]
    inner_optimizer = NelderMead()
    results = optimize(Objective, lower, upper, Initial_parameters, Fminbox(inner_optimizer))
    return results.minimizer
end
