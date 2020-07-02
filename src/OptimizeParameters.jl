function OptimizeParameters(Objective,Initial_parameters) #Avoid too much overhead for the user.
    lower = zeros(length(Initial_parameters))
    upper = [Inf for i in 1:length(Initial_parameters)]
    inner_optimizer = NelderMead()
    results = optimize(Objective, lower, upper, Initial_parameters, Fminbox(inner_optimizer))
    return results.minimizer
end
