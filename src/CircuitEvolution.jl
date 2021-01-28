function initializecircuit(head=8,terminals="RCLP")
    karva = generatekarva(head,terminals) 
    parameters = karva_parameters(karva)
    return Circuit(karva,parameters,nothing)
end

function initializepopulation(size=20,head=8,terminals="RCLP")
    return [initializecircuit(head,terminals) for i in 1:size]
end

function initializevariedpopulation(size=30,head=8)
    subpopulation_size = floor(Int(size/3))
    RCs = initializepopulation(subpopulation_size,head,"RC")
    RCLs = initializepopulation(subpopulation_size,head,"RCL")
    RCLPs = initializepopulation(subpopulation_size,head,"RCLP")
    return vcat(RCs,RCLs,RCLPs)
end

function circuit_offspring(circuit_1,circuit_2,terminals = "RCLP")
    offspring = ""
    rand_ = rand()
    if rand_ > 0.5
        offspring = one_point_crossover(circuit_1,circuit_2)
    elseif rand_ > 0.1
        offspring = two_point_crossover(circuit_1,circuit_2)
    else
        offspring = rand([circuit_1,circuit_2])
    end
    return mutate(offspring,terminals)
end

function circuitfitness(circuit,measurements,frequencies)
    tree = get_circuit_tree(circuit)
    circfunc,params,upper,param_inds = func_and_params_for_optim(tree)
    objective = objectivefunction(circfunc,measurements,frequencies)
    optparams,fitness = optimizeparameters(objective,params,upper) 
    return deflatten_parameters(optparams,tree,param_inds), fitness
end

function generate_offspring(population,fitnesses) 
    population_size = length(population)
    elite_size = Int(ceil(0.1*population_size))
    mating_pool_size = population_size-elite_size
    progenitors = tournamentselection(population,fitnesses,elite_size)
    elites = population[1:elite_size] 
    offspring = Array{Circuit}(undef, mating_pool_size)
    for e in 1:mating_pool_size
        offspring[e] = circuit_offspring(population[progenitors[e]],population[progenitors[mating_pool_size-e+1]])
    end
    vcat(elites,offspring)
end

function evaluate_fitness(population,measurements,frequencies)
    fitnesses = Array{Float64}(undef, length(population)) 
    for (e,individual) in enumerate(population)
        optParams, objective = circuitfitness(individual,measurements,frequencies)
        fitnesses[e] = objective
    end
    return fitnesses
end

function circuitevolution(measurements,frequencies,generations=1,population_size=30,terminals = "RCLP",head=8) 
    population = initializepopulation(population_size,head,terminals)
    fitnesses = evaluate_fitness(population,measurements,frequencies)
    sortingpermutation = sortperm(fitnesses)
    population = population[sortingpermutation]
    fitnesses = fitnesses[sortingpermutation]
    for i in 1:generations
        population = generate_offspring(population,fitnesses)
        fitnesses = evaluate_fitness(population,measurements,frequencies)
        sortingpermutation = sortperm(fitnesses)
        population = population[sortingpermutation]
        fitnesses = fitnesses[sortingpermutation]
    end
    return population, fitnesses
end