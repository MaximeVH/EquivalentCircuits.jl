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

function simplifypopulation!(population,terminals="RCPL")
    for circuit in population
        simplifycircuit!(circuit,terminals)
    end
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
    optparams = max.(min.(optparams,upper),0) #workaround of optim.jl's bug
    return deflatten_parameters(optparams,tree,param_inds), fitness, param_inds 
end

function generate_offspring(population,terminals="RCLP")
    population_size = length(population)
    elite_size = Int(ceil(0.1*population_size))
    mating_pool_size = population_size-elite_size
    progenitors = tournamentselection(population,mating_pool_size,elite_size)
    elites = population[1:elite_size]
    offspring = Array{Circuit}(undef, mating_pool_size)
    for e in 1:mating_pool_size
        offspring[e] = circuit_offspring(progenitors[e],progenitors[mating_pool_size-e+1],terminals)
    end
    vcat(elites,offspring)
end

function evaluate_fitness!(population,measurements,frequencies)
    params = []
    param_inds = []
    for circuit in population 
        params,circuit.fitness,param_inds = circuitfitness(circuit,measurements,frequencies)
        circuit.parameters[param_inds] = params
    end 
end

function circuitevolution(filepath::String;generations::Real=10,population_size=30,terminals = "RCLP",head=8,initial_population = nothing)
    meansurement_file = readdlm(filepath,',')
    reals = meansurement_file[:,1]
    imags = meansurement_file[:,2]
    frequencies = meansurement_file[:,3]
    measurements = reals + imags*im
    # population = initializepopulation(population_size,head,terminals) #initializevariedpopulation(population_size,head)
    if isnothing(initial_population)    
        population = initializepopulation(population_size,head,terminals) #initializevariedpopulation(population_size,head)
    else
        population = initial_population
    end
    simplifypopulation!(population) 
    evaluate_fitness!(population,measurements,frequencies)
    sort!(population)
    generation = 0
    convergence_threshold = 5e-4 
    elite = minimum(population)
    min_fitness = elite.fitness
    while (min_fitness > convergence_threshold) && (generation<=generations)
        population = generate_offspring(population,terminals)
        simplifypopulation!(population)
        evaluate_fitness!(population,measurements,frequencies)
        sort!(population)
        elite = minimum(population).fitness < elite.fitness ? minimum(population) : elite
        min_fitness = elite.fitness
        generation += 1
    end
    replace_redundant_cpes!(population)
    population = removeduplicates(sort!(vcat(population,elite)))
    for i in 1:3
        population[i] = removeredundancy(population[i],measurements,frequencies)
    end
    return readablecircuit.(population[1:5])
end

function circuitevolution(measurements,frequencies;generations::Real=10,population_size=30,terminals = "RCLP",head=8,initial_population=nothing)
    if isnothing(initial_population)  
        population = initializepopulation(population_size,head,terminals) #initializevariedpopulation(population_size,head)
    else
        population = initial_population
    end
        simplifypopulation!(population) 
        evaluate_fitness!(population,measurements,frequencies)
        sort!(population)
        generation = 0
        convergence_threshold = 5e-4 
        elite = minimum(population)
        min_fitness = elite.fitness
        while (min_fitness > convergence_threshold) && (generation<=generations)
            population = generate_offspring(population,terminals)
            simplifypopulation!(population)
            evaluate_fitness!(population,measurements,frequencies)
            sort!(population)
            elite = minimum(population).fitness < elite.fitness ? minimum(population) : elite
            min_fitness = elite.fitness
            generation += 1
        end
        replace_redundant_cpes!(population)
        population = removeduplicates(sort!(vcat(population,elite)))
        for i in 1:3
            population[i] = removeredundancy(population[i],measurements,frequencies)
        end
        return readablecircuit.(population[1:5])
    end

function visualizesolutions(measurements,frequencies,population)
    fig = scatter(real(measurements),-imag(measurements), label = "impedance measurements",markershape = :diamond,markersize = 8, title = "Top evolved circuits",legend=:outertopright,size = (1000, 500))
    for n in 1:5
        a = simulateimpedance_noiseless(population[n],frequencies)
        scatter!(real(a),-imag(a),label = "$(n).  $(readablecircuit(population[n]))")
    end
    return fig
end

function visualizesolutions(circuit::Circuit,frequencies,population)
    measurements = simulateimpedance_noiseless(circuit,frequencies)
    fig = scatter(real(measurements),-imag(measurements), label = "ground truth: $(readablecircuit(circuit))",markershape = :diamond,markersize = 8, title = "Top evolved circuits",legend=:outertopright,size = (1000, 500))
    for n in 1:3
        a = simulateimpedance_noiseless(population[n],frequencies)
        scatter!(real(a),-imag(a),label = "$(n).  $(readablecircuit(population[n]))")
    end
    return fig
end

