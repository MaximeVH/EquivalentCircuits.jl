"""
CircuitEvolution(measurements,frequencies,generations=1,population_size=100)

For a given set of impedance measurements and their corresponding input frequencies, a genetic programming
algorithm with a population size of population_size is run for generations iterations. 50% truncation selection
    is applied in each generation so that the initial population size is preserved. The resulting population
is returned.
"""
# function CircuitEvolution(measurements,frequencies,generations=1,population_size=30,terminals = "LRCP")
#     population = Generate_population(population_size,8,3,terminals)
#     population_fitness = evaluate_fitness(population,measurements,frequencies)
#     offspring = Array{LRC_Circuit}(undef, population_size)
#     offspring_fitness = Array{Float32}(undef, population_size)
#     for i in 1:generations
#         offspring = generate_offspring(population)
#         offspring_fitness = evaluate_fitness(offspring,measurements,frequencies)
#         population, population_fitness = TruncationSelection(population,offspring,population_fitness,offspring_fitness)
#     end
#     return population, population_fitness
# end


function CircuitEvolution2(measurements,frequencies,generations=1,population_size=30,n_terminals=4)
    population = Generate_population(population_size,6,3,n_terminals)
    population_fitness = evaluate_fitness(population,measurements,frequencies)
    offspring = Array{LRC_Circuit}(undef, population_size)
    offspring_fitness = Array{Float32}(undef, population_size)
    for i in 1:generations
        offspring = generate_offspring(population)
        offspring_fitness = evaluate_fitness(offspring,measurements,frequencies)
        population, population_fitness = TruncationSelection(population,offspring,population_fitness,offspring_fitness)
    end
    return population, population_fitness
end

function CircuitEvolution3(measurements,frequencies,generations=1,population_size=30,n_terminals=4)
    population = Generate_population(population_size,6,3,n_terminals)
    population_fitness = evaluate_fitness2(population,measurements,frequencies)
    offspring = Array{LRC_Circuit}(undef, population_size)
    offspring_fitness = Array{Float32}(undef, population_size)
    for i in 1:generations
        offspring = generate_offspring(population,n_terminals)
        offspring_fitness = evaluate_fitness2(offspring,measurements,frequencies)
        population, population_fitness = TruncationSelection(population,offspring,population_fitness,offspring_fitness)
    end
    return population, population_fitness
end




generations=1;population_size=30;n_terminals=4
population = Generate_population(population_size,6,3,n_terminals)
population_fitness = evaluate_fitness2(population,measurements,frequencies)
offspring = Array{LRC_Circuit}(undef, population_size)
offspring_fitness = Array{Float32}(undef, population_size)

# for i in 1:generations
offspring = generate_offspring(population,4)
offspring_fitness = evaluate_fitness2(offspring,measurements,frequencies)
# end



fitnesses = Array{Float32}(undef, length(population))
for (e,individual) in enumerate(offspring)
    OptParams, objective = CircuitFitness2(individual,measurements,frequencies)
    individual.parameters[individual.parameter_indices] = OptParams
    fitnesses[e] = objective
end

i += 1
OptParams, objective = CircuitFitness2(offspring[8],measurements,frequencies)

OptParams, objective = CircuitFitness(offspring[26],measurements,frequencies)


offspring[14]

for i in 1:30
    println(offspring[i])
end