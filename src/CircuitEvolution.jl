"""
CircuitEvolution(measurements,frequencies,generations=1,population_size=100)

For a given set of impedance measurements and their corresponding input frequencies, a genetic programming
algorithm with a population size of population_size is run for generations iterations. 50% truncation selection
    is applied in each generation so that the initial population size is preserved. The resulting population
is returned.
"""
function CircuitEvolution(measurements,frequencies,generations=1,population_size=30,terminals = "LRC")
    population = Generate_population(population_size,8,3,terminals)
    population_fitness = evaluate_fitness(population,measurements,frequencies)
    offspring = Array{LRC_Circuit}(undef, population_size)
    offspring_fitness = Array{Float64}(undef, population_size)
    for i in 1:generations
        offspring = generate_offspring(population)
        offspring_fitness = evaluate_fitness(offspring,measurements,frequencies)
        population, population_fitness = TruncationSelection(population,offspring,population_fitness,offspring_fitness)
    end
    return population, population_fitness
end
