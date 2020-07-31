"""
CircuitEvolution(measurements,frequencies,generations=1,population_size=100)

For a given set of impedance measurements and their corresponding input frequencies, a genetic programming
algorithm with a population size of population_size is run for generations iterations. 50% truncation selection
    is applied in each generation so that the initial population size is preserved. The resulting population
is returned.
"""
function CircuitEvolution(measurements,frequencies,generations=1,population_size=100)
    Population = InitializePopulation(population_size)
    Offspring = Array{RC_Circuit}(undef, population_size)
    Progenitors = [Population[rand(1:population_size,2)] for i in 1:population_size]
    #Generate Offspring.
    for i in 1:population_size
        println(i)
        a = true
        Parents = Progenitors[i]
        while a == true
               Child = CircuitCrossover(Parents[1],Parents[2])
               Mutant = CircuitMutation(Child)
               if length(Mutant.parameters)>1 #We don't want univariate offspring.
                   Offspring[i] = Mutant
                   a = false
               end
    end
    end
    Population_fitnesses = [CircuitFitness(i,measurements,frequencies) for i in Population]
    Offspring_fitness = [CircuitFitness(i,measurements,frequencies) for i in Offspring]
    Total_population = vcat(Population,Offspring)
    Total_fitness = vcat(Population_fitnesses,Offspring_fitness)
    Selected = sortperm(Total_fitness)[1:population_size]
    Next_generation = Total_population[Selected]
    return Next_generation
end
