function CircuitEvolution(Population,measurements,frequencies,fitnesses=Nothing,mutation_rate=0.2)
    Progenitors = [Population[rand(1:N,2)] for i in 1:N]
    Parents = Array{RC_Circuit}(undef, 2)
    Offspring = Array{RC_Circuit}(undef, N)

    for i in 1:N
               Parents = Progenitors[i]
               Child = CircuitCrossover(Parents[1],Parents[2])
               Mutant = CircuitMutation(Child)
               Offspring[i] = Mutant
    end

    if fitnesses == Nothing
        Population_fitnesses = [CircuitFitness(i,measurements,frequencies) for i in Population]
    else
        Population_fitnesses = fitnesses
    end

    Offspring_fitness = [CircuitFitness(i,measurements,frequencies) for i in Offspring]

    Total_population = vcat(Population,Offspring)
    Total_fitness = vcat(Population_fitnesses,Offspring_fitness)

    Next_generation = Total_population[sortperm(Total_fitness)[1:100]]
    Next_generation_fitness = Total_fitness[sortperm(Total_fitness)[1:100]]

    return
end
