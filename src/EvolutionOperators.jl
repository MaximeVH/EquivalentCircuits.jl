function one_point_crossover(circuit1,circuit2)
    karva1 = circuit1.karva
    karva2 = circuit2.karva
    crossover_point = rand(1:(length(karva1)-1)) 
    new_karva = karva1[1:crossover_point]*karva2[crossover_point+1:end]
    new_parameters =  vcat(circuit1.parameters[1:crossover_point],circuit2.parameters[crossover_point+1:end])
    return Circuit(new_karva,new_parameters,nothing)
end

function two_point_crossover(circuit1,circuit2)
    karva1 = circuit1.karva
    karva2 = circuit2.karva
    crossover_points = sort(rand(1:length(karva1)-1,2)) #if the two crossover points are the same, no crossing over occurs and the first circuit is returned.
    new_karva = karva1[1:crossover_points[1]]*karva2[crossover_points[1]+1:crossover_points[2]]*karva1[crossover_points[2]+1:end] 
    new_parameters = vcat(circuit1.parameters[1:crossover_points[1]],circuit2.parameters[crossover_points[1]+1:crossover_points[2]],circuit1.parameters[crossover_points[2]+1:end])
    return Circuit(new_karva,new_parameters,nothing)
end

function mutate(circuit,terminals = "RCLP")
    karva = circuit.karva ; parameters = copy(circuit.parameters)
    total_length = length(karva)
    head_length = (total_length-1)/2 ; tail_length = head_length+1
    operations = "+-"; new_element = ""
    position = rand(1:total_length)
    mutable_karva = collect(karva)
    if position == 1
        new_element = rand(operations)
        mutable_karva[position] = new_element
    elseif position <= head_length
        new_element = rand(operations*terminals)
        mutable_karva[position] = new_element
        parameters[position] = replacementparameter(new_element)
    elseif position > head_length
        new_element = rand(terminals)
        mutable_karva[position] = new_element
        parameters[position] = replacementparameter(new_element)
    end
    return Circuit(String(mutable_karva),parameters,nothing)
end

function replacementparameter(element)
    ranges= Dict('R'=>1000,'C'=>0.001,'L'=>1,'P'=>[100,1],'+'=>0,'-'=>0)
    return rand().*ranges[element]
end

function multipoint_mutation(circuit,N,terminals="RCLP")
    for i in 1:N
        circuit = mutate(circuit,terminals)
    end
    return circuit
end

function tournamentselection(population,fitnesses,elite_size)
    selected = []
    population_size = length(population)
    mating_pool_size = population_size - elite_size
    tournament_size = elite_size
    for i in 1:mating_pool_size
        tournament = rand(1:population_size,tournament_size)
        max_,ind = findmin(fitnesses[tournament])
        push!(selected,tournament[ind])
    end
    return selected
end