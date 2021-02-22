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
    crossover_points = sort(rand(1:length(karva1)-1,2))
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
    ranges= Dict('R'=>1000,'C'=>0.001,'L'=>1,'P'=>[100,0.99999],'+'=>0,'-'=>0)
    return rand().*ranges[element]
end

function multipoint_mutation(circuit,N,terminals="RCLP")
    for i in 1:N
        circuit = mutate(circuit,terminals)
    end
    return circuit
end

function transposition(circuit)
    parameters = copy(circuit.parameters)
    coding_length = length(circuit.karva)
    startsite = rand(2:coding_length-1)
    max_transposonlength = coding_length-startsite
    transposonlength = rand(1:max_transposonlength)
    targetsite = rand(2:startsite)
    mutable_karva = collect(circuit.karva)
    mutable_karva[targetsite:targetsite+transposonlength] = mutable_karva[startsite:startsite+transposonlength]
    parameters[targetsite:targetsite+transposonlength] = parameters[startsite:startsite+transposonlength]
    return  Circuit(String(mutable_karva),parameters,nothing)
end

function tournamentselection(population,mating_pool_size,tournament_size)
    selected = []
    for i in 1:mating_pool_size
        tournament = rand(population,tournament_size)
        push!(selected,minimum(tournament))
    end
    return selected
end

function truncationselection(population,truncation_size)
    sort!(population)
    return population[1:truncation_size]
end