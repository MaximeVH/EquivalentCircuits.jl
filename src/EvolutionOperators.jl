"""
    one_point_crossover(circuit1,circuit2)

Recombine two parent circuits, p1 and p2, by replacing a randomly selected
component or subcircuit of p1 by a randomly selected subcircuit (that can also be nothing)
of p2. The output is a new LRC_Circuit object.
"""
function one_point_crossover(circuit1,circuit2)
    karva1 = circuit1.karva
    karva2 = circuit2.karva
    crossover_point = rand(1:(length(karva1)-1))
    new_karva = karva1[1:crossover_point]*karva2[crossover_point+1:end]
    new_parameters =  vcat(circuit1.parameters[1:crossover_point],circuit2.parameters[crossover_point+1:end])
    ET,inds = Karva_to_ET(new_karva)
    new_circuit = ET_to_expression(number_ET(ET))
    return LRC_Circuit(new_karva,new_circuit,new_parameters,inds)
end

"""
    two_point_crossover(circuit1,circuit2)

Recombine two parent circuits, p1 and p2, by replacing a randomly selected
component or subcircuit of p1 by a randomly selected subciiiiircuit (that can also be nothing)
of p2. The output is a new LRC_Circuit object.
"""
function two_point_crossover(circuit1,circuit2)
    karva1 = circuit1.karva
    karva2 = circuit2.karva
    crossover_points = sort(rand(1:length(karva1)-1,2))
    new_karva = karva1[1:crossover_points[1]]*karva2[crossover_points[1]+1:crossover_points[2]]*karva1[crossover_points[2]+1:end]
    new_parameters = vcat(circuit1.parameters[1:crossover_points[1]],circuit2.parameters[crossover_points[1]+1:crossover_points[2]],circuit1.parameters[crossover_points[2]+1:end])
    ET,inds = Karva_to_ET(new_karva)
    new_circuit = ET_to_expression(number_ET(ET))
    return LRC_Circuit(new_karva,new_circuit,new_parameters,inds)
end
"""
    mutate(circuit,terminals = "LRC")

A given circuit p1 is subjected to mutation with a probability of mutation_probability.
mutation entails the generation of a new component/subcircuit that is randomly inserted in the
parent circuit.
"""
function mutate(circuit,terminals = "LRCP")
    #You could adjust later so that only operations can be replaced in first position.
    karva = circuit.karva
    parameters = circuit.parameters
    # parameters = circuit.parameters
    total_length = length(karva)
    head_length = (total_length-1)/2
    tail_length = head_length+1
    operations = "+-"
    # terminals = "LRC"
    position = rand(1:total_length)
    new_terminal_element = "";new_head_element=""
    origional_element = karva[position]
    mutable_karva = collect(karva)
    if position > head_length
        new_terminal_element = rand(terminals)
        mutable_karva[position] = new_terminal_element
        if new_terminal_element != 'P'
            parameters[position] = rand()*Dict('R'=>1000,'C'=> 0.001,'L'=> 1)[new_terminal_element]
        else
            parameters[position] = [rand()*100,rand()*1]
        end
    elseif position <= head_length
        new_head_element = rand(operations*terminals)
        mutable_karva[position] = new_head_element
        if new_head_element != 'P'
            parameters[position] = rand()*Dict('R'=>1000,'C'=> 0.001,'L'=> 1,'+'=> 0,'-'=> 0)[new_head_element]
        else
            parameters[position] = [rand()*100,rand()*1]
        end
    end
    new_karva = String(mutable_karva)
    # new_circuit = karva_to_circuit(new_karva)
    circuit.karva = new_karva
    circuit.parameters = parameters
    ET,inds = Karva_to_ET(new_karva)
    circuit.circuit = ET_to_expression(number_ET(ET))
    circuit.parameter_indices = inds
    return circuit
end
"""
    CircuitFitness(circuit,measurements,frequencies)

Outputs the the residual least squares objective between the
impedance measurements and the simulated impedance values of a given circuit with the input frequencies of the measurements.
the parameters of the circuit are first optimised with respect to the measurement values.
The output fitness value is a measure of how well a given circuit configuration is capable
of fitting the experimental data.
"""
function circuit_offspring(circuit_1,circuit_2)
    offspring = ""
    rand_ = rand()
    if rand_ > 0.5
        offspring = one_point_crossover(circuit_1,circuit_2)
    elseif rand_ > 0.1
        offspring = two_point_crossover(circuit_1,circuit_2)
    else
        offspring = rand([circuit_1,circuit_2])
    end
    return mutate(offspring)
end

function CircuitFitness_inclusive(circuit,measurements,frequencies)
    CircFunc = CircuitFunction_inclusive_Final(circuit.circuit)
    ObjFunc = ObjectiveFunction(CircFunc,measurements,frequencies)
    OptParams = OptimizeParameters(ObjFunc,to_flat_array(circuit.parameters[circuit.parameter_indices]))
    return deflatten_parameters(circuit.circuit,OptParams), ObjFunc(OptParams)
end

function generate_offspring(population)
    population_size = length(population)
    progenitors = [population[rand(1:population_size,2)] for i in 1:population_size]
    offspring = Array{LRC_Circuit}(undef, population_size)
    for (e,parents) in enumerate(progenitors)
        offspring[e] = circuit_offspring(parents[1],parents[2])
    end
    return offspring
end

function evaluate_fitness(population,measurements,frequencies)
    fitnesses = Array{Float32}(undef, length(population))
    for (e,individual) in enumerate(population)
        OptParams, objective = CircuitFitness_inclusive(individual,measurements,frequencies)
        individual.parameters[individual.parameter_indices] = OptParams
        fitnesses[e] = objective
    end
    return fitnesses
end

function TruncationSelection(population,offspring,population_fitness,offspring_fitness)
    Total_population = vcat(population,offspring)
    Total_fitness = vcat(population_fitness,offspring_fitness)
    Selected = sortperm(Total_fitness)[1:length(population)]
    Next_generation = Total_population[Selected]
    Next_generation_fitness = Total_fitness[Selected]
    return Next_generation,Next_generation_fitness
end

function multi_point_mutation(circuit,N)
    circuit = LRC_Circuit(circuit.karva,circuit.circuit,circuit.parameters,circuit.parameter_indices)
    for i in 1:N
        circuit = mutate(circuit)
    end
    return circuit
end

function segment_insertion(circuit)
    circuit = LRC_Circuit(circuit.karva,circuit.circuit,circuit.parameters,circuit.parameter_indices)
    karva = circuit.karva
    parameters = circuit.parameters
    total_length = length(karva)
    head_length = (total_length-1)/2
    tail_length = head_length+1
    operations = "+-"
    possible_insertions = 2:head_length
    possible_starts = 1:total_length
    start = rand(possible_starts)
    instertion = rand(possible_insertions)
    length_of_instertion = rand(1:(head_length-instertion))
    insertion_karva = circuit.karva[start:start+length_of_instertion]
    instertion_parameters = circuit.parameters[start:start+length_of_instertion]
    circuit.karva[instertion:instertion+length_of_instertion] = insertion_karva
    circuit.parameters[instertion:instertion+length_of_instertion] = instertion_parameters
    return circuit
end
