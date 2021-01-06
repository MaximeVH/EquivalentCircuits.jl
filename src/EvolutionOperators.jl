"""
    one_point_crossover(circuit1,circuit2)

Recombine two parent circuits, p1 and p2, by replacing a randomly selected
component or subcircuit of p1 by a randomly selected subcircuit (that can also be nothing)
of p2. The output is a new LRC_Circuit object.
"""
function one_point_crossover(circuit1,circuit2)
    karva1 = circuit1.karva
    karva2 = circuit2.karva
    crossover_point = rand(2:(length(karva1)-1)) #changed from 1
    new_karva = karva1[1:crossover_point]*karva2[crossover_point+1:end]
    new_parameters =  vcat(circuit1.parameters[1:crossover_point],circuit2.parameters[crossover_point+1:end])
    ET,inds = Karva_to_ET(new_karva)
    return LRC_Circuit(new_karva,ET,new_parameters,inds)
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
    crossover_points = sort(rand(2:length(karva1)-1,2))
    new_karva = karva1[1:crossover_points[1]]*karva2[crossover_points[1]+1:crossover_points[2]]*karva1[crossover_points[2]+1:end]
    new_parameters = vcat(circuit1.parameters[1:crossover_points[1]],circuit2.parameters[crossover_points[1]+1:crossover_points[2]],circuit1.parameters[crossover_points[2]+1:end])
    ET,inds = Karva_to_ET(new_karva)
    return LRC_Circuit(new_karva,ET,new_parameters,inds)
end
"""
    mutate(circuit,terminals = "LRC")

A given circuit p1 is subjected to mutation with a probability of mutation_probability.
mutation entails the generation of a new component/subcircuit that is randomly inserted in the
parent circuit.
"""
function mutate(circ,n_terms = 6) #adjust later to allow for a tailored selection of terminals.
    #You could adjust later so that only operations can be replaced in first position.
    @assert n_terms <= 26
    circuit = LRC_Circuit(circ.karva,circ.circuit,circ.parameters,circ.parameter_indices)
    karva = circuit.karva
    parameters = circuit.parameters
    total_length = length(karva)
    head_length = (total_length-1)/2
    tail_length = head_length+1
    operations = "+-"
    terminals = "abcdefghijklmnopqrstuvwxyz"[1:n_terms]
    position = rand(2:total_length)
    new_terminal_element = "";new_head_element=""
    origional_element = karva[position]
    mutable_karva = collect(karva)
    if position > head_length
        new_terminal_element = rand(terminals)
        mutable_karva[position] = new_terminal_element
        parameters[position] = get_replacement_parameters(new_terminal_element)
    elseif position <= head_length
        new_head_element = rand(operations*terminals)
        mutable_karva[position] = new_head_element
        if occursin(new_head_element, "+-")
            parameters[position] = 0
        else
            parameters[position] = get_replacement_parameters(new_head_element)
        end
    end
    new_karva = String(mutable_karva)
    circuit.karva = new_karva
    circuit.parameters = parameters
    ET,inds = Karva_to_ET(new_karva)
    circuit.circuit = ET
    circuit.parameter_indices = inds
    return circuit
end

function get_replacement_parameters(k)
    ranges=Dict('a'=>1000,'b'=>0.001,'c'=>1,'d'=>[100,1])
    params = []
    if k in ['a','b','c']
        push!(params,rand()*ranges[k])
    elseif k == 'd'
        push!(params,[rand()*ranges['d'][1],rand()*ranges['d'][2]])
    elseif k in ['e','h'] #RC
        push!(params,[rand()*ranges['a'],rand()*ranges['b']])
    elseif k in ['f','i'] #RL
        push!(params,[rand()*ranges['a'],rand()*ranges['c']])
    elseif k in ['k','m'] #CL
        push!(params,[rand()*ranges['b'],rand()*ranges['c']])
    elseif k in ['g','j'] #RP 
        push!(params,[rand()*ranges['a'],[rand()*ranges['d'][1],rand()*ranges['d'][2]]])
    elseif k in ['l','n'] #CP 
        push!(params,[rand()*ranges['b'],[rand()*ranges['d'][1],rand()*ranges['d'][2]]])
    elseif k == 'p' #LP 
        push!(params,[rand()*ranges['c'],[rand()*ranges['d'][1],rand()*ranges['d'][2]]])
    end
return params[1]
end
"""
    CircuitFitness(circuit,measurements,frequencies)

Outputs the the residual least squares objective between the
impedance measurements and the simulated impedance values of a given circuit with the input frequencies of the measurements.
the parameters of the circuit are first optimised with respect to the measurement values.
The output fitness value is a measure of how well a given circuit configuration is capable
of fitting the experimental data.
"""
# function circuit_offspring(circuit_1,circuit_2)
#     offspring = ""
#     rand_ = rand()
#     if rand_ > 0.5
#         offspring = one_point_crossover(circuit_1,circuit_2)
#     elseif rand_ > 0.1
#         offspring = two_point_crossover(circuit_1,circuit_2)
#     else
#         offspring = rand([circuit_1,circuit_2])
#     end
#     return mutate(offspring)
# end

function circuit_offspring(circuit_1,circuit_2,n_terms)
    offspring = ""
    rand_ = rand()
    if rand_ > 0.5
        offspring = one_point_crossover(circuit_1,circuit_2)
    elseif rand_ > 0.1
        offspring = two_point_crossover(circuit_1,circuit_2)
    else
        offspring = rand([circuit_1,circuit_2])
    end
    return mutate(offspring,n_terms)
end

function CircuitFitness(circuit,measurements,frequencies)
    CircFunc = Solve_tree(circuit.circuit)
    ObjFunc = ObjectiveFunction(CircFunc,measurements,frequencies)
    OptParams = OptimizeParameters(ObjFunc,flatten_parameters(circuit.parameters[circuit.parameter_indices]))
    return deflatten_parameters(circuit.circuit,OptParams), ObjFunc(OptParams)
end

function CircuitFitness2(circuit,measurements,frequencies)
    CircFunc = Solve_tree(circuit.circuit)
    ObjFunc = ObjectiveFunction(CircFunc,measurements,frequencies)
    init_parameters = get_parameters(circuit)
    upper_bound = get_parameter_upper_bound(circuit,init_parameters)
    OptParams = OptimizeParameters2(ObjFunc,init_parameters,upper_bound)
    return deflatten_parameters(circuit.circuit,OptParams), ObjFunc(OptParams)
end

function generate_offspring(population,n_terms)
    population_size = length(population)
    progenitors = [population[rand(1:population_size,2)] for i in 1:population_size]
    offspring = Array{LRC_Circuit}(undef, population_size)
    for (e,parents) in enumerate(progenitors)
        offspring[e] = circuit_offspring(parents[1],parents[2],n_terms)
    end
    return offspring
end

function evaluate_fitness(population,measurements,frequencies)
    fitnesses = Array{Float32}(undef, length(population))
    for (e,individual) in enumerate(population)
        OptParams, objective = CircuitFitness(individual,measurements,frequencies)
        individual.parameters[individual.parameter_indices] = OptParams
        fitnesses[e] = objective
    end
    return fitnesses
end

function evaluate_fitness2(population,measurements,frequencies)
    fitnesses = Array{Float32}(undef, length(population))
    for (e,individual) in enumerate(population)
        OptParams, objective = CircuitFitness2(individual,measurements,frequencies)
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
circuit = LRC_Circuit(circuit1.karva,circuit1.circuit,circuit1.parameters,circuit1.parameter_indices)
karva = collect(circuit.karva)
parameters = circuit.parameters
total_length = length(karva)
head_length = Int((total_length-1)/2)
tail_length = Int(head_length+1)
operations = "+-"
possible_insertions = 2:head_length
possible_starts = 1:total_length
start = rand(possible_starts)
insertion = rand(possible_insertions)#
length_of_insertion = rand(1:(head_length-insterion))
insertion_karva = karva[start:start+length_of_insertion]
insertion_parameters = circuit.parameters[start:start+length_of_insertion]
karva[insertion:insertion+length_of_insertion] = insertion_karva
circuit.karva=  String(karva)
circuit.parameters[insertion:insertion+length_of_insertion] = insertion_parameters
return circuit
end
