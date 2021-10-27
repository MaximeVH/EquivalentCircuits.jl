mutable struct Circuit 
    karva::String
    parameters::Array{Any,1}
    fitness::Union{Nothing,AbstractFloat}
end

struct TreeNode{I<:Integer,P<:Real}
    Index::I 
    ParentIndex::Union{Nothing,I}
    Child1Index::Union{Nothing,I}
    Child2Index::Union{Nothing,I}
    Type::Char
    Parameter::Union{P,Array{P,1}}
end

isless(a::TreeNode, b::TreeNode) = isless(a.Index, b.Index)
isless(a::Circuit, b::Circuit) = isless(a.fitness, b.fitness) 

function number_circuit(circuit)
    l_positions = findall(isterminal ,circuit) 
    numbers = [string(i) for i in 1:length(l_positions)]
    offsets = cumsum([length(n) for n in numbers])
    l_positions[2:end] .+= offsets[1:end-1]
    for (p,n) in zip(l_positions,numbers)
        circuit = circuit[1:p]*n*circuit[p+1:end]
    end
    return circuit
end

isoperation(chr) = occursin(chr,"-+")
isterminal(chr) = occursin(chr,"RCLP")

function generatekarva(head,terminals="RCLP")
    operators =  "+-+-"
    all_elements = operators*terminals
    karva = rand(operators)*String(rand(all_elements,head-1))*String(rand(terminals,head+1))
    return karva
end

function karva_parameters(karva) 
    parameters = Array{Any}(undef, length(karva))
    ranges = Dict('R'=>4000,'C'=> 0.0001,'L'=> 1,'+'=> 0,'-'=> 0,'P'=>[4000,1], 'W'=>4000)
    for (e,i) in enumerate(karva)
            if e == 'P'
                Element_values[n] = [ranges[e][1]*rand(),ranges[e][2]*rand()]
            else
            parameters[e] = rand()*ranges[i]
        end
    end
    return parameters
end

function karva_to_tree(karva,parameters)
    tree = []
    Root_node = TreeNode(1,nothing,2,3,karva[1],parameters[1])
    push!(tree,Root_node)
    child1_karva_index = 2 ; child2_karva_index = 3
    ChildNode = 4
    free_spaces = 2
    for i in 1:length(karva)
        free_spaces -= 1
        K = karva[i]
        if isoperation(K)
            free_spaces += 2
            if isoperation(karva[child1_karva_index])
                push!(tree,TreeNode(child1_karva_index,i,ChildNode,ChildNode+1,karva[child1_karva_index],parameters[child1_karva_index]))
                child1_karva_index += 2
                ChildNode += 2
            else
                push!(tree,TreeNode(child1_karva_index,i,nothing,nothing,karva[child1_karva_index],parameters[child1_karva_index]))
                child1_karva_index += 2
            end

            if isoperation(karva[child2_karva_index])
                push!(tree,TreeNode(child2_karva_index,i,ChildNode,ChildNode+1,karva[child2_karva_index],parameters[child2_karva_index]))
                child2_karva_index += 2
                ChildNode += 2
            else
                push!(tree,TreeNode(child2_karva_index,i,nothing,nothing,karva[child2_karva_index],parameters[child2_karva_index]))
                child2_karva_index += 2
            end
        end
        if free_spaces == 1
            break 
        end
    end
    return tree
end

function karva_to_tree(karva)
    parameters = karva_parameters(karva)
    tree = []
    Root_node = TreeNode(1,nothing,2,3,karva[1],parameters[1])
    push!(tree,Root_node)
    child1_karva_index = 2 ; child2_karva_index = 3
    ChildNode = 4
    free_spaces = 2
    for i in 1:length(karva)
        free_spaces -= 1
        K = karva[i]
        if isoperation(K)
            free_spaces += 2
            if isoperation(karva[child1_karva_index])
                push!(tree,TreeNode(child1_karva_index,i,ChildNode,ChildNode+1,karva[child1_karva_index],parameters[child1_karva_index]))
                child1_karva_index += 2
                ChildNode += 2
            else
                push!(tree,TreeNode(child1_karva_index,i,nothing,nothing,karva[child1_karva_index],parameters[child1_karva_index]))
                child1_karva_index += 2
            end

            if isoperation(karva[child2_karva_index])
                push!(tree,TreeNode(child2_karva_index,i,ChildNode,ChildNode+1,karva[child2_karva_index],parameters[child2_karva_index]))
                child2_karva_index += 2
                ChildNode += 2
            else
                push!(tree,TreeNode(child2_karva_index,i,nothing,nothing,karva[child2_karva_index],parameters[child2_karva_index]))
                child2_karva_index += 2
            end
        end
        if free_spaces == 1
            break
        end
    end
    return tree
end

get_circuit_tree(circuit::Circuit) = karva_to_tree(circuit.karva,circuit.parameters)

listcomponents(tree) = tree[findall(x->isterminal(x.Type),tree)]
listoperations(tree) = tree[findall(x->isoperation(x.Type),tree)]

function removecomponent(tree,node_to_remove)
    Parent = tree[node_to_remove.ParentIndex]
    new_tree_legnth = (length(tree)-2)
    sibling = Parent.Child1Index == node_to_remove.Index ? tree[Parent.Child2Index] : tree[Parent.Child1Index]
    all_indices = collect(1:new_tree_legnth)
    old_to_new_location_dict = Dict(Parent.Index=>Parent.ParentIndex,sibling.Index => Parent.Index)
    children_dict = Dict()
    new_tree = Array{TreeNode}(undef, new_tree_legnth)
    new_tree[Parent.Index] = tree[sibling.Index]
    edited_old_tree = reverse(tree[findall(x-> !(x.Index in [Parent.Index,sibling.Index,node_to_remove.Index]),tree)])
    splice!(all_indices, Parent.Index) 
    rev_inds = reverse(all_indices)
    for (i,j) in zip(rev_inds,edited_old_tree)
        new_tree[i] = j
        old_to_new_location_dict[j.Index] = i
    end
    correct_final_tree = Array{TreeNode}(undef, new_tree_legnth)
    for (e,node) in enumerate(new_tree)
        if isterminal(node.Type)
            correct_final_tree[e] = TreeNode(e,old_to_new_location_dict[node.ParentIndex],nothing,nothing,node.Type,node.Parameter)
            haskey(children_dict,old_to_new_location_dict[new_tree[e].ParentIndex]) ? push!(children_dict[old_to_new_location_dict[node.ParentIndex]],e) : children_dict[old_to_new_location_dict[node.ParentIndex]] = [e]
        else
            if node.Index != 1
                haskey(children_dict,old_to_new_location_dict[node.ParentIndex]) ? push!(children_dict[old_to_new_location_dict[node.ParentIndex]],e) : children_dict[old_to_new_location_dict[node.ParentIndex]] = [e]
            end
        end
    end
    for e in 1:new_tree_legnth
        if !isterminal(new_tree[e].Type)
            if e == 1
                correct_final_tree[e] = TreeNode(e,nothing,2,3,new_tree[e].Type,new_tree[e].Parameter)
            else
                correct_final_tree[e] = TreeNode(e,old_to_new_location_dict[new_tree[e].ParentIndex],children_dict[e][1],children_dict[e][2],new_tree[e].Type,new_tree[e].Parameter)
            end
        end
    end
    return correct_final_tree
end

function descendents(Tree,Node) 
    if isterminal(Node.Type)
        return [Node]
    else
        descendent1 = Tree[Node.Child1Index]
        descendent2 = Tree[Node.Child2Index]
        if all(isterminal.([descendent1.Type,descendent2.Type]))
            return [Node,descendent1,descendent2]
        else 
            return sort(vcat(Node,descendents(Tree,descendent1),descendents(Tree,descendent2)))
        end
    end
end

function readablecircuit(circuit::Circuit)
    tree = karva_to_tree(circuit.karva,circuit.parameters)
    return tree_to_circuit(tree)[1]
end

function tree_to_karva(tree,head=8,terminals = "RCPL")
    coding_karva = join([node.Type for node in tree])
    remaining_length = 2*head+1 - length(coding_karva)
    karva = coding_karva*join(rand(terminals,remaining_length))
    return karva
end

function tree_to_circuit(tree) 
    nodecount = length(tree)
    essential_info = [[node.ParentIndex,node.Type,[node.Parameter]] for node in tree]
    for i in 1:((length(essential_info)-1)/2)
        parent1,type1,params1 = pop!(essential_info)
        parent2,type2,params2 = pop!(essential_info)
        essential_info[parent1][2] = essential_info[parent1][2] == '+' ? type1*'-'*type2 : '['*type1*','*type2*']'
        essential_info[parent1][3] = vcat(params1,params2)
    end
return  number_circuit(essential_info[1][2]),  essential_info[1][3] 
end

denumber_circuit(Circuit) = replace(Circuit, r"[0-9]" => "")

head_length(karva) = findlast(x->isoperation(x),karva)

function get_karva_elements_and_parameters(circuit,circuit_parameters)
    circuit_bare = foldl(replace,["["=>"","]"=>"","-"=>"+",","=>"-"],init = denumber_circuit(circuit))
    karva_params = Array{Any}(undef,length(circuit_bare))
    counter = 1
    for (e,i) in enumerate(circuit_bare)
        if isterminal(i)
            karva_params[e] = circuit_parameters[counter]
                counter +=1
        else
            karva_params[e] = 0
        end
    end
    operation_indexes = findall(karva_params .== 0)
    terminal_indexes = findall(karva_params .!= 0)
    return circuit_bare, karva_params, operation_indexes,terminal_indexes
end

function exhaustive_search(circuit,circuit_parameters)
    # Get simulated measurement from the circuit
    circuit_bare, karva_params, operation_indexes,terminal_indexes = get_karva_elements_and_parameters(circuit,circuit_parameters)
    simulated_measurement = get_target_impedance(circuit,circuit_parameters)
   
    # Get all possible array permutations (the next few parts can be swapped for a heuristic if there are too many permutations).
    index_permutations = collect(permutations(1:length(circuit_bare)))
    # Remove permutations that lead to synctactically incorrect circuits.
    permitted_index_permutations = index_permutations[findall(x-> (x[1] in operation_indexes)&&(x[end] in terminal_indexes)&&(x[end-1] in terminal_indexes),index_permutations)]
    # Generate the arrays of potential solution.
    potential_karvas = [circuit_bare[inds] for inds in permitted_index_permutations]
    potential_paramsets = [karva_params[inds] for inds in permitted_index_permutations] 
    potential_trees = [karva_to_tree(k,p) for (k,p) in zip(potential_karvas,potential_paramsets)]
    potential_circuits = [tree_to_circuit(t)[1] for t in potential_trees]

    # Limit the potential solutions to those with a correct circuit configuration.
    correct_circuits = findall(potential_circuits .== circuit)
    potential_karvas = potential_karvas[correct_circuits]
    potential_paramsets = potential_paramsets[correct_circuits]
    potential_trees = potential_trees[correct_circuits]
   
    # Limit the potential solutions to those with a correct function evaluation.
    potential_functions = [tree_to_function(T) for T in potential_trees];
    function_evaluations = [f(params,1000) for (f,params) in zip(potential_functions,potential_paramsets)] 
    correct_evaluations = findall(function_evaluations .== simulated_measurement)
    potential_karvas = potential_karvas[correct_evaluations]
    potential_paramsets = potential_paramsets[correct_evaluations]
    potential_trees = potential_trees[correct_evaluations]

    # Obtain solution(s) with the sortest head. Maybe randomly select one of them.
    potential_heads = head_length.(potential_karvas)
    potential_head_minimum = min(potential_heads...)
    minimal_heads = findall(potential_heads .== potential_head_minimum)
    karva_solutions = potential_karvas[minimal_heads]
    parameter_solutions = potential_paramsets[minimal_heads]

    return karva_solutions[1] , parameter_solutions[1] 
end

function genetic_algorithm(circuit,circuit_parameters,pop_size = 50) #Maybe try to reimplement this with the Evolutionary package.
    # Preliminaries
    circuit_bare, karva_params, operation_indexes,terminal_indexes = get_karva_elements_and_parameters(circuit,circuit_parameters)
    target_impedance = get_target_impedance(circuit,circuit_parameters)
    # Generate initial population.
    Candidate_index_sets = generate_candidate_index_sets(pop_size,length(karva_params),operation_indexes,terminal_indexes)
    # Evaluate initial population's fitness.
    fitnesses = [indexset_fitness(Candidate_index_sets[i],circuit_bare,karva_params,target_impedance) for i in 1:pop_size]
    # Get the sorted fitnesses and candidate indexsets sorted by fitness.
    order = sortperm(fitnesses)
    Candidate_index_sets = Candidate_index_sets[order]
    fitnesses = fitnesses[order]
    best_possible_fitness = count(z->z==0,karva_params)
    if fitnesses[1] == best_possible_fitness 
        return circuit_bare[Candidate_index_sets[1]], karva_params[Candidate_index_sets[1]]
    end
    max_fitness = fitnesses[1]
    unchanged_generation_counter = 0
    while fitnesses[1] != best_possible_fitness && !(unchanged_generation_counter >= 20 && max_fitness < 10) 
        Candidate_index_sets,fitnesses = next_generation(Candidate_index_sets,fitnesses,circuit_bare,karva_params,operation_indexes,terminal_indexes,target_impedance,0.05)
        println(max_fitness)
        if fitnesses[1] >= max_fitness
            unchanged_generation_counter += 1
        else
            max_fitness = fitnesses[1]
            unchanged_generation_counter = 0
        end
        if unchanged_generation_counter >= 30 #restart if you get stuck in a local minimum. as long as fitness>10, the algorithm doesnt really have any direction, maybe an additional distance measure should be considered.
            fitnesses = [indexset_fitness(Candidate_index_sets[i],circuit_bare,karva_params,target_impedance) for i in 1:pop_size]
            order = sortperm(fitnesses)
            Candidate_index_sets = Candidate_index_sets[order]
            fitnesses = fitnesses[order]
            best_possible_fitness = count(z->z==0,karva_params)
        end
    end
    return circuit_bare[Candidate_index_sets[1]], karva_params[Candidate_index_sets[1]]
end

function generate_indexset(size,operation_indexes,terminal_indexes)
    initial_operation_constraint = false
    terminal_tail_constraint = false
    indexset = randperm(size)
    while !initial_operation_constraint || !terminal_tail_constraint
        initial_operation_constraint = false
        terminal_tail_constraint = false
        if indexset[1] in operation_indexes
            initial_operation_constraint = true
            if indexset[end-1] in terminal_indexes && indexset[end] in terminal_indexes
                terminal_tail_constraint = true
            else indexset = randperm(size)
            end
        else indexset = randperm(size)
        end

    end
    return indexset
end

function generate_candidate_index_sets(population_size,size,operation_indexes,terminal_indexes) 
    return [generate_indexset(size,operation_indexes,terminal_indexes) for i in 1:population_size]
end

function indexset_fitness(indexset,circuit_bare,karva_params,target_impedance)
    candidate_solution = circuit_bare[indexset]
    candidate_parameters = karva_params[indexset]
    candidate_function = tree_to_function(karva_to_tree(candidate_solution,candidate_parameters))
    candidate_impedance = candidate_function(candidate_parameters,1000)

    fitness = candidate_impedance == target_impedance ? head_length(candidate_solution) : 10
 
    return fitness
end

function next_generation(Candidate_index_sets,fitnesses,circuit_bare,karva_params,operation_indexes,terminal_indexes,target_impedance,mutationRate=0.05)
    #SELECTION
    individual_size = length(karva_params)
    pop_size = length(fitnesses)
    elite_size = Int(ceil(pop_size*0.1))
    children = Candidate_index_sets[1:elite_size] #The elites are automatically accepted into the new generation (should they be exempt from mutation as well?)
    Mating_pool = get_parent_inds_from_tournamant(pop_size,fitnesses,elite_size)
    Mating_pool_size = length(Mating_pool)
    #BREEDING
    for i in 1:Mating_pool_size
        parent1,parent2 = Candidate_index_sets[[Mating_pool[i],Mating_pool[Mating_pool_size-i+1]]]
        push!(children,ordered_restricted_crossover(individual_size,parent1,parent2,operation_indexes,terminal_indexes))
    end
    #MUTATION
    children = [restricted_swap_mutation(child,mutationRate,operation_indexes,terminal_indexes) for child in children]
    #FITNESS EVAULATION
    child_fitnesses = [indexset_fitness(children[i],circuit_bare,karva_params,target_impedance) for i in 1:pop_size]
    order = sortperm(child_fitnesses)
    
    return children[order], child_fitnesses[order]
end

function get_parent_inds_from_tournamant(pop_size,fitnesses,elite_size)
    parents = []
    mating_pool_size = pop_size-elite_size
    tournament_size = Int(ceil(0.1*pop_size))
    for i in 1:mating_pool_size
    tournament = rand(1:pop_size,tournament_size)
    max_,ind = findmin(fitnesses[tournament])
    push!(parents,tournament[ind])
    end
    return parents
end

function restricted_swap_mutation(Candidate,mutationrate,operation_indicess,terminal_indicess)
    individual = copy(Candidate)
    individual_length = length(individual)
    swapWith = 0
    swapped = 0
    for swapped in 1:individual_length
        if rand() < mutationrate
            if swapped == 1
                swapWith = rand(findall(x->x in operation_indicess,individual))
            elseif swapped in [individual_length-1,individual_length]
                swapWith = rand(findall(x->x in terminal_indicess,individual))
            else
                if swapped in findall(x->x in operation_indicess,individual)
                    swapWith = rand(1:individual_length-2)
                else
                    swapWith = rand(2:individual_length)
                end
                
            end
            gene1 = individual[swapped]
            gene2 = individual[swapWith]
            individual[swapped] = gene2
            individual[swapWith] = gene1
        end
    end
    return individual
end

function ordered_crossover(size,parent1,parent2)
    geneA = Int(ceil(rand()*size))
    geneB = Int(ceil(rand()*size))
    StartGene = min(geneA,geneB)
    EndGene = max(geneA,geneB)
    child = Int[]
    parent1_range = StartGene:EndGene
    parent1 = parent1[StartGene:EndGene]
    parent2 = parent2[findall(x->!(x in parent1),parent2)]
    for i in 1:size
        if i in parent1_range
            push!(child,popfirst!(parent1))
        else
            push!(child,popfirst!(parent2))
        end
    end
    return child
end

function ordered_restricted_crossover(size,parent1,parent2,operation_indexes,terminal_indexes)
    valid_crossover = false
    child = ordered_crossover(size,parent1,parent2)
    while !valid_crossover
        if child[1] in operation_indexes && child[end-1] in terminal_indexes && child[end] in terminal_indexes
            valid_crossover = true
        else
            child = ordered_crossover(size,parent1,parent2)
        end
    end
    return child
end

function circuit_to_karva(circuit,circuit_parameters,head=8)
    karva = ""
    params = []
    if length(foldl(replace,["["=>"","]"=>"","-"=>"+",","=>"-"], init = denumber_circuit(circuit))) <= 9
        karva,params = exhaustive_search(circuit,circuit_parameters)
    else
        karva,params = genetic_algorithm(circuit,circuit_parameters)
    end
    amount_to_extend = (2*head+1) - length(karva)
    extension_karva = String(rand("RCLP",amount_to_extend))
    extension_params = karva_parameters(extension_karva)
    return karva*extension_karva , vcat(params,extension_params)
end

function get_tree_parameters(tree)
    return [node.Parameter for node in tree]
end

function tree_to_circuit_with_inds(tree) 
    nodecount = length(tree)
    essential_info = [[node.ParentIndex,node.Type,node.Index,[node.Parameter]] for node in tree]
    for i in 1:((length(essential_info)-1)/2)
        parent1,type1,index1,params1 = pop!(essential_info)
        parent2,type2,index2,params2 = pop!(essential_info)
        essential_info[parent1][2] = essential_info[parent1][2] == '+' ? type1*'-'*type2 : '['*type1*','*type2*']'
        essential_info[parent1][4] = vcat(params1,params2)
        essential_info[parent1][3] = vcat(index1,index2)
    end
return  number_circuit(essential_info[1][2]),  essential_info[1][4] , essential_info[1][3]
end