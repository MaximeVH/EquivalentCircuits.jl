function get_removal_parameters(tree) #in order not to miss hierarchically higher series operations: first do an upward search of the tree going up parents with a while loop as long as only serial operators (so no parallel operator) are encountered.
    components_to_remove_params = [] #consider an additional sorting by parent (necesssary?).
    parallel_detection = false
    components = listcomponents(tree)
    if length(tree) == 3
        return components_to_remove_params
    end
    for component in components
        component_type = component.Type
        Parent = tree[component.ParentIndex]
        parentoperation = Parent.Type
        if !(component.Parameter in components_to_remove_params)
            sibling = Parent.Child1Index == component.Index ? tree[Parent.Child2Index] : tree[Parent.Child1Index]
            sibling_decendents = descendents(tree,sibling)
            descendant_counter = 1
            while !parallel_detection && descendant_counter <= length(sibling_decendents)
                if sibling_decendents[descendant_counter].Type == '-' 
                    parallel_detection = true
                elseif sibling_decendents[descendant_counter].Type == component_type
                    push!(components_to_remove_params,sibling_decendents[descendant_counter].Parameter)
                end
                descendant_counter += 1
            end
        end
        parallel_detection = false
    end
    if length(components_to_remove_params) + 1 == length(components)
        return []
    else
        return components_to_remove_params
    end
end

function simplifycircuit(tree)
    if length(tree) == 3
        return tree
    end
    simplified = false
    while !simplified
        removal_parameters = get_removal_parameters(tree)
        if isempty(removal_parameters)
            simplified = true
        else
            for P in removal_parameters
                if length(tree) > 3
                component = tree[findfirst(x->x.Parameter == P,tree)]
                tree = removecomponent(tree,component)
                end
            end
        end
    end
    return tree
end

function simplifycircuit(circuit::Circuit,terminals="RCPL")
    if count(isoperation,circuit.karva[1:3]) == 1
        return circuit
    end
    circuit_coding_length = length(circuit.karva)
    initialtree = karva_to_tree(circuit.karva,circuit.parameters)
    tree = simplifycircuit(initialtree)
    if tree == initialtree
        return circuit
    end
    rest_karva = circuit.karva[length(initialtree):end]
    rest_parameters = circuit.parameters[length(initialtree):end]
    extra_karva = join(rand(terminals,circuit_coding_length-length(tree)-length(rest_karva)))
    extra_parameters = karva_parameters(extra_karva) 
    return Circuit(join([node.Type for node in tree])*rest_karva*extra_karva,vcat(get_tree_parameters(tree),rest_parameters,extra_parameters),nothing)
end

function simplifycircuit!(circuit::Circuit,terminals="RCPL") 
    if count(isoperation,circuit.karva[1:3]) != 1
        circuit_coding_length = length(circuit.karva)
        initialtree = karva_to_tree(circuit.karva,circuit.parameters)
        tree = simplifycircuit(initialtree)
        if tree != initialtree
            rest_karva = circuit.karva[length(initialtree):end]
            rest_parameters = circuit.parameters[length(initialtree):end]
            extra_karva = join(rand(terminals,circuit_coding_length-length(tree)-length(rest_karva)))
            extra_parameters = karva_parameters(extra_karva) 
            circuit.karva, circuit.parameters =  join([node.Type for node in tree])*rest_karva*extra_karva, vcat(get_tree_parameters(tree),rest_parameters,extra_parameters)
        end
    end
end

function replace_redundant_cpes!(circuit::Circuit)
    cpe_positions = findall(x->x=='P',circuit.karva)
    mutable_karva = split(circuit.karva,"")
    for position in cpe_positions
        if circuit.parameters[position][2] > 0.99 #threshold to be adjusted. Capacitor replacement.
            mutable_karva[position] = "C"
            circuit.parameters[position] = 1/(circuit.parameters[position][1])
        elseif circuit.parameters[position][2] < 0.01 #threshold to be adjusted. Resistor replacement.
            mutable_karva[position] = "R"
            circuit.parameters[position] = circuit.parameters[position][1]
        end
    end
    circuit.karva = join(mutable_karva);
end

#Extra simplification to be included: if the second CPE parameter is sufficiently close to 0.5: it can be replaced by a Warburg element.

function replace_redundant_cpes!(population::Array{Circuit,1})
    for circuit in population
        replace_redundant_cpes!(circuit)
    end
end