function subcircuits(circuit,terminals = "RCLP")
    subcircs = Circuit[]; karva = ""; parameters = []; 
    sublength = 0
    coding_length = length(circuit.karva)
    tree = karva_to_tree(circuit.karva,circuit.parameters)
    saplings = subtrees(tree)
    for subtree in saplings
        karva = tree_to_karva(subtree,div((coding_length-1),2),terminals)
        parameters = get_tree_parameters(subtree)
        sublength = length(parameters)
        push!(subcircs,Circuit(karva,vcat(parameters,karva_parameters(karva[sublength+1:end])),nothing))
    end
    return subcircs
end

function subtrees(tree) 
    trees = []
    components = listcomponents(tree)
    for component in components
        push!(trees,removecomponent(tree,component))
    end
    return trees
end

function redundacy_testing(circuit,measurements,frequencies,terminals = "RCLP",cutoffratio = 0.80)
    if count(isoperation,circuit.karva[1:3]) < 2
        return circuit
    end
    subcircs = subcircuits(circuit,terminals)
    evaluate_fitness!(subcircs,measurements,frequencies)
    candidate = minimum(subcircs)
    fitnessratio = circuit.fitness/candidate.fitness
    if fitnessratio >= cutoffratio
        return candidate
    else
        return circuit
    end
end

function get_subcircuits(circuit,measurements,frequencies,terminals = "RCLP")
    subcircs = subcircuits(circuit,terminals)
    evaluate_fitness!(subcircs,measurements,frequencies) 
    return subcircs
end

function removeredundancy(circuit,measurements,frequencies,terminals = "RCPL",cutoff = 0.80)
    if count(isoperation,circuit.karva[1:3]) == 1
        return circuit
    end
    subcircs = get_subcircuits(circuit,measurements,frequencies,terminals)
    candidate = minimum(subcircs)
    redundancy = true
    while redundancy && count(isoperation,candidate.karva[1:3])>1
        if circuit.fitness/candidate.fitness > cutoff
            circuit = candidate
        else 
            redundancy = false
        end
        subcircs = get_subcircuits(circuit,measurements,frequencies,terminals)
        candidate = minimum(subcircs)
    end
    return circuit
end