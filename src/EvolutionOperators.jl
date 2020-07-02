function CircuitCrossover(p1,p2,max_depth=10)
    child = deepcopy(p1.expression)
    crosspoint = sample(p2.expression)
    typ = return_type(CircuitGrammar,crosspoint.ind)
    d_subtree = depth(crosspoint)
    # println(Quickstring(crosspoint))
    d_max = max_depth + 1 - d_subtree
    if d_max > 0 && contains_returntype(child,CircuitGrammar,typ,d_max)
        loc = sample(NodeLoc,child,typ,CircuitGrammar,d_max)
        insert!(child,loc,deepcopy(crosspoint))
    end
    InitializeCircuit(child)
end

function CircuitMutation(p1,mutation_probability=0.2)
    child = deepcopy(p1.expression)
    if rand() < mutation_probability
        loc = sample(NodeLoc,child)
        typ = return_type(CircuitGrammar,get(child,loc).ind)
        alt_subtree = CircuitExpression(typ)
        insert!(child,loc,alt_subtree)
    end
    return InitializeCircuit(child)
end

function CircuitFitness(circuit,measurements,frequencies)
    CircFunc = CircuitFunction(circuit.string)
    ObjFunc = ObjectiveFunction(CircFunc,measurements,frequencies)
    OptParams = OptimizeParameters(Objective,Initial_parameters)
    return ObjFunc(OptParams)
end
