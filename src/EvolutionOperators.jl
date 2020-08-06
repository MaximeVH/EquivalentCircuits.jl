"""
    CircuitCrossover(p1,p2,max_depth=10)

Recombine two parent circuits, p1 and p2, by replacing a randomly selected
component or subcircuit of p1 by a randomly selected subcircuit (that can also be nothing)
of p2. The output is a new RC_Circuit object.
"""
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
"""
    CircuitMutation(p1,mutation_probability=0.2)

A given circuit p1 is subjected to mutation with a probability of mutation_probability.
mutation entails the generation of a new component/subcircuit that is randomly inserted in the
parent circuit.
"""
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
"""
    CircuitFitness(circuit,measurements,frequencies)

Outputs the the residual least squares objective between the
impedance measurements and the simulated impedance values of a given circuit with the input frequencies of the measurements.
the parameters of the circuit are first optimised with respect to the measurement values.
The output fitness value is a measure of how well a given circuit configuration is capable
of fitting the experimental data.
"""
function CircuitFitness(circuit,measurements,frequencies)
    CircFunc = CircuitFunction2(circuit.string)
    ObjFunc = ObjectiveFunction(CircFunc,measurements,frequencies)
    OptParams = OptimizeParameters(ObjFunc,circuit.parameters)
    return ObjFunc(OptParams)
end
