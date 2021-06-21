using EquivalentCircuits, DelimitedFiles, Test 
using EquivalentCircuits: subcircuits, flatten, readablecircuit, simplifycircuit!, replace_redundant_cpes!
using EquivalentCircuits: isterminal, isoperation, generatekarva, karva_parameters,karva_to_tree, Circuit
using EquivalentCircuits: tree_to_function, circuitfunction 

data = readdlm(raw"C:\Users\Dell\Documents\EquivalentCircuits.jl\example_measurements.csv",',')
measurements = data[:,1] .+ data[:,2]im
frequencies = data[:,3]
library = loadpopulation("Circuitlibrary.csv");
library_fit = "R1-[C2,R3-[C4,R5]]"

head = 8
encoding = generatekarva(head)
parameters = karva_parameters(encoding)
example_encoding = "+-+RRR-PCLLLCRPPP"
example_parameters = karva_parameters(encoding)
exam_params = [0.0, 0.0, 0.0, 1543.907617339479, 770.7386241387493, 1273.4835727121672, 0.0, [3499.829595663477, 0.8749573989158692],
 1.6880836678204815e-5, 0.04940667507884311, 0.6677953592773065, 0.42043991655286717, 9.192012270021564e-5, 3694.724527430516,
  [618.972522236354, 0.1547431305590885], [764.360425865072, 0.191090106466268], [3817.8583766637335, 0.9544645941659333]]
exam_params_redundant_CPE = [0.0, 0.0, 0.0, 1543.907617339479, 770.7386241387493, 1273.4835727121672, 0.0, [3499.829595663477, 1],
 1.6880836678204815e-5, 0.04940667507884311, 0.6677953592773065, 0.42043991655286717, 9.192012270021564e-5, 3694.724527430516,
  [618.972522236354, 0.1547431305590885], [764.360425865072, 0.191090106466268], [3817.8583766637335, 0.9544645941659333]]
example_tree = karva_to_tree(example_encoding,example_parameters)
example_circuit = Circuit(example_encoding,example_parameters,nothing)
example_circuit_redundant_CPE = EquivalentCircuits.Circuit(example_encoding,exam_params_redundant_CPE,nothing)
example_function = tree_to_function(example_tree)
Readable_example_circuit = readablecircuit(example_circuit)
example_function2 = circuitfunction(Readable_example_circuit)
circuit = EquivalentCircuits.Circuit(encoding,parameters,nothing)
readable_circuit = readablecircuit(circuit)
tree = karva_to_tree(encoding,parameters)


@testset "EquivalentCircuits.jl" begin
    # Check the circuitevolution function.
    library_fit == circuitevolution("example_measurements.csv",initial_population = library)[1]
    # Evaluate properties of the generated circuit encoding.
    EquivalentCircuits.isoperation(encoding[1]) == true
    # Only terminals in the encoding's tail.
    all(EquivalentCircuits.isterminal.(collect(encoding[head+1:end]))) == true
    # conversion of circuit object to user-readable circuit.
    readablecircuit(example_circuit) == "[C1,P2]-R3-[R4,R5]"
    # Check the output of circuit simulation.
    example_function(exam_params,100) == 1789.4845259396134 - 10.86703970004648im
    # Replace redundant CPEs and simplify
    replace_redundant_cpes!(example_circuit_redundant_CPE) # converts the redundant CPE to an equivalent capacitor.
    simplifycircuit!(example_circuit_redundant_CPE) # simplifies all parralelly or serially connected similar components.
    # Circuit simplification should lead to a single resistor element R3.
    readablecircuit(example_circuit_redundant_CPE) == "C1-R2"
    # Parameteroptimisation : basic checks of solution lengths and bounds.
    optparams = parameteroptimisation("[C1,P2]-R3",measurements,frequencies)
    length(optparams) == 3
    size(optparams[2])[1] == 2
    C1,P2_1,P2_2,R3 = flatten(optparams)
    0<C1<10
    0<P2_1<1.0e9
    0<P2_2<1
    0<R3<1.0e9
    # Subtrees length checking. There ought to be as many circuits as coding terminal elements.
    length(subcircuits(example_circuit)) == 3 

end
