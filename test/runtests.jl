using EquivalentCircuits, DelimitedFiles, Test 
using EquivalentCircuits: subcircuits, flatten, readablecircuit, simplifycircuit!, replace_redundant_cpes!
using EquivalentCircuits: isterminal, isoperation, generatekarva, karva_parameters,karva_to_tree, Circuit
using EquivalentCircuits: tree_to_function, circuitfunction 

measurements = [5919.90 - 15.79im, 5919.58 - 32.68im, 5918.18 - 67.58im, 5912.24 - 139.49im,
 5887.12 - 285.74im, 5785.04 - 566.88im, 5428.94 - 997.19im, 4640.21 - 1257.83im, 3871.84 - 978.97im,
  3537.68 - 564.96im, 3442.94 - 315.40im, 3418.14 - 219.69im, 3405.51 - 242.57im, 3373.90 - 396.07im,
   3249.67 - 742.03im, 2808.42 - 1305.92im, 1779.41 - 1698.97im, 701.96 - 1361.47im, 208.29 - 777.65im, 65.93 - 392.51im]

frequencies = [0.10, 0.21, 0.43, 0.89, 1.83, 3.79, 7.85, 16.24, 33.60, 69.52, 143.84, 297.64, 615.85, 1274.27, 2636.65,
 5455.59, 11288.38, 23357.21, 48329.30, 100000.00]

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
    @test library_fit == circuitevolution(measurements,frequencies,initial_population = library)[1]
    # Evaluate properties of the generated circuit encoding.
    @test EquivalentCircuits.isoperation(encoding[1]) == true
    # Only terminals in the encoding's tail.
    @test all(EquivalentCircuits.isterminal.(collect(encoding[head+1:end]))) == true
    # conversion of circuit object to user-readable circuit.
    @test readablecircuit(example_circuit) == "[C1,P2]-R3-[R4,R5]"
    # Check the output of circuit simulation.
    @test example_function(exam_params,100) == 1789.4845259396134 - 10.86703970004648im
    # Replace redundant CPEs and simplify
    replace_redundant_cpes!(example_circuit_redundant_CPE) # converts the redundant CPE to an equivalent capacitor.
    simplifycircuit!(example_circuit_redundant_CPE) # simplifies all parralelly or serially connected similar components.
    # Circuit simplification should lead to a single resistor element R3.
    @test readablecircuit(example_circuit_redundant_CPE) == "C1-R2"
    # Parameteroptimisation : basic checks of solution lengths and bounds.
    optparams = parameteroptimisation("[C1,P2]-R3",measurements,frequencies)
    @test length(optparams) == 3
    @test size(optparams[2])[1] == 2
    C1,P2_1,P2_2,R3 = flatten(optparams)
    @test 0<C1<10
    @test 0<P2_1<1.0e9
    @test 0<P2_2<1
    @test 0<R3<1.0e9
    # Subtrees length checking. There ought to be as many circuits as coding terminal elements.
    @test length(subcircuits(example_circuit)) == 5 

end