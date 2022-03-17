using EquivalentCircuits, DelimitedFiles, Test 
using EquivalentCircuits: subcircuits, flatten, readablecircuit, simplifycircuit!, replace_redundant_cpes!
using EquivalentCircuits: isterminal, isoperation, generatekarva, karva_parameters,karva_to_tree, Circuit
using EquivalentCircuits: tree_to_function, circuitfunction 

measurements = [5919.90 - 15.79im, 5919.58 - 32.68im, 5918.18 - 67.58im, 5912.24 - 139.49im,
 5887.12 - 285.74im, 5785.04 - 566.88im, 5428.94 - 997.19im, 4640.21 - 1257.83im, 3871.84 - 978.97im,
  3537.68 - 564.96im, 3442.94 - 315.40im, 3418.14 - 219.69im, 3405.51 - 242.57im, 3373.90 - 396.07im,
   3249.67 - 742.03im, 2808.42 - 1305.92im, 1779.41 - 1698.97im, 701.96 - 1361.47im, 208.29 - 777.65im, 65.93 - 392.51im];

frequencies = [0.10, 0.21, 0.43, 0.89, 1.83, 3.79, 7.85, 16.24, 33.60, 69.52, 143.84, 297.64, 615.85, 1274.27, 2636.65,
 5455.59, 11288.38, 23357.21, 48329.30, 100000.00];

library = [Circuit("+-R+C-RRCLCRCLLLR", Any[0, 0, 22.377192750081925, 0, 4.042529368461352e-9, 0, 3399.0384392519263, 2498.7384040467605, 4.0061261087766165e-6, 0.7948425599214597, 6.126315997167497e-5, 1373.5802658510875, 3.6191586501947424e-5, 0.17338854108986923, 0.46734655676131776, 0.24812830993871104, 1500.925039002129], 8.293095935567338e-6)
Circuit("+++--+RCRPRLRRCLL", Any[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 16.952525176348484, 4.019951476143508e-9, 3379.008162497164, [209604.83569111855, 0.9610007184206273], 2549.666570971039, 3.2866976510279723e-7, 2.418161729014113, 1654.6779402286252, 7.295488055251058e-5, 0.2994773175317633, 0.2905678390917996], 1.7509954992998087e-5)
Circuit("+-RRPCCRLLPPLLRCP", Any[0, 0, 0.00041729614400979383, 4595.3298212413, [6.1253600051612616e7, 0.8859174267613648], 8.481342404807135e-6, 8.005423544288912e-5, 1537.5259065364694, 0.7866046048091446, 0.4825771386985225, [3075.672528367832, 0.768918132091958], [421.20526643818135, 0.10530131660954534], 0.23435266803982535, 0.4782133397334183, 657.0140358658155, 5.755479554270154e-5, [524.214144301653, 0.13105353607541326]], 0.02705946970454961)];

library_fit = "R1-[C2,R3-[C4,R5]]";

head = 8;
encoding = generatekarva(head);
parameters = karva_parameters(encoding);
example_encoding = "+-+RRR-PCLLLCRPPP";
example_parameters = karva_parameters(encoding);
exam_params = [0.0, 0.0, 0.0, 1543.907617339479, 770.7386241387493, 1273.4835727121672, 0.0, [3499.829595663477, 0.8749573989158692],
 1.6880836678204815e-5, 0.04940667507884311, 0.6677953592773065, 0.42043991655286717, 9.192012270021564e-5, 3694.724527430516,
  [618.972522236354, 0.1547431305590885], [764.360425865072, 0.191090106466268], [3817.8583766637335, 0.9544645941659333]];
exam_params_redundant_CPE = [0.0, 0.0, 0.0, 1543.907617339479, 770.7386241387493, 1273.4835727121672, 0.0, [3499.829595663477, 1],
 1.6880836678204815e-5, 0.04940667507884311, 0.6677953592773065, 0.42043991655286717, 9.192012270021564e-5, 3694.724527430516,
  [618.972522236354, 0.1547431305590885], [764.360425865072, 0.191090106466268], [3817.8583766637335, 0.9544645941659333]];
example_tree = karva_to_tree(example_encoding,example_parameters);
example_circuit = Circuit(example_encoding,example_parameters,nothing);
example_circuit_redundant_CPE = Circuit(example_encoding,exam_params_redundant_CPE,nothing);
example_function = tree_to_function(example_tree);
Readable_example_circuit = readablecircuit(example_circuit);
example_function2 = circuitfunction(Readable_example_circuit);
circuit = Circuit(encoding,parameters,nothing);
readable_circuit = readablecircuit(circuit);
tree = karva_to_tree(encoding,parameters);


@testset "EquivalentCircuits.jl" begin
    # Check the circuit_evolution function.
    @test library_fit == circuit_evolution(measurements,frequencies,initial_population = library).Circuit
    # Evaluate properties of the generated circuit encoding.
    @test isoperation(encoding[1]) == true
    # Only terminals in the encoding's tail.
    @test all(isterminal.(collect(encoding[head+1:end]))) == true
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
    @test length(optparams) == 4
    C1,P2_1,P2_2,R3 = optparams
    @test 0<C1<10
    @test 0<P2_1<1.0e9
    @test 0<P2_2<1
    @test 0<R3<1.0e9
    # Subtrees length checking. There ought to be as many circuits as coding terminal elements.
    @test length(subcircuits(example_circuit)) == 5 

end