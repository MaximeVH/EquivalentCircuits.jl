using EquivalentCircuits
using Test

@testset "EquivalentCircuits.jl" begin
example_circuit = "[R1-R2-[R3-R4,R5]-C6,R7]-R8-[R9-C10,R11]-R12"
example_frequencies = [1 10 50 100 500 1000 5000 10000 50000 100000]
example_generated = Generate_circuit(10,8,"RC")
example_function = CircuitFunction_inclusive_Final(example_generated.circuit)
example_parameters = InitializeParameters_inclusive(example_generated.circuit)
example_measurements = SimulateImpedance(example_function,example_parameters,example_frequencies)


@test example_generated isa LRC_Circuit
@test length(InitializeParameters(example_circuit)) == count(x->x=='R'||x=='C'||x=="L"||x=="P", example_circuit)
@test example_function(example_parameters,example_frequencies[1]) isa Complex
@test CircuitEvolution(example_measurements,example_frequencies) isa Vector{RC_Circuit}

end
