using EquivalentCircuits
using Test

data = readdlm(raw"C:\Users\Dell\Documents\EquivalentCircuits.jl\example_measurements.csv",',')
measurements = data[:,1] .+ data[:,2]im
frequencies = data[:,3]
library = loadpopulation("Circuitlibrary.csv");
library_fit =["R1-[C2,R3-[C4,R5]]"
"P1-L2-[R3,C4]-[R5,C6]"
"[R1,C2-[R3,C4]]"
"[R1-[P2,R3],C4]"
"R1-[C2,R3-P4]"]

@testset "EquivalentCircuits.jl" begin
    library_fit == circuitevolution("example_measurements.csv",initial_population = library)
end