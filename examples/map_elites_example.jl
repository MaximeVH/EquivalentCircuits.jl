using Plots

# Example measurements and corresponding frequencies.
measurements = [5919.9 - 15.7, 5918.1 - 67.5im, 5887.1 - 285.7im, 5428.9 - 997.1im, 3871.8 - 978.9im, 
3442.9 - 315.4im, 3405.5 - 242.5im, 3249.6 - 742.0im, 1779.4 - 1698.9im,  208.2 - 777.6im, 65.9 - 392.5im];
frequencies = [0.10, 0.43, 1.83, 7.85, 33.60, 143.84, 615.85,  2636.65, 11288.38, 48329.30, 100000.00];

# P is the matrix that keeps track of the circuit archive's fitness for different circuit complexities.
# X is the archive itself, where position (i,j), if occupied, is the circuit containing i integer components and j fractional components.
# I is the number of iterations (usually low) required to find a circuit under the convergence threshold.
@time P,X,I = circuit_map_elites_R(measurements,frequencies,population_size=40,iterations=60)

# Two functions to visualize the found circuits according to their complexity and fitness
annotated_archive(P,X)
plot_MAP_archive(P)

# Gets the population of circuits
population_ = population_from_archive(P,X)

# Nyquist plots of the well-performing candidate cirucits, the cutoff argument is defaultly set at 1e-3
Archive_nyquist(measurements,frequencies,P,X)

# Continue searching for more circuits over N generations (recommended).
N = 50
@time P,X = circuit_map_elites_continue_R(N,P,X,measurements,frequencies)

# Get the resulting population of circuit instances as an array.
population_ = population_from_archive(P,X)

# Adjusted version of MAP-Elites with multiple solutions (pe argument) at each feature description. 
# This version takes longer but provides more solution candidates.
@time P3,X3,I3 = circuit_map_elites_R_m(measurements,frequencies,population_size=40,iterations=60,pe=3)
@time P3,X3 = circuit_map_elites_continue_R_m(50,P3,X3,measurements,frequencies)

# The same plotting functions are supported for the resulting outputs.
annotated_archive(P3,X3)
plot_MAP_archive(P3)
Archive_nyquist(measurements,frequencies,P3,X3)
