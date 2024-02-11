"""
EquivalentCircuit(circuitstring::String, Parameters::NamedTuple)

A structure used as the output of the evolutionary algorithm. It's fields are 'circuitstring', which is the string representation of a circuit (e.g. "R1-[C2,R3]-P4")
and 'Parameters', a NamedTuple with the circuit's components and their corresponding parameter values.

"""
struct EquivalentCircuit
    circuitstring::String
    Parameters::NamedTuple
end

Base.display(circuit::EquivalentCircuit) = println(circuit.circuitstring)

function parametertuple(circuit,parameters)
    elements = foldl(replace,["["=>"","]"=>"","-"=>"",","=>""],init = denumber_circuit(circuit))
    name_strings = []
    number = 1
    for e in elements
        if e == 'P'
            push!(name_strings,"P$(number)w")
            push!(name_strings,"P$(number)n")
        else
            push!(name_strings,"$(e)$(number)")
        end
        number += 1
    end
    name_symbols = Tuple(map(Symbol, name_strings))
    return NamedTuple{name_symbols}(Tuple(parameters))
end

# function initializecircuit(head=8,terminals="RCLP")
#     karva = generatekarva(head,terminals)
#     parameters = karva_parameters(karva)
#     return Circuit(karva,parameters,nothing)
# end

function initializecircuit(bounds,head=8,terminals="RCLP")
    karva = generatekarva(head,terminals)
    parameters = karva_parameters_(karva,bounds)
    return Circuit(karva,parameters,nothing)
end



function initializepopulation(bounds,size=20,head=8,terminals="RCLP")
    return [initializecircuit(bounds,head,terminals) for i in 1:size]
end

# function initializevariedpopulation(size=30,head=8)
#     subpopulation_size = floor(Int(size/3))
#     RCs = initializepopulation(subpopulation_size,head,"RC")
#     RCLs = initializepopulation(subpopulation_size,head,"RCL")
#     RCLPs = initializepopulation(subpopulation_size,head,"RCLP")
#     return vcat(RCs,RCLs,RCLPs)
# end

function simplifypopulation!(population,terminals="RCPL")
    for circuit in population
        simplifycircuit!(circuit,terminals)
    end
end

function circuit_offspring(circuit_1,circuit_2,terminals = "RCLP")
    offspring = ""
    rand_ = rand()
    if rand_ > 0.5
        offspring = one_point_crossover(circuit_1,circuit_2)
    elseif rand_ > 0.1
        offspring = two_point_crossover(circuit_1,circuit_2)
    else
        offspring = rand([circuit_1,circuit_2])
    end
    return mutate(offspring,terminals)
end

function circuitfitness(circuit,measurements,frequencies,bounds)
    tree = get_circuit_tree(circuit)
    circfunc,params,lower,upper,param_inds = func_and_params_for_optim(tree,bounds)
    objective = objectivefunction_noweight(circfunc,measurements,frequencies)
    params = max.(min.(params,upper),lower.+10^-15)
    optparams,fitness = optimizeparameters(objective,params,lower,upper)
    optparams = max.(min.(optparams,upper),lower.+10^-15) #workaround of optim.jl's bug
    return deflatten_parameters(optparams,tree,param_inds), fitness, param_inds
end

function circuitfitness(circuit,measurements,frequencies)
    tree = get_circuit_tree(circuit)
    circfunc,params,upper,param_inds = func_and_params_for_optim(tree)
    objective = objectivefunction_noweight(circfunc,measurements,frequencies)
    optparams,fitness = optimizeparameters(objective,params,upper)
    optparams = max.(min.(optparams,upper),0) #workaround of optim.jl's bug
    return deflatten_parameters(optparams,tree,param_inds), fitness, param_inds
end

function generate_offspring(population,terminals="RCLP")
    population_size = length(population)
    elite_size = Int(ceil(0.1*population_size))
    mating_pool_size = population_size-elite_size
    progenitors = tournamentselection(population,mating_pool_size,elite_size)
    elites = population[1:elite_size]
    offspring = Array{Circuit}(undef, mating_pool_size)
    for e in 1:mating_pool_size
        offspring[e] = circuit_offspring(progenitors[e],progenitors[mating_pool_size-e+1],terminals)
    end
    vcat(elites,offspring)
end

function get_parameter_bounds(tree,bounds)
    lowers = [bounds[node.Type][1] for node in tree]
    uppers=  [bounds[node.Type][2] for node in tree]
    return lowers,uppers
end

function get_karva_bounds(karva,bounds)
    lowers = [bounds[k][1] for k in karva]
    uppers=  [bounds[k][2] for k in karva]
    return lowers,uppers
end

function evaluate_fitness!(population,measurements,frequencies,bounds)
    params = []
    param_inds = []
    for circuit in population
        params,circuit.fitness,param_inds = @suppress circuitfitness(circuit,measurements,frequencies,bounds)
        circuit.parameters[param_inds] = params
    end
end

function evaluate_fitness!(population,measurements,frequencies)
    params = []
    param_inds = []
    for circuit in population
        params,circuit.fitness,param_inds = circuitfitness(circuit,measurements,frequencies)
        circuit.parameters[param_inds] = params
    end
end

function population_from_string(circuitstring_list, head=8, size=30, terminals="RCLP")
    if typeof(circuitstring_list) == Vector{Tuple{String, Vector{Float64}}}
        circuit_pool = [circuit_to_karva(ci[1],ci[2],head=head,terminals=terminals) for ci in circuitstring_list]
    elseif typeof(circuitstring_list) == Vector{String}
    circuit_pool = [circuit_to_karva(ci,head=head,terminals=terminals) for ci in circuitstring_list]
    end
    circuitpopulation = Vector{Circuit}(undef,size)
    for c in 1:size
        if c <= length(circuit_pool)
            circuitpopulation[c] = circuit_pool[c]
        else
            circuitpopulation[c] = mutate(rand(circuit_pool),terminals)
        end
    end
    return circuitpopulation
end

"""
circuit_evolution(measurements::Array{Complex{Float64},1},frequencies::Array{Float64,1}; <keyword arguments>)

Identify an equivalent electrical circuit that fits a given set of electrochemical impedance spectroscopy measurements.

The inputs are a circuit (e.g. "R1-[C2,R3]-P4") an array of complex-valued impedance measurements and their corresponding
frequencies. The output is an EquivalentCircuit object containing a field circuitstring, which is a string denoting the identified circuit
and a field Parameters, which is a NamedTuple of the circuit's components with their corresponding parameter values.

# Arguments
- `generations::Integer=10`: the maximum number of iterations of the evolutionary algorithm.
- `population_size::Integer=30`: the number of individuals in the population during each iteration.
- `terminals::String="RCLP"`: the circuit components that are to be included in the circuit identification.
- `head::Integer=8`: a hyperparameter than controls the maximum considered complexity of the circuits.
- `cutoff::Float64=0.8`: a hyperparameter that controls the circuit complexity by removing redundant components.
- `initial_population::Array{Circuit,1}=nothing`:the option to provide an initial population of circuits
(obtained by using the loadpopulation function) with which the algorithm starts.
 Alternatively, users can provide a custom list of circuits which can either be a list of one or more circuit strings or a list of tuples
 where each tuple has the circuit string as first value and the parameters as second value.
 - `convergence_threshold::Float64=5e-4`: Convergence threshold for circuit identification. Increase this to reduce strictness, e.g., in case of noisy measurements.
- `bounds`::Dict{Char, Vector}: Optional custom bounds for the circuit component parameter values.
# Example
```julia

julia> using EquivalentCircuits, Random

julia> Random.seed!(25);

julia> measurements = [5919.9 - 15.7, 5918.1 - 67.5im, 5887.1 - 285.7im, 5428.9 - 997.1im, 3871.8 - 978.9im,
3442.9 - 315.4im, 3405.5 - 242.5im, 3249.6 - 742.0im, 1779.4 - 1698.9im,  208.2 - 777.6im, 65.9 - 392.5im];

julia> frequencies = [0.10, 0.43, 1.83, 7.85, 33.60, 143.84, 615.85,  2636.65, 11288.38, 48329.30, 100000.00];

julia> circuit_evolution(measurements, frequencies , generations= 15, terminals = "RC")
[R1,C2]-[C3,R4]-R5
```
"""
function circuit_evolution(measurements,frequencies;generations::Real=10,population_size=30,terminals = "RCLP",head=8,cutoff=0.8,initial_population=nothing,convergence_threshold=5e-4,bounds=nothing,quiet=false)
    quiet && disable_logging(Logging.Warn)
    parameter_bounds = Dict('R'=>[0,1.0e9],'C'=>[0,10],'L'=>[0,5],'P'=>[[0,0],[1.0e9,1]],'W'=>[0,1.0e9],'+'=>[0,0],'-'=>[0,0])
    if typeof(bounds) == Dict{Char, Vector}
        for key in keys(bounds)
            parameter_bounds[key] = bounds[key]
        end
    end
    # Either initialize a new population, or work with a provided initial population.
    # @assert typeof(initial_population) in [nothing,Vector{Circuit},Vector{String},Vector{Tuple{String, Vector{Float64}}}] "The initial population must be a list."
    if isnothing(initial_population)
        population = initializepopulation(parameter_bounds,population_size,head,terminals) #initializevariedpopulation(population_size,head)
    elseif typeof(initial_population) == Vector{Circuit}
        population = initial_population
    elseif typeof(initial_population) in [Vector{String},Vector{Tuple{String, Vector{Float64}}}]
        population = population_from_string(initial_population, head, population_size, terminals)
    end
    # Theoretical simplification of the initial population.
        simplifypopulation!(population,terminals)
        evaluate_fitness!(population,measurements,frequencies,parameter_bounds)
        sort!(population)
        generation = 0
        # Keep track of the fittest individual circuit.
        elite = minimum(population)
        min_fitness = elite.fitness
        while (min_fitness > convergence_threshold) && (generation<=generations)
            population = generate_offspring(population,terminals)
            simplifypopulation!(population,terminals)
            evaluate_fitness!(population,measurements,frequencies,parameter_bounds)
            sort!(population)
            elite = minimum(population).fitness < elite.fitness ? minimum(population) : elite
            min_fitness = elite.fitness
            population[1] = redundacy_testing(population[1],measurements,frequencies,parameter_bounds,terminals,cutoff)
            generation += 1
        end
        replace_redundant_cpes!(population)
        population = removeduplicates(sort!(vcat(population,elite)))
        for i in 1:3
            population[i] = removeredundancy(population[i],measurements,frequencies,parameter_bounds,terminals,cutoff)
        end
        # extract the converged circuits.
    population = filter(p -> p.fitness ≤ convergence_threshold, population)
    # in case of no converged circuits => alternate output print statement "Algorithm did not converge"
    if length(population) == 0
        @info "Algorithm did not converge"
    else
        best_circuit = readablecircuit(population[1])
        return EquivalentCircuit(best_circuit,parameteroptimisation(best_circuit,measurements,frequencies)) #readablecircuit.(population[1:top_n])
    end
end

"""
    circuit_evolution(filepath::String; <keyword arguments>)

 Identify an equivalent electrical circuit that fits a given set of electrochemical impedance spectroscopy measurements.

 The inputs are a circuit (e.g. "R1-[C2,R3]-P4") and a filepath to a CSV file containing the three following columns:
 the real part of the impedance, the imaginary part of the impedance, and the frequencies corresponding to the measurements.
 The output is NamedTuple of the circuit's components with their corresponding parameter values.

 # Arguments
- `filepath::String`: filepath of the CSV file containing the impedance measurements.
- `generations::Integer=10`: the maximum number of iterations of the evolutionary algorithm.
- `population_size::Integer=30`: the number of individuals in the population during each iteration.
- `terminals::String="RCLP"`: the circuit components that are to be included in the circuit identification.
- `head::Integer=8`: a hyperparameter than controls the maximum considered complexity of the circuits.
- `cutoff::Float64=0.8`: a hyperparameter that controls the circuit complexity by removing redundant components.
- `initial_population::Array{Circuit,1}=nothing`:the option to provide an initial population of circuits
(obtained by using the loadpopulation function) with which the algorithm starts.
 Alternatively, users can provide a custom list of circuits which can either be a list of one or more circuit strings or a list of tuples
 where each tuple has the circuit string as first value and the parameters as second value.
 - `convergence_threshold::Float64=5e-4`: Convergence threshold for circuit identification. Increase this to reduce strictness, e.g., in case of noisy measurements.
 - `bounds`::Dict{Char, Vector}: Optional custom bounds for the circuit component parameter values.

"""
function circuit_evolution(filepath::String;generations::Real=10,population_size=30,terminals = "RCLP",head=8,cutoff=0.8,initial_population = nothing,convergence_threshold=5e-4)
    # Read the measurement file.
    meansurement_file = readdlm(filepath,',')
    # convert the measurement data into usable format.
    reals = meansurement_file[:,1]
    imags = meansurement_file[:,2]
    frequencies = meansurement_file[:,3]
    measurements = reals + imags*im
    parameter_bounds = Dict('R'=>[0,1.0e9],'C'=>[0,10],'L'=>[0,5],'P'=>[[0,0],[1.0e9,1]],'W'=>[0,1.0e9],'+'=>[0,0],'-'=>[0,0])
    if typeof(bounds) == Dict{Char, Vector}
        for key in keys(bounds)
            parameter_bounds[key] = bounds[key]
        end
    end
    # Either initialize a new population, or work with a provided initial population.
    # @assert typeof(initial_population) in [nothing,Vector{Circuit},Vector{String},Vector{Tuple{String, Vector{Float64}}}] "The initial population must be a list."
    if isnothing(initial_population)
        population = initializepopulation(parameter_bounds,population_size,head,terminals) #initializevariedpopulation(population_size,head)
    elseif typeof(initial_population) == Vector{Circuit}
        population = initial_population
    elseif typeof(initial_population) in [Vector{String},Vector{Tuple{String, Vector{Float64}}}]
        population = population_from_string(initial_population, head, population_size, terminals)
    end
    # Theoretical simplification of the initial population.
        simplifypopulation!(population,terminals)
        evaluate_fitness!(population,measurements,frequencies,parameter_bounds)
        sort!(population)
        generation = 0
        # Keep track of the fittest individual circuit.
        elite = minimum(population)
        min_fitness = elite.fitness
        while (min_fitness > convergence_threshold) && (generation<=generations)
            population = generate_offspring(population,terminals)
            simplifypopulation!(population,terminals)
            evaluate_fitness!(population,measurements,frequencies,parameter_bounds)
            sort!(population)
            elite = minimum(population).fitness < elite.fitness ? minimum(population) : elite
            min_fitness = elite.fitness
            population[1] = redundacy_testing(population[1],measurements,frequencies,parameter_bounds,terminals,cutoff)
            generation += 1
        end
        replace_redundant_cpes!(population)
        population = removeduplicates(sort!(vcat(population,elite)))
        for i in 1:3
            population[i] = removeredundancy(population[i],measurements,frequencies,parameter_bounds,terminals,cutoff)
        end
        # extract the converged circuits.
    population = filter(p -> p.fitness ≤ convergence_threshold, population)
    # in case of no converged circuits => alternate output print statement "Algorithm did not converge"
    if length(population) == 0
        @info "Algorithm did not converge"
    else
        best_circuit = readablecircuit(population[1])
        return EquivalentCircuit(best_circuit,parameteroptimisation(best_circuit,measurements,frequencies)) #readablecircuit.(population[1:top_n])
    end
end


"""
    circuit_evolution_batch(args...; kwargs...)

Similar to the `circuit_evolution` function, but to generate multiple candidate circuits.

# Arguments

- `measurements::Array{Complex{Float64},1}`: Array of complex-valued impedance measurements.
- `frequencies::Array{Float64,1}`: Frequencies corresponding to the measurements.

# Keywords

- `generations::Integer=10`: Maximum number of iterations of the evolutionary algorithm.
- `population_size::Integer=30`: Number of individuals in the population during each iteration.
- `terminals::String="RCLP"`: Allowed components to be included in the circuit.
- `head::Integer=8`: Controls the maximum considered complexity of the circuits.
- `cutoff::Float64=0.8`: Controls the circuit complexity by removing redundant components.
- `initial_population::Array{Circuit,1}=nothing`: The option to provide an initial
    population of circuits (obtained by using the loadpopulation function) with which the
    algorithm starts. Alternatively, users can provide a custom list of circuits which can
    either be a list of one or more circuit strings or a list of tuples where each tuple
    has the circuit string as first value and the parameters as second value.
- `convergence_threshold::Float64=5e-4`: Convergence threshold for circuit identification.
    Increase to reduce strictness, e.g., in case of noisy measurements.
- `bounds`::Dict{Char, Vector}: Upper/lower bounds for the circuit component parameter values.
- `numprocs::Integer=num_physical_cores()`: Number of processes to use for parallelisation.
- `iters`: Number of candidate circuits to generate.

# Returns

- `results::Vector{EquivalentCircuit}`: A vector of EquivalentCircuit objects.

"""
function circuit_evolution_batch(
    measurements,
    frequencies;
    generations::Real=10,
    population_size=30,
    terminals="RCLP",
    head=8,
    cutoff=0.8,
    initial_population=nothing,
    convergence_threshold=5e-4,
    bounds=nothing,
    numprocs=num_physical_cores(),
    iters,
)
    # Set up workers
    _workers = addprocs(min(numprocs, iters))
    import_module_on_workers(_workers)

    # Run the circuit evolution in parallel on all workers
    @info "Starting circuit evolution on $(length(_workers)) workers"
    results = pmap(
        _ -> circuit_evolution(
            measurements,
            frequencies;
            generations=generations,
            population_size=population_size,
            terminals=terminals,
            head=head,
            cutoff=cutoff,
            initial_population=initial_population,
            convergence_threshold=convergence_threshold,
            bounds=bounds,
        ),
        WorkerPool(_workers),
        1:iters,
    )

    # Remove workers
    rmprocs(_workers)

    return results
end
