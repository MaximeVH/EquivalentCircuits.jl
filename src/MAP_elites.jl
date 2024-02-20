## See examples/map_elites_example.jl for usage examples
## TO DO: more code documentation
## Version with serially connected resistor below

archive_to_population(P,χ) = sort(χ[P .> 0])

function random_selection(χ,P)
    return selection = χ[rand(findall(!iszero,P))]
end

function feature_descriptor(x′::Circuit)
    coding_karva = get_coding_karva(x′)
    component_count = count(isterminal,coding_karva)
    operation_count = count(isoperation,coding_karva)
    return CartesianIndex((component_count, operation_count))
end

isinteger_order(chr) = occursin(chr,"RCL")
isfractional_order(chr) = occursin(chr,"PW")

function fractional_descriptor(x′::Circuit)
    coding_karva = get_coding_karva(x′)
    integer_order_count = count(isinteger_order,coding_karva)
    fractional_order_count = count(isfractional_order,coding_karva)
    return CartesianIndex((integer_order_count+1, fractional_order_count+1))
end

function get_coding_karva(x′::Circuit)
    return join([c.Type for c in karva_to_tree(x′.karva)])
end

function evaluate_circuit_fitness!(circuit,measurements,frequencies,bounds)
    params,circuit.fitness,param_inds = circuitfitness(circuit,measurements,frequencies,bounds);
    circuit.parameters[param_inds] = params;
end

function evaluate_population_fitness!(population,measurements,frequencies,bounds)
    for circuit in population
    params,circuit.fitness,param_inds = circuitfitness(circuit,measurements,frequencies,bounds);
    circuit.parameters[param_inds] = params;
    end
end

function coding_length(x′::Circuit)
    return length(get_coding_karva(x′))
end

function mutate_circuit(circuit::Circuit,terminals = "RCLP")
    karva = circuit.karva ; parameters = copy(circuit.parameters)
    total_length = length(karva)
    head_length = (total_length-1)/2 ; tail_length = head_length+1
    operations = "+-"; new_element = ""
    position = rand(1:coding_length(circuit))
    mutable_karva = collect(karva)
    if position == 1
        new_element = rand(operations)
        mutable_karva[position] = new_element
    elseif position <= head_length
        new_element = rand(operations*terminals)
        mutable_karva[position] = new_element
        parameters[position] = replacementparameter(new_element)
    elseif position > head_length
        new_element = rand(terminals)
        mutable_karva[position] = new_element
        parameters[position] = replacementparameter(new_element)
    end
    return Circuit(String(mutable_karva),parameters,nothing)
end

## MAP elites with whole population consideration per generation
function circuit_map_elites(measurements,frequencies;population_size= 30,iterations=10,feat_desc=fractional_descriptor,convergence_threshold = 5e-5,bounds = nothing,terminals="RCLP",head=8,initial_population=nothing)
    parameter_bounds = Dict('R'=>[0,1.0e9],'C'=>[0,10],'L'=>[0,5],'P'=>[[0,0],[1.0e9,1]],'W'=>[0,1.0e9],'+'=>[0,0],'-'=>[0,0])
    if typeof(bounds) == Dict{Char, Vector}
        for key in keys(bounds)
            parameter_bounds[key] = bounds[key]
        end
    end
    if isnothing(initial_population)  
        population = initializepopulation(parameter_bounds,population_size,head,terminals) 
    elseif typeof(initial_population) == Vector{Circuit}
        population = initial_population
    elseif typeof(initial_population) in [Vector{String},Vector{Tuple{String, Vector{Float64}}}]
        population = population_from_string(initial_population, head, population_size, terminals)
    end
    simplifypopulation!(population,terminals)
    @suppress evaluate_population_fitness!(population,measurements,frequencies,parameter_bounds)
    init_descriptions = [feat_desc(circ) for circ in population]
    init_fitnesses = [circ.fitness for circ in population]
    min_fit, iter = 1,0;
    P = fill(0.0,(12,12));
    χ = Matrix{Circuit}(undef,(12,12));
    for (desc,fit,circ) in zip(init_descriptions,init_fitnesses,population) 
        if iszero(P[desc]) || P[desc] > fit
            P[desc] = fit
            χ[desc] = circ
        end
    end
    while (iter ≤ iterations) && (min_fit ≥ convergence_threshold)
        iter += 1
        for i in 1:population_size
            x = random_selection(χ,P)
            x′ = simplifycircuit(mutate_circuit(x,terminals))
            b′ = feat_desc(x′) #als input.
            @suppress evaluate_circuit_fitness!(x′,measurements,frequencies,parameter_bounds)
            p′ = x′.fitness
            if iszero(P[b′]) || P[b′] > p′
                P[b′] = p′
                χ[b′] = x′
            end
        end
        min_fit = minimum(P[P .> 0])
    end
    return P,χ,iter
end

function add_to_archive(circuit,P,χ,measurements,frequencies,bounds=nothing,feat_desc=fractional_descriptor)
    parameter_bounds = Dict('R'=>[0,1.0e9],'C'=>[0,10],'L'=>[0,5],'P'=>[[0,0],[1.0e9,1]],'W'=>[0,1.0e9],'+'=>[0,0],'-'=>[0,0])
    if typeof(bounds) == Dict{Char, Vector}
        for key in keys(bounds)
            parameter_bounds[key] = bounds[key]
        end
    end
    x′ = circuit
    b′ = feat_desc(x′) #als input.
    @suppress evaluate_circuit_fitness!(x′,measurements,frequencies,parameter_bounds)
    p′ = x′.fitness
    if iszero(P[b′]) || P[b′] > p′
        P[b′] = p′
        χ[b′] = x′
    end
    return P, χ
end

function add_to_archive_multiple(circuitlist,P,χ,measurements,frequencies,bounds=nothing,feat_desc=fractional_descriptor)
    parameter_bounds = Dict('R'=>[0,1.0e9],'C'=>[0,10],'L'=>[0,5],'P'=>[[0,0],[1.0e9,1]],'W'=>[0,1.0e9],'+'=>[0,0],'-'=>[0,0])
    if typeof(bounds) == Dict{Char, Vector}
        for key in keys(bounds)
            parameter_bounds[key] = bounds[key]
        end
    end
    for circuit in circuitlist
    x′ = circuit
    b′ = feat_desc(x′) #als input.
    @suppress evaluate_circuit_fitness!(x′,measurements,frequencies,parameter_bounds)
    p′ = x′.fitness
    if iszero(P[b′]) || P[b′] > p′
        P[b′] = p′
        χ[b′] = x′
    end
end
    return P, χ
end

function add_serial_resistor!(circuit::Circuit,resistance=10)
    new_karva = "+R" * circuit.karva 
    pushfirst!(circuit.parameters,resistance)
    pushfirst!(circuit.parameters,0)
    circuit.karva = new_karva
end

## MAP elites including recombination, mutation and recombination probabilities are also optional input arguments.
function circuit_map_elites_reco(measurements,frequencies;population_size= 30,iterations=10,feat_desc=fractional_descriptor,convergence_threshold = 5e-5,bounds = nothing,terminals="RCLP",head=8,initial_population=nothing)
    parameter_bounds = Dict('R'=>[0,1.0e9],'C'=>[0,10],'L'=>[0,5],'P'=>[[0,0],[1.0e9,1]],'W'=>[0,1.0e9],'+'=>[0,0],'-'=>[0,0])
    if typeof(bounds) == Dict{Char, Vector}
        for key in keys(bounds)
            parameter_bounds[key] = bounds[key]
        end
    end
    if isnothing(initial_population)  
        population = initializepopulation(parameter_bounds,population_size,head,terminals) 
    elseif typeof(initial_population) == Vector{Circuit}
        population = initial_population
    elseif typeof(initial_population) in [Vector{String},Vector{Tuple{String, Vector{Float64}}}]
        population = population_from_string(initial_population, head, population_size, terminals)
    end
    simplifypopulation!(population,terminals)
    @suppress evaluate_population_fitness!(population,measurements,frequencies,parameter_bounds)
    init_descriptions = [feat_desc(circ) for circ in population]
    init_fitnesses = [circ.fitness for circ in population]
    min_fit, iter = 1,0;
    P = fill(0.0,(12,12));
    χ = Matrix{Circuit}(undef,(12,12));
    for (desc,fit,circ) in zip(init_descriptions,init_fitnesses,population) 
        if iszero(P[desc]) || P[desc] > fit
            P[desc] = fit
            χ[desc] = circ
        end
    end
    while (iter ≤ iterations) && (min_fit ≥ convergence_threshold)
        iter += 1
        for e in 1:population_size
            x′ = simplifycircuit(circuit_offspring(random_selection(χ,P),random_selection(χ,P),terminals))
            b′ = feat_desc(x′) #als input.
                @suppress evaluate_circuit_fitness!(x′,measurements,frequencies,parameter_bounds)
                p′ = x′.fitness
                if iszero(P[b′]) || P[b′] > p′
                    P[b′] = p′
                    χ[b′] = x′
                end
        end
        min_fit = minimum(P[P .> 0])
    end
    return P,χ,iter
end

## Allow fixed subcircuit part to be serially connected.
function population_from_archive(P,χ)
   return χ[findall(!iszero,P)]
end

function plot_MAP_archive(p::Array{Float64, 2})
    X = []
    Y = []
    Z = []
    for x in 1:size(p)[1]
        for y in 1:size(p)[2]
            push!(X,x-1)
            push!(Y,y-1)
            push!(Z,p[x,y])
        end
    end
    non_zero_inds = findall(!iszero,Z)
    Z_ = Z[non_zero_inds]
    X_ = X[non_zero_inds]
    Y_ = Y[non_zero_inds]
    fig = scatter(X_,Y_,log10.(Z_),legend=false,color=:green)
    xlabel!("Integer-order count")
    ylabel!("Fractional-order count")
    zlabel!("Log fitness")
    title!("Feature-performance plot")
    return fig
end

function annotated_archive(p,x) #find a way to avoid the overlay plotting, maybe with sortperm
    non_zero_individuals = findall(!iszero,p)
    Z = p[non_zero_individuals]
    X = [non_zero_individuals[i][1] for i in eachindex(non_zero_individuals)]
    Y = [non_zero_individuals[i][2] for i in eachindex(non_zero_individuals)]
    new_x = X .+ Y .-2
    a = scatter(new_x ,log10.(Z),markersize =1.5,markerstrokewidth=0,color=:blue,legend=false)
    for i in eachindex(new_x)
        annotate!(new_x[i],log10(Z[i])+0.1,text(readablecircuit(x[non_zero_individuals][i]),4));
    end
    xlims!(1.5,10)
    xlabel!("Circuit complexity")
    ylabel!("log10 Fitness")
        return a
end

function plot_MAP_archive_next(p::Array{Float64, 2})
    X = []
    Y = []
    Z = []
    for x in 1:size(p)[1]
        for y in 1:size(p)[2]
            push!(X,x-1)
            push!(Y,y-1)
            push!(Z,p[x,y])
        end
    end
    non_zero_inds = findall(!iszero,Z)
    Z_ = Z[non_zero_inds]
    X_ = X[non_zero_inds]
    Y_ = Y[non_zero_inds]
    fig = scatter(X_,Y_,log10.(Z_),legend=false,color=:green)
    xlabel!("Integer-order count")
    ylabel!("Fractional-order count")
    title!("Feature-performance plot")
    return fig
end


function plot_MAP_archive(p::Array{Float64, 2}, cam = (40, 30))
    X = []
    Y = []
    Z = []
    for x in 1:size(p)[1]
        for y in 1:size(p)[2]
            push!(X,x)
            push!(Y,y)
            push!(Z,p[x,y])
        end
    end
    Z_ = Z[findall(!iszero,Z)]
    X_ = X[findall(!iszero,Z)]
    Y_ = Y[findall(!iszero,Z)]
    fig = scatter(X_,Y_,log10.(Z_),legend=false,camera = cam,color=:green)  
    xlabel!("Integer-order count")
    ylabel!("Fractional-order count")
    zlabel!("log10 fitness")
    title!("Feature-performance plot")
    return fig
end


####################### Version of MAP-Elites with serial resistors #######################

function circuit_map_elites_R(measurements,frequencies;population_size= 30,iterations=10,feat_desc=fractional_descriptor,convergence_threshold = 5e-5,bounds = nothing,terminals="RCLP",head=8,initial_population=nothing)
    parameter_bounds = Dict('R'=>[0,1.0e9],'C'=>[0,10],'L'=>[0,5],'P'=>[[0,0],[1.0e9,1]],'W'=>[0,1.0e9],'+'=>[0,0],'-'=>[0,0])
    if typeof(bounds) == Dict{Char, Vector}
        for key in keys(bounds)
            parameter_bounds[key] = bounds[key]
        end
    end
    if isnothing(initial_population)  
        population = initializepopulation_R(parameter_bounds,population_size,head,terminals) 
    elseif typeof(initial_population) == Vector{Circuit}
        population = initial_population
    elseif typeof(initial_population) in [Vector{String},Vector{Tuple{String, Vector{Float64}}}]
        population = population_from_string(initial_population, head, population_size, terminals)
    end
    simplifypopulation!(population,terminals)
    @suppress evaluate_population_fitness!(population,measurements,frequencies,parameter_bounds)
    init_descriptions = [feat_desc(circ) for circ in population]
    init_fitnesses = [circ.fitness for circ in population]
    min_fit, iter = 1,0;
    P = fill(0.0,(12,12));
    χ = Matrix{Circuit}(undef,(12,12));
    for (desc,fit,circ) in zip(init_descriptions,init_fitnesses,population) 
        if iszero(P[desc]) || P[desc] > fit
            P[desc] = fit
            χ[desc] = circ
        end
    end
    while (iter ≤ iterations) && (min_fit ≥ convergence_threshold)
        iter += 1
        for e in 1:population_size
            x′ = simplifycircuit(circuit_offspring_R(random_selection(χ,P),random_selection(χ,P),terminals))
            b′ = feat_desc(x′) #als input.
                @suppress evaluate_circuit_fitness!(x′,measurements,frequencies,parameter_bounds)
                p′ = x′.fitness
                if iszero(P[b′]) || P[b′] > p′
                    P[b′] = p′
                    χ[b′] = x′
                end
        end
        min_fit = minimum(P[P .> 0])
    end
    return P,χ,iter
end

function circuit_map_elites_continue_R(iterations,P,χ,measurements,frequencies,terminals="RCLP",feat_desc=fractional_descriptor,bounds=nothing) #continue evolving starting from a given archive
    parameter_bounds = Dict('R'=>[0,1.0e9],'C'=>[0,10],'L'=>[0,5],'P'=>[[0,0],[1.0e9,1]],'W'=>[0,1.0e9],'+'=>[0,0],'-'=>[0,0])
    if typeof(bounds) == Dict{Char, Vector}
        for key in keys(bounds)
            parameter_bounds[key] = bounds[key]
        end
    end
    iter = 0
    population_size = length(population_from_archive(P,χ))
    while (iter ≤ iterations) 
        iter += 1
        for e in 1:population_size
            x′ = simplifycircuit(circuit_offspring_R(random_selection(χ,P),random_selection(χ,P),terminals))
            b′ = feat_desc(x′) #als input.
                @suppress evaluate_circuit_fitness!(x′,measurements,frequencies,parameter_bounds)
                p′ = x′.fitness
                if iszero(P[b′]) || P[b′] > p′
                    P[b′] = p′
                    χ[b′] = x′
                end
        end
    end
    return P,χ
end

function initializepopulation_R(bounds,size=20,head=8,terminals="RCLP")
    return [initializecircuit_R(bounds,head,terminals) for i in 1:size]
end

function initializecircuit_R(bounds,head=8,terminals="RCLP")
    karva = generatekarva_R(head,terminals) 
    parameters = karva_parameters_(karva,bounds)
    return Circuit(karva,parameters,nothing)
end

function generatekarva_R(head,terminals="RCLP")
    operators =  "+-+-"
    all_elements = operators*terminals
    karva = "+R"*String(rand(all_elements,head-2))*String(rand(terminals,head+1))
    return karva
end

function one_point_crossover_R(circuit1,circuit2)
    karva1 = circuit1.karva
    karva2 = circuit2.karva
    crossover_point = rand(3:(length(karva1)-1)) 
    new_karva = karva1[1:crossover_point]*karva2[crossover_point+1:end]
    new_parameters =  vcat(circuit1.parameters[1:crossover_point],circuit2.parameters[crossover_point+1:end])
    return Circuit(new_karva,new_parameters,nothing)
end

function two_point_crossover_R(circuit1,circuit2)
    karva1 = circuit1.karva
    karva2 = circuit2.karva
    crossover_points = sort(rand(3:length(karva1)-1,2))
    new_karva = karva1[1:crossover_points[1]]*karva2[crossover_points[1]+1:crossover_points[2]]*karva1[crossover_points[2]+1:end] 
    new_parameters = vcat(circuit1.parameters[1:crossover_points[1]],circuit2.parameters[crossover_points[1]+1:crossover_points[2]],circuit1.parameters[crossover_points[2]+1:end])
    return Circuit(new_karva,new_parameters,nothing)
end


function circuit_offspring_R(circuit_1,circuit_2,terminals = "RCLP")
    offspring = ""
    rand_ = rand()
    if rand_ > 0.5
        offspring = one_point_crossover_R(circuit_1,circuit_2)
    elseif rand_ > 0.1
        offspring = two_point_crossover_R(circuit_1,circuit_2)
    else
        offspring = rand([circuit_1,circuit_2])
    end
    return mutate_R(offspring,terminals)
end

function mutate_R(circuit,terminals = "RCLP")
    karva = circuit.karva ; parameters = copy(circuit.parameters)
    total_length = length(karva)
    head_length = (total_length-1)/2 ; tail_length = head_length+1
    operations = "+-"; new_element = ""
    position = rand(3:total_length) #alternative: rand(3:coding_length(circuit))
    mutable_karva = collect(karva)
    if position == 1
        new_element = rand(operations)
        mutable_karva[position] = new_element
    elseif position <= head_length
        new_element = rand(operations*terminals)
        mutable_karva[position] = new_element
        parameters[position] = replacementparameter(new_element)
    elseif position > head_length
        new_element = rand(terminals)
        mutable_karva[position] = new_element
        parameters[position] = replacementparameter(new_element)
    end
    return Circuit(String(mutable_karva),parameters,nothing)
end

function Archive_nyquist(measurements,frequencies,P,χ,cutoff=1e-3)
# plot ground truth
    archive_nyquist_plot = scatter(real(measurements),-imag(measurements),label="Measured",legend = :outertopright)
# get population from archive 
    circuit_population = filter_fitness(population_from_archive(P,χ),cutoff)
    circuit_population_unique = make_unique(circuit_population)
    circuit_population_sorted = sort_by_complexity(circuit_population_unique)
# plot nyquist for each circuit in archive
    for circ in circuit_population_sorted
        simulations = simulateimpedance_noiseless(circ,frequencies)
        plot!(real(simulations),-imag(simulations),label=readablecircuit(circ))
    end
    xlabel!("Real{Z}/Ω")
    ylabel!("-Imag{Z}/Ω")
    return archive_nyquist_plot
end


function sort_by_complexity(circuit_population)
    return circuit_population[sortperm(coding_length.(circuit_population))]
end

function filter_fitness(circuit_population,cutoff)
    filter(x->x.fitness<cutoff,circuit_population)
end

### Version with multiple entries (n) at each point in the feature space.

function circuit_map_elites_R_m(measurements,frequencies;population_size= 30,iterations=10,feat_desc=fractional_descriptor,convergence_threshold = 5e-5,bounds = nothing,terminals="RCLP",head=8,initial_population=nothing, pe = 3)
    parameter_bounds = Dict('R'=>[0,1.0e9],'C'=>[0,10],'L'=>[0,5],'P'=>[[0,0],[1.0e9,1]],'W'=>[0,1.0e9],'+'=>[0,0],'-'=>[0,0])
    if typeof(bounds) == Dict{Char, Vector}
        for key in keys(bounds)
            parameter_bounds[key] = bounds[key]
        end
    end
    if isnothing(initial_population)  
        population = initializepopulation_R(parameter_bounds,population_size,head,terminals) 
    elseif typeof(initial_population) == Vector{Circuit}
        population = initial_population
    elseif typeof(initial_population) in [Vector{String},Vector{Tuple{String, Vector{Float64}}}]
        population = population_from_string(initial_population, head, population_size, terminals)
    end
    simplifypopulation!(population,terminals)
    @suppress evaluate_population_fitness!(population,measurements,frequencies,parameter_bounds)
    init_descriptions = [feat_desc(circ) for circ in population]
    init_fitnesses = [circ.fitness for circ in population]
    min_fit, iter = 1,0;
    P = fill(0.0,(12,12,pe)); #Change dimensions
    χ = Array{Circuit}(undef,(12,12,pe));
    for (desc,fit,circ) in zip(init_descriptions,init_fitnesses,population) 
        archive_update!(desc,fit,circ,P,χ)
    end
    while (iter ≤ iterations) && (min_fit ≥ convergence_threshold)
        iter += 1
        for e in 1:population_size
            x′ = simplifycircuit(circuit_offspring_R(random_selection(χ,P),random_selection(χ,P),terminals)) #Random Selection also needs to be updated
            b′ = feat_desc(x′) #als input.
                @suppress evaluate_circuit_fitness!(x′,measurements,frequencies,parameter_bounds)
                p′ = x′.fitness
                archive_update!(b′,p′,x′,P,χ)
        end
        min_fit = minimum(P[P .> 0])
    end
    return P,χ,iter
end

function archive_update!(desc,fit,circ,P,χ)
    # Check if archive contains at least one zero entry at feature descriptor location.
    firstzero_index = findfirst(x -> x == 0, P[desc[1],desc[2],:])
    if !isnothing(firstzero_index)
        P[desc[1],desc[2],firstzero_index] = fit
        χ[desc[1],desc[2],firstzero_index] = circ
    elseif maximum(P[desc[1],desc[2],:]) > fit
        P[desc[1],desc[2],findmax(P[desc[1],desc[2],:])[2]] = fit
        χ[desc[1],desc[2],findmax(P[desc[1],desc[2],:])[2]] = circ
    end
end


function circuit_map_elites_continue_R_m(iterations,P,χ,measurements,frequencies,terminals="RCLP",feat_desc=fractional_descriptor,bounds=nothing) #continue evolving starting from a given archive
    parameter_bounds = Dict('R'=>[0,1.0e9],'C'=>[0,10],'L'=>[0,5],'P'=>[[0,0],[1.0e9,1]],'W'=>[0,1.0e9],'+'=>[0,0],'-'=>[0,0])
    if typeof(bounds) == Dict{Char, Vector}
        for key in keys(bounds)
            parameter_bounds[key] = bounds[key]
        end
    end
    iter = 0
    population_size = length(population_from_archive(P,χ))
    while (iter ≤ iterations) 
        iter += 1
        for e in 1:population_size
            x′ = simplifycircuit(circuit_offspring_R(random_selection(χ,P),random_selection(χ,P),terminals))
            b′ = feat_desc(x′) #als input.
                @suppress evaluate_circuit_fitness!(x′,measurements,frequencies,parameter_bounds)
                p′ = x′.fitness
            archive_update!(b′,p′,x′,P,χ)
        end
    end
    return P,χ
end



function plot_MAP_archive(p)

    map_arch = plot_MAP_archive(p[:,:,1])
    for i in 2:size(p)[3]
        map_arch = plot_MAP_archive!(p[:,:,i])
    end
    return map_arch
end

function Archive_nyquist_m(measurements,frequencies,P,χ,cutoff=1e-3)
    # plot ground truth
        archive_nyquist_plot = scatter(real(measurements),-imag(measurements),label="Measured",legend = :outertopright)
    # get population from archive 
        circuit_population = filter_fitness(population_from_archive(P,χ),cutoff)
        circuit_population_sorted = sort_by_complexity(circuit_population)
    # plot nyquist for each circuit in archive
        for circ in circuit_population_sorted
            simulations = simulateimpedance_noiseless(circ,frequencies)
            plot!(real(simulations),-imag(simulations),label=readablecircuit(circ))
        end
        xlabel!("Real{Z}/Ω")
        ylabel!("-Imag{Z}/Ω")
        return archive_nyquist_plot
    end


    function population_from_archive_m(P,χ)
        return χ[findall(!iszero,P)]
     end
     

    function plot_MAP_archive!(p::Array{Float64, 2}, cam = (40, 30))
        X = []
        Y = []
        Z = []
        for x in 1:size(p)[1]
            for y in 1:size(p)[2]
                push!(X,x)
                push!(Y,y)
                push!(Z,p[x,y])
            end
        end
        Z_ = Z[findall(!iszero,Z)]
        X_ = X[findall(!iszero,Z)]
        Y_ = Y[findall(!iszero,Z)]
        fig = scatter!(X_,Y_,log10.(Z_),legend=false,camera = cam,color=:green)  
        xlabel!("Integer-order count")
        ylabel!("Fractional-order count")
        zlabel!("log10 fitness")
        title!("Feature-performance plot")
        return fig
    end
    
    function make_unique(circuitpopulation)
        newpop = []
        newpop_strings = []
        for circ in circuitpopulation
            if !(readablecircuit(circ) in newpop_strings)
                push!(newpop_strings,readablecircuit(circ))
                push!(newpop,circ)
            end
        end
        return newpop
    end
