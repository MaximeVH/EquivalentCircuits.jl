function savepopulation(fname,population)
    fname = fname[end-3:end] == ".csv" ? fname : fname*".csv"
    array = Array{Any,2}(undef,length(population),2+length(population[1].karva))
        for (n,circuit) in enumerate(population)
            array[n,:] = vcat([circuit.karva,circuit.fitness == nothing ? Inf : circuit.fitness],circuit.parameters)
        end
        writedlm(fname,array,';')
end

"""
    loadpopulation(filepath::String)

Load a library of equivalent electrical circuits, to be provided as initial population for the circuitevolution function.

"""
function loadpopulation(filepath)
    population_array = readdlm(filepath,';')
    population_size = size(population_array,1)
    population = Array{Circuit}(undef,population_size)
    cpe_indices = findall(x->typeof(x)==SubString{String},population_array[:,:])
    population_array[cpe_indices[size(population_array,1)+1:end]] = eval.(Meta.parse.(population_array[cpe_indices[population_size+1:end]]))
    for individual in 1:population_size
        population[individual] = Circuit(population_array[individual,1],vcat(population_array[individual,3:end]),population_array[individual,2] == Inf ? nothing : population_array[individual,2])
    end
    return population
end

function readmeasurements(file)
    meansurement_file = readdlm(file,',')
    reals = meansurement_file[:,1]
    imags = meansurement_file[:,2]
    frequencies = meansurement_file[:,3]
    measurements = reals + imags*im
    return measurements,frequencies
end

"""
initialize_circuitlibrary(circuit,parameters)

Initializes a circuit library in the current directory. the circuit
argument should be provided in circuit notation form (e.g. "R1-[C2,R3-[C4,R5]]"), 
and the parameters should be provided as an array (e.g. [20,4e-9,3400,4e-6,2500]).
"""

function initialize_circuitlibrary(circuit,parameters=nothing,type = "N",source = "no information provided.")
    @assert !isfile("Circuitlibrary.csv") "A circuitlibrary already exists in the current directory."
    if isnothing(parameters)
        parameters = circuitparameters(circuit)
    end
    circuit_code,circuit_parameters = circuit_to_karva(circuit,parameters)
    circuitlibrary = [Circuit(circuit_code,circuit_parameters,nothing)]
    circuitmetadata = Array{Any}(undef,1,4)
    circuitmetadata[1:4] = [circuit,parameters,type,source]
    savepopulation("Circuitlibrary.csv",circuitlibrary)
    writedlm("circuitmetadata.csv",circuitmetadata,';')
end

"""
add_to_circuitlibrary(circuit,parameters)

Adds a new circuit to an already initialized circuit library. the circuit
argument should be provided in circuit notation form (e.g. "R1-[C2,R3-[C4,R5]]"), 
and the parameters should be provided as an array (e.g. [20,4e-9,3400,4e-6,2500]).
"""
function add_to_circuitlibrary(circuit,parameters=nothing,type = "N",source = "no information provided.") #path to circuit library should be able to be provided in function arguments.
    @assert isfile("CircuitLibrary.csv") "A circuitlibrary has not yet been initialized in the current directory."
    if isnothing(parameters)
        parameters = circuitparameters(circuit)
    end
    #convert the circuit to its basic encoding and parameter values.
    circuit_code,circuit_parameters = circuit_to_karva(circuit,parameters)
    circuitobject = Circuit(circuit_code,circuit_parameters,nothing)
    #elongate with random terminals until required length, dictated by the head length.
    #load the circuitlibrary and circuitlibrary metadata files.
    circuitlibrary = loadpopulation("CircuitLibrary.csv")
    metadatalibrary = readdlm("circuitmetadata.csv",';')
    #append the circuit.
    push!(circuitlibrary,circuitobject)
    circuitmetadata = Array{Any}(undef,1,4)
    circuitmetadata[1:4] = [circuit,parameters,type,source]
    metadatalibrary = vcat(metadatalibrary,circuitmetadata)
    #save the circuitlibrary files.
    savepopulation("Circuitlibrary.csv",circuitlibrary)
    writedlm("Circuitmetadata.csv",metadatalibrary,';')
end


function add_to_circuitlibrary_limitedparams(circuit,limitedparams,type = "N",source = "no information provided.") #when only a part of the parameters are known.
    @assert isfile("CircuitLibrary.csv") "A circuitlibrary has not yet been initialized in the current directory."
    parameters = circuitparameters(circuit,limitedparams)
    #convert the circuit to its basic encoding and parameter values.
    circuit_code,circuit_parameters = circuit_to_karva(circuit,parameters)
    circuitobject = Circuit(circuit_code,circuit_parameters,nothing)
    #elongate with random terminals until required length, dictated by the head length.
    #load the circuitlibrary and circuitlibrary metadata files.
    circuitlibrary = loadpopulation("CircuitLibrary.csv")
    metadatalibrary = readdlm("circuitmetadata.csv",';')
    #append the circuit.
    push!(circuitlibrary,circuitobject)
    circuitmetadata = Array{Any}(undef,1,4)
    circuitmetadata[1:4] = [circuit,parameters,type,source]
    metadatalibrary = vcat(metadatalibrary,circuitmetadata)
    #save the circuitlibrary files.
    savepopulation("Circuitlibrary.csv",circuitlibrary)
    writedlm("Circuitmetadata.csv",metadatalibrary,';')
end

function add_encoding_to_circuitlibrary(circuit_code,type = "N",source = "no information provided.") 
    circuit_parameters = karva_parameters(circuit_code)
    additional_codings = join(rand("RCLP",17-length(circuit_code))) #assuming a headlength of 8.
    addidional_parameters = karva_parameters(additional_codings)
    karva = circuit_code*additional_codings
    parameters = vcat(circuit_parameters,addidional_parameters)
    circuitlibrary = loadpopulation("CircuitLibrary.csv")
    metadatalibrary = readdlm("circuitmetadata.csv",';')
    circuitobject = Circuit(karva,parameters,nothing)
    push!(circuitlibrary,circuitobject)
    circuitmetadata = Array{Any}(undef,1,4)
    tree =  karva_to_tree(circuit_code,circuit_parameters)
    circuit,parameters_ = tree_to_circuit(tree)
    circuitmetadata[1:4] = [circuit,parameters_,type,source]
    metadatalibrary = vcat(metadatalibrary,circuitmetadata)
     #save the circuitlibrary files.
     savepopulation("Circuitlibrary.csv",circuitlibrary)
     writedlm("Circuitmetadata.csv",metadatalibrary,';')
end

function add_encoding_to_circuitlibrary_with_parameters(circuit_code,parameters,type = "N",source = "no information provided.") 
    additional_codings = join(rand("RCLP",17-length(circuit_code))) #assuming a headlength of 8.
    addidional_parameters = karva_parameters(additional_codings)
    karva = circuit_code*additional_codings
    parameters = vcat(parameters,addidional_parameters)
    circuitlibrary = loadpopulation("CircuitLibrary.csv")
    metadatalibrary = readdlm("circuitmetadata.csv",';')
    circuitobject = Circuit(karva,parameters,nothing)
    push!(circuitlibrary,circuitobject)
    circuitmetadata = Array{Any}(undef,1,4)
    tree =  karva_to_tree(circuit_code,circuit_parameters)
    circuit,parameters_ = tree_to_circuit(tree)
    circuitmetadata[1:4] = [circuit,parameters_,type,source]
    metadatalibrary = vcat(metadatalibrary,circuitmetadata)
     #save the circuitlibrary files.
     savepopulation("Circuitlibrary.csv",circuitlibrary)
     writedlm("Circuitmetadata.csv",metadatalibrary,';')
end


function circuitparameters(circuit)
    elements = foldl(replace,["["=>"","]"=>"","-"=>"",","=>""],init = denumber_circuit(circuit))
    return karva_parameters2(elements)
end

function circuitparameters(circuit,limited_params)
    elements = foldl(replace,["["=>"","]"=>"","-"=>"",","=>""],init = denumber_circuit(circuit))
    return karva_parameters3(elements,limited_params)
end

function karva_parameters2(karva) 
    parameters = Array{Float64}(undef, length(karva))
    ranges = Dict('R'=>4000,'C'=> 0.0001,'L'=> 1,'+'=> 0,'-'=> 0,'P'=>[4000,1])
    for (e,i) in enumerate(karva)
            if e == 'P'
                Element_values[n] = [ranges[e][1]*rand(),ranges[e][2]*rand()]
            else
            parameters[e] = rand()*ranges[i]
        end
    end
    return parameters
end

function karva_parameters3(karva,limited_params) 
    parameters = Array{Any}(undef, length(karva))
    ranges = Dict('R'=>4000,'C'=> 0.0001,'L'=> 1,'+'=> 0,'-'=> 0,'P'=>[4000,1])
    for (e,(i,l)) in enumerate(zip(karva,limited_params))
            if l != 0
                parameters[e] = l
            else
                if i == 'P'
                    Element_values[e] = [ranges[i][1]*rand(),ranges[i][2]*rand()]
                else
                parameters[e] = rand()*ranges[i]
            end
        end
    end
    return parameters
end

function nyquist_parameter_optimisation(circuit,measurements,frequencies,first_method=:de_rand_1_bin)
    elements = foldl(replace,["["=>"","]"=>"","-"=>"",","=>""],init = denumber_circuit(circuit))
    initial_parameters = flatten(karva_parameters(elements));
    circfunc = circuitfunction(circuit)
    objective = objectivefunction(circfunc,measurements,frequencies)
    lower = zeros(length(initial_parameters))
    upper = get_parameter_upper_bound(circuit)

    ### First step ###
    SR = Array{Tuple{Float64,Float64},1}(undef,length(initial_parameters))
    for (e,(l,u)) in enumerate(zip(lower,upper))
        SR[e] = (l,u)
    end
    res = bboptimize(objective; SearchRange = SR, Method = first_method,MaxSteps=170000,TraceMode = :silent);
    initial_parameters = best_candidate(res)
    fitness_1 = best_fitness(res)
    ### Second step ###
    inner_optimizer = NelderMead()
    results = optimize(objective, lower, upper, initial_parameters, Fminbox(inner_optimizer), Optim.Options(time_limit = 50.0)); #20.0
    fitness_2 = results.minimum
    best = results.minimizer

    parameters = fitness_2 < fitness_1 ? best : initial_parameters
    fitness = max(fitness_2,fitness_1)  
    return parameters,fitness
end

function element_count(circuit)
    return sum([occursin(t,"RCLPW") for t in circuit])
end

function circuit_literaturesearch(measurements,frequencies,terminals,domain,fitness_threshold,max_complexity) #both the circuitlibrary and corresponding metadata are required here.
    #Filtering of the circuit library, using the input terminal and domain info.
    circuit = ["P1-R2-[C3,R4]","R1-[C2,R3]","P1-[R2,C3]-P4","R1-[C2,R3]-R4-P5","[R1,C2-R3-[C4,R5]]","[R1,[C2,R3]-R4]","[R1,C2-[R3,C4-R5]]","R1-[P2,R3]",          "[R1-C2,R3]",        "W1-[R2,C3]-R4",       "R1-[P2,R3]-[P4,R5]",  "[R1,[R2,R3-C4]]",  "[R1,C2]",           "[R1,C2-[R3-C4,R5]]","[R1-[P2,R3],C4]","[C1-R2-C3,C4]","R1-[P2,R3]","R1-[R2,C3]","R1-[R3-W4,C2]","C1-R2-C3","R1-[C2,[P3,R4]-R5]","R1-[R2-W3,C4]","R1-[W2-P3,R4]","R1-[W2-R3,P4]","R1-[P2,R3-W4]","R1-[C2-C3-C4,R5-W6]","R1-P2-[R3,C4]","R1-[P2,R3-P4]","R1-[P2,R3-W4]-[P5,R6]","R1-[C2,R3-W4]-[C5,R6]","R1-[C2,R3]-P4","R1-C2-[R3,P4]","R1-[R2,C3]-[R4,C5]-W6","[C1-R2,[C3-R4,R5]-[C6-R7,R8]]-R9-L10","[R1-[R2,C3],C4]","R1-[[R2,C3],L4,R5]","R1-[C2,R3]-[C4,W5]","R1-[R2-P3,P4]","R1-[W2,[R3,C4-W5]]","R1-[[C2,R3],R4-L5]","P1-R2-L3-[R4,P5]-[R6,P7]","[C1,R2-C3]","R1-[P2,R3-[P4,R5]]","R1-[C2,R3]-[C4,R5]","R1-[P2,R3]-[L4,R5]-P6", "[C1,[W2,R3]-W4]","R1-[C2,R3-W4]"]
    domains = ["Animals","Animals","Animals","Animals","Animals","Plants","Plants","Plants","Plants","Plants","Plants","Plants","Plants","Plants","Biosensors","Biosensors","Biosensors","Biosensors","Biosensors","Biosensors","Biosensors","Biosensors","Biosensors","Biosensors","Biosensors","Biosensors","Biosensors","Biosensors","Biosensors","Biosensors","Biosensors","Biosensors","Batteries","Batteries","Batteries","Batteries","Batteries","Batteries","Batteries","Batteries","Batteries","Materials","Materials","Materials","Materials","Materials","Materials"]
    @assert domain ∈ unique(domains) "Invalid domain provided."
    domain_inds = findall(x-> x == domain,domains)
    circuit = circuit[domain_inds]
    complexity_inds = findall(x-> x <= max_complexity,element_count.(circuit))
    circuit = circuit[complexity_inds]
    terminal_complement = join(["R","W","P","L","C"][[!occursin(x,terminals) for x in ["R","W","P","L","C"]]])
    terminal_inds = findall(x->!occursin(terminal_complement,x),circuit)
    circuit = circuit[terminal_inds]
    parameters = []
    fitnesses = []
    for c in circuit
        params,fitness =  nyquist_parameter_optimisation(c,measurements,frequencies)
        push!(parameters,params)
        push!(fitnesses,fitness)
    end
    fitness_inds = findall(x -> x <= fitness_threshold, fitnesses)
    fitnesses = fitnesses[fitness_inds]
    circuit = circuit[fitness_inds]
    if isempty(circuit)
        println("No available circuits adequately fit the given measurement data.")
    else
        sp = sortperm(fitnesses)
        return circuit[sp]
    end
end