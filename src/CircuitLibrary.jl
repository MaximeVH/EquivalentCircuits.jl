function savepopulation(fname,population)
    fname = fname[end-3:end] == ".csv" ? fname : fname*".csv"
    array = Array{Any,2}(undef,length(population),2+length(population[1].karva))
        for (n,circuit) in enumerate(population)
            array[n,:] = vcat([circuit.karva,circuit.fitness == nothing ? Inf : circuit.fitness],circuit.parameters)
        end
        writedlm(fname,array,';')
end

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

function initialize_circuitlibrary(circuit,parameters,type = "N",source = "no information provided.")
    @assert !isfile("Circuitlibrary.csv") "A circuitlibrary already exists in the current directory."
    circuit_code,circuit_parameters = circuit_to_karva(circuit,parameters)
    circuitlibrary = [Circuit(circuit_code,circuit_parameters,nothing)]
    circuitmetadata = Array{Any}(undef,1,4)
    circuitmetadata[1:4] = [circuit,parameters,type,source]
    savepopulation("Circuitlibrary.csv",circuitlibrary)
    writedlm("circuitmetadata.csv",circuitmetadata,';')
end

function initialize_circuitlibrary(circuit,type = "N",source = "no information provided.")
    @assert !isfile("Circuitlibrary.csv") "A circuitlibrary already exists in the current directory."
    parameters = circuitparameters(circuit)
    circuit_code,circuit_parameters = circuit_to_karva(circuit,parameters)
    circuitlibrary = [Circuit(circuit_code,circuit_parameters,nothing)]
    circuitmetadata = Array{Any}(undef,1,4)
    circuitmetadata[1:4] = [circuit,parameters,type,source]
    savepopulation("Circuitlibrary.csv",circuitlibrary)
    writedlm("circuitmetadata.csv",circuitmetadata,';')
end

function add_to_circuitlibrary(circuit,parameters,type = "N",source = "no information provided.") #path to circuit library should be able to be provided in function arguments.
    @assert isfile("CircuitLibrary.csv") "A circuitlibrary has not yet been initialized in the current directory."
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

function add_to_circuitlibrary(circuit,type = "N",source = "no information provided.")
    @assert isfile("CircuitLibrary.csv") "A circuitlibrary has not yet been initialized in the current directory."
    parameters = circuitparameters(circuit)
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

function add_encoding_to_circuitlibrary(circuit_code,circuit_parameters,type = "N",source = "no information provided.") 
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