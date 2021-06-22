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