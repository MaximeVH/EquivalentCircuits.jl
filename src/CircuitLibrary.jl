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

function get_literature_info()
    animals = ["P1-R2-[C3,R4]","R1-[C2,R3]","P1-[R2,C3]-P4","R1-[C2,R3]-R4-P5","[R1,C2-R3-[C4,R5]]"]
    plants = ["[R1,[C2,R3]-R4]",       "[R1,C2-[R3,C4-R5]]",  "R1-[P2,R3]",            "[R1-C2,R3]",          "W1-[R2,C3]-R4" ,        "R1-[P2,R3]-[P4,R5]",    "[R1,[R2,R3-C4]]",     "[R1,C2]",            "[R1,C2-[R3-C4,R5]]",  "R1-[R2,C3]-[R4-P5,[C6,[P7,[R8,R9-C10]]"]
    biosensors = ["[C1-R2-C3,C4]","R1-[P2,R3]","R1-[R2,C3]","R1-[R2-W3,C4]","R1-[P2,P3]","C1-R2-C3","R1-[C2,[P3,R4]-R5]","R1-[R2-W3,C4]","R1-[W2-P3,R4]","R1-[P2,R3-W4]","R1-[C2-C3-C4,R5-W6]", "R1-P2-[R3,C4]","R1-[P2,R3-P4]","R1-[P2,R3-W4]-[P5,R6]","R1-[C2,R3-W4]-[C5,R6]","R1-[C2,R3]-P4","R1-C2-[R3,P4]"]
    batteries = ["R1-[R2,C3]-[R4,C5]-W6","[R1-[R2,C3],C4]","R1-[[R2,C3],L4,R5]","R1-[C2,R3]-[C4,W5]","R1-[R2-P3,P4]","R1-[W2,[R3,C4-W5]]","R1-[[C2,R3],R4-L5]","L1-R2-[P3,R4]-[P5,R6]-[P7,R8]-W9","L1-R2-[P3,R4]-[P5,R6]-[P7,R8]","L1-R2-[P3,R4]-[P5,R6]","L1-R2-[P3,R4]","R1-[P2,R3-P4]","P1-R2-L3-[R4,P5]-[R6,P7]","R1-[C2,R3]-[C4,R5-W6]","R1-[C2,R3]-[P4,R5-W6]","L1-R2-[P3,R4-W6]","R1-[R2,P3]-[R4,P5]-[R6,P7]-P8","R1-[R2,C3]-[R4,C5]-W6-C7","R1-[R2,C3]-[R4,C5]","R1-[R2,P3]-[R4-W5,P6]-[R7,P8]-P9","R1-[L2,R3]-[C4,R5]-[C5,R6]","L1-R2-[C3,R4]-P5-C6","L1-R2-[P3,R4]-[P5,R6-W7]"]
    fuel_cells = ["R1-[P2,R3-W4]","[C1-R2,[C3-R4,R5]-[C6-R7,R8]]-R9-L10" ,"[P1,R2]-R3-[P4,R5]" ,"[C1,R2]-R3-[C4,R5]" ,"[P1,R2]-R3-[P4,R5]-L6" ,"R1-[P2,R3]-[[R4-L5,P6],R7]" ,"R1-L2-[P3,[R4,R5-L6]]" ,"R1-[C2,R3-[C4-R5]]-[C6,R7]","R1-[R2,C3]-[R4,C5]-[R6,C7]","R1-[R2,L3]-[R4,P5]-[R6,P7]-[R6,C9]","R1-[R2,C3]-[R4,C5]-[R6,C7]-L8" ,"R1-[R2,C3]-[R4,C5]-L6 ","R1-[P2,R3]-[P4,R5-[C6,R7]]-[C8,[R9,R10-L11]]" ,"R1-[P2,R3]-[P4,R5-W6]","R1-[P2,R3]-[P4,R5]-[W6,L7]","R1-[C2,R3-[C4,R5-P6]]","R1-[P2,R3]-[P4,R5]-[P6,R7]-L8","L1-R2-[C3,R4-[p5,R6]]" ,"[C1,[R2,R3-L4]]","[P1,R2]-R3-[P4,R5-W6]","R1-[C2,[R3-W4,R5-C6]]","L1-R2-[C3,R4-[C5,R6-[L7,R8]]" ,"R1-[C2,R3-[C4,R5-W6]]"]
    supercapacitors = ["L1-R2-[C3,P4-R5]","R1-[C2,R3-W4]","R1-P2-W3","R1-[P2,R3]-W4","R1-[P2,R3]-[P4,R5]","L1-R2-[P3,R4]","L1-R2-[P3,R4]-W5","R1-[C2,R3-W4]-C5","R1-[R2,P3]-P4","R1-[P2,R3-W4]-C5","R1-[P2,R3-[P4,R5-[P6,R7-P8]]]","R1-[R2,C3-W4]-[R5,P6]","R1-[P2,R3-[P4,R5-W6]]","R1-[P2,R3]-[P4,R5-W6]","R1-[R2,P3]-W4","R1-[R2,P3]-[R4,P5]-W6","R1-[R2,P3]-[R4-W5,P6]-P7","L1-R2-[P3,R4]-[P5,R6]-[C7,R8]","R1-[C2,R3-P4-W5]","R1-[C2,R3-P4]","R1-[P2,R3-W4]-[P5,R6]","R1-[P2,R3-W4]-P5","R1-[P2,R3-W4]","R1-[C2,P3-W4]-[C5,R6]","R1-C2","R1-[R2,P3]-[R4,C5]","R1-[R2-P3,P4] ","R1-[C2,W3]-[C4,R5]","L1-R2-[R3,P4]","R1-[R2,P3]-[R4,P5]-[R5,P6]","R1-[R2,C3]","R1-P2","R1-[R2,P3]-W4","[P1,W2-R3]-[C4,R5]-[R6,P7]","R1-[P2,R3-W4]-[R5,C6]","R1-[C2,R3]-[P4,R5]-[C6,R7]-[C8,R9]","R1-[C2,R3]-[C4,R5]-[P6,R7-W8]"]
    material_science = ["[C1,R2-C3]","R1-[P2,R3-[P4,R5]]","R1-[P2,R3]-[P4,R5]","R1-[P2,R3]-[L4,R5]-P6","R1-[[C2,R3-[C4,R5]],L6-R7]","[C1,[W2,R3]-W4]","R1-[C2,R3-W4]"]
    circuits = vcat(animals,plants,biosensors,batteries,fuel_cells,supercapacitors,material_science)
    categories = vcat(split(repeat("Animals ",length(animals))," ")[1:end-1],
    split(repeat("Plants ",length(plants))," ")[1:end-1],
    split(repeat("Biosensors ",length(biosensors))," ")[1:end-1],
    split(repeat("Batteries ",length(batteries))," ")[1:end-1],
    split(repeat("Fuel_cells ",length(fuel_cells))," ")[1:end-1],
    split(repeat("Supercapacitors ",length(supercapacitors))," ")[1:end-1],
    split(repeat("Materials ",length(material_science))," ")[1:end-1])
    animal_DOIs = ["10.1016/j.elecom.2020.106742","10.1007/s10439-010-0127-y","10.1007/s10439-016-1559-9","10.1016/j.bios.2015.10.027","10.1016/S0302-4598(98)00079-8"]
    plant_DOIs = ["10.1093/jxb/20.2.177","10.1093/jxb/42.11.1465","10.1007/s11694-017-9545-y","10.1016/j.lwt.2016.03.019","10.1016/j.postharvbio.2019.110978", "10.1093/jexbot/51.353.2095","10.1016/j.protcy.2016.03.024","10.1111/j.1365-2621.2011.02636.x","10.1016/S0925-5214(99)00056-3","10.1007/s13197-019-03590-3"]
    biosensor_DOIs = ["10.1016/j.talanta.2012.02.056","10.1016/j.bios.2014.09.037","10.1016/j.bios.2020.112709","10.1021/nl080366q","10.1016/j.bios.2012.03.001","10.1016/j.snb.2014.01.098","10.1016/j.bios.2012.03.001","10.1016/j.bios.2018.08.064","10.1016/j.snb.2012.06.004","10.1016/j.snb.2016.07.121","10.1016/j.aca.2010.12.018","10.1016/j.snb.2015.11.059","10.1016/j.bios.2016.01.057","10.1016/j.elecom.2009.12.017","10.1016/j.bios.2012.10.054","10.1021/ac9001854","10.1007/s11694-010-9101-5"]
    batteries_DOIs = ["10.1016/j.egyr.2020.03.029","10.1016/0013-4686(95)00026-B","10.1149/1.2221040","10.1016/0378-7753(90)85020-D","10.3390/en12234507","10.1016/0378-7753(93)80086-5","10.1149/1.1467940","10.1007/s40243-015-0052-y","10.1016/j.electacta.2020.135787","10.1016/j.electacta.2020.135787","10.1016/j.electacta.2020.135787","10.1007/s40243-015-0052-y","10.1016/j.jpowsour.2010.12.102","10.1109/ITEC-India48457.2019.ITECINDIA2019-232","10.3906/kim-1910-72","10.1109/TPEL.2017.2780184","10.1016/j.electacta.2020.136944","10.1016/j.est.2016.09.001","10.1109/TPEL.2020.2977274","10.1016/j.ssi.2021.115680","10.1016/j.est.2017.08.004","10.1016/j.jpowsour.2021.229454","10.1002/er.5350"] #There's one missing here.
    fuel_cell_DOIs = ["10.1016/j.ijhydene.2017.07.094","10.1109/APEC.2004.1295834","10.1016/j.jpowsour.2014.10.138","10.1016/j.jpowsour.2018.08.082","10.1016/j.apenergy.2020.115718","10.1007/s40243-015-0052-y","10.1016/j.energy.2020.118185","10.1115/1.2931462","10.1109/TIE.2015.2393292","10.1016/j.ijhydene.2011.04.076","10.1016/j.jpowsour.2014.01.049","10.1016/j.jpowsour.2010.02.054","10.1016/j.electacta.2016.01.128","10.3390/fractalfract5010021","10.3390/fractalfract5010021","10.1016/j.apenergy.2019.113396","10.1016/j.electacta.2008.04.047","10.1016/j.electacta.2011.03.056","10.1016/S0378-7753(99)00331-6","10.1016/j.ijhydene.2013.01.132","10.1016/j.electacta.2019.04.102","10.1016/j.ijhydene.2007.03.040","10.1016/j.electacta.2020.136036"]
    supercapacitor_DOIs = ["10.1016/j.electacta.2020.137412","10.1002/er.5634","10.1109/IFEEC47410.2019.9015091","10.1109/IFEEC47410.2019.9015091","10.1016/j.electacta.2019.134706","10.1007/s10965-020-02183-5","10.1016/j.jpowsour.2018.07.108","10.1016/j.jpcs.2018.04.044","10.1002/elan.202060239","10.1016/j.est.2021.102738","10.1016/j.est.2021.103494","10.1007/s10854-019-02686-y","10.1016/j.electacta.2018.02.116","10.1016/j.jmst.2016.05.008","10.3390/en14133807","10.3390/en14133807","10.3390/en14133807","10.3390/en14133807","10.1016/S0013-4686(02)00091-9","10.1016/j.jelechem.2019.01.022","10.1016/j.matpr.2020.07.248","10.1016/j.matpr.2020.07.248","10.1016/j.electacta.2018.03.091","10.1016/j.electacta.2017.12.167","10.1016/j.surfin.2020.100524","10.1109/TR.2018.2869212","10.1007/s12034-020-02139-x","10.1016/j.est.2021.102328","10.3390/en14041139","10.1016/j.electacta.2021.137746","10.3390/molecules24081452","10.1134/S1811238218020194","10.3390/nano11051062","10.3390/ma11010048","10.1016/j.electacta.2018.03.091","10.1016/j.polymer.2020.122954","10.1109/TPEL.2018.2810889"]
    material_science_DOIs = ["10.1016/S0008-8846(02)00720-2","10.1016/j.corsci.2008.02.010","10.3390/ma7010218","10.1016/S0379-6779(01)00667-1","10.1016/j.electacta.2013.12.124","10.3390/coatings9040254","10.1016/j.msec.2007.10.081"]
    DOIs = vcat(animal_DOIs,plant_DOIs,biosensor_DOIs,batteries_DOIs,fuel_cell_DOIs,supercapacitor_DOIs,material_science_DOIs)
 return circuits, categories, DOIs
end

"""
circuit_search(measurements::Array{Complex{Float64},1},frequencies::Array{Float64,1},domain::String; <keyword arguments>)

Evaluates the compatibility of a given set of measurement data with circuits from the literature, for a specific application domain.

The inputs are a circuit (e.g. "R1-[C2,R3]-P4"), an array of complex-valued impedance measurements, their corresponding
frequencies, and an application domain (e.g. "Batteries"). The output is an array of fitting circuit configurations along with the Digital Object Identifiers (DOI)
of the accompanying literature.

# Arguments
- `terminals::String="RCLPW"`: the circuit components that are to be included in the circuit literature search.
- `max_complexity::Integer=20`: a hyperparameter than controls the maximum considered complexity of the circuits.
- `fitness_threshold::10e-5`: The objective function threshold under which a circuit is considered to fit the measurements.
"""
function circuit_search(measurements,frequencies,domain;terminals="RCLPW",fitness_threshold=10e-5,max_complexity=20)
    circuit,domains,dois = get_literature_info()
    @assert domain ∈ unique(domains) "Invalid domain provided."
    domain_inds = findall(x-> x == domain,domains)
    circuit = circuit[domain_inds]
    dois = dois[domain_inds]
    complexity_inds = findall(x-> x <= max_complexity,element_count.(circuit))
    circuit = circuit[complexity_inds]
    dois = dois[complexity_inds]
    terminal_complement = join(["R","W","P","L","C"][[!occursin(x,terminals) for x in ["R","W","P","L","C"]]])*"A"
    terminal_inds = findall(x->!occursin(terminal_complement,x),circuit)
    circuit = circuit[terminal_inds]
    dois = dois[terminal_inds]
    parameters = []
    fitnesses = []
    for c in circuit
        params,fitness =  nyquist_parameter_optimisation(c,measurements,frequencies);
        push!(parameters,params)
        push!(fitnesses,fitness)
    end
    fitness_inds = findall(x -> x <= fitness_threshold, fitnesses)
    fitnesses = fitnesses[fitness_inds]
    circuit = circuit[fitness_inds]
    dois = dois[fitness_inds]
    if isempty(circuit)
        println("No available circuits adequately fit the given measurement data.")
    else
        sp = sortperm(fitnesses)
        return hcat(circuit[sp],dois[sp])
    end
end