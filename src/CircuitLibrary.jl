using DelimitedFiles

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