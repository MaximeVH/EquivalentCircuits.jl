function circuitfunction(Circuit)
    b_new = true
    b_debug = true
    if b_debug println("Circuit, in: ", Circuit) end
    for (f, t) in zip(["-", "[",  ",",      "]"],
                      ["+", "((", ")^-1+(", ")^-1)^-1"])
        Circuit= replace(Circuit, f=>t)
    end
    if b_debug println("Circuit, after replace: ", Circuit) end
    for I in 2:-1:1
        Es = eachmatch(Regex("([CLRPW])([0-9]){$(I)}"),Circuit)
            for e in Es
                match = e.match
                if match[1] == 'C'
                    Circuit = replace(Circuit,match=>"(1/(2im*π*f*"*"T"*"))")
                    if b_debug println("circuitfunction, Circuit 'C': ", Circuit) end
                elseif match[1] == 'R'
                    Circuit = replace(Circuit,match=>"T")
                elseif match[1] == 'L'
                    Circuit = replace(Circuit,match=>"(2im*π*f*"*"T"*")")
                elseif match[1] == 'P'
                    if b_debug println("P: Circuit; in: ", Circuit, ", match: ", match) end
                    Circuit_P_new  = replace(Circuit,match=>"(1/(T*(2im*π*f)^N))")
                    Circuit_P_orig = replace(Circuit,match=>"T*(2*π*f)^(-N)"*"*(cos((π*N)*0.5)-sin((π*T)*0.5)im)")
                    if b_new
                        Circuit = Circuit_P_new
                        if b_debug println("new: circuit ") end
                    else
                    # original:
                        Circuit = Circuit_P_orig
                        if b_debug println("orig: circuit") end
                    end
                    if b_debug println("circuitfunction, Circuit 'P': ", Circuit_P_orig) end
                    if b_debug println("circuitfunction, new     'P': ", Circuit_P_new) end
                elseif match[1] == 'W'
                    Circuit = replace(Circuit,match=>"T*(2*π*f)^(-0.5)"*"*(cos((π*0.5)*0.5)-sin((π*0.5)*0.5)im)")
            end
        end
    end
    new_circuit = ""
    counter = 1
    for i in Circuit
        if i == 'T'
            new_circuit = new_circuit*"T["*string(counter)*"]"
            counter += 1
        elseif i == 'N'
            new_circuit = new_circuit*"T["*string(counter)*"]"
        else
            new_circuit = new_circuit*i
        end
    end

    Circuit_expression = Meta.parse(new_circuit)
    return mk_function([:T,:f],[], Circuit_expression)
end

function tree_to_function(tree_array)
    b_new = true
    dictionary_of_calculations_orig = Dict('R'=>"T",'C'=>"(1/(2im*π*f*"*"T"*"))",'L'=>"(2im*π*f*"*"T"*")",'P'=>"T[1]*(2*π*f)^(-T[2])"*"*(cos((π*T[2])*0.5)-sin((π*T[2])*0.5)im)") 
    dictionary_of_calculations_new = Dict('R'=>"T",'C'=>"(1/(2im*π*f*"*"T"*"))",'L'=>"(2im*π*f*"*"T"*")",'P'=>"(1/(T[1]*(2im*π*f)^(T[2])))"*"*(cos((π*T[2])*0.5)-sin((π*T[2])*0.5)im)") 

    if b_new
        dictionary_of_calculations = dictionary_of_calculations_new
        println("orig: ", dictionary_of_calculations_orig)
        println("new:  ", dictionary_of_calculations_new)
    else
        dictionary_of_calculations = dictionary_of_calculations_orig
    end
    tree = copy(tree_array) 
    nodecount = length(tree)
    essential_info = [[node.ParentIndex,node.Type,node.Index] for node in tree]
    for i in 1:((length(essential_info)-1)/2)
        parent1,type1,index1 = pop!(essential_info)
        parent2,type2,index2 = pop!(essential_info)
        type1 = haskey(dictionary_of_calculations,type1) ? replace(dictionary_of_calculations[type1],"T"=>"T["*string(index1)*']') : type1
        type2 = haskey(dictionary_of_calculations,type2) ? replace(dictionary_of_calculations[type2],"T"=>"T["*string(index2)*']') : type2
        operation = essential_info[parent1][2]
        if operation == '+'
            essential_info[parent1][2] = type1*'+'*type2
        elseif operation == '-'
            essential_info[parent1][2] = "((("*type1*")^-1)+("*type2*")^-1)^-1"
        end
    end
    CircuitExpression = Meta.parse(essential_info[1][2])
return mk_function([:T,:f],[],CircuitExpression)
end

karva_to_function(karva::String) = tree_to_function(karva_to_tree(karva))
karva_to_function(circuit::Circuit) = karva_to_function(circuit.karva)

function flatten(params)
    new_array = Float64[]
    for param in params
        if length(param) == 1
            push!(new_array,param)
        else
            push!(new_array,param[1])
            push!(new_array,param[2])
        end
    end
    return new_array
end

function get_target_impedance(circuit,circuit_parameters,frequency=1000)
    circfunc = circuitfunction(circuit)
    target_impedance = circfunc(flatten(circuit_parameters),frequency)
   return target_impedance
end