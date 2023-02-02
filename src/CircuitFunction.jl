"""
circuitfunction(circuitstring::String)

Generate a function that calculates the impedance of a given string representation of a circuit.

The input is the string representation of a circuit (e.g. "R1-[C2,R3]-P4").
The output is a function that calculates the impedance for the circuit for a given set of parameters and frequency value.

"""
function circuitfunction(circuitstring)
    for (f, t) in zip(["-",  "[",    ","         ,"]"        ],
                      ["+",  "((",   ")^-1+("    ,")^-1)^-1" ])
        circuitstring=replace(circuitstring, f=>t)
    end
    for I in 2:-1:1
        Es = eachmatch(Regex("([CLRPW])([0-9]){$(I)}"),circuitstring)
            for e in Es
                match = e.match
                if match[1] == 'C'
                    circuitstring = replace(circuitstring,match=>"(1/(2im*π*f*"*"T"*"))")
                elseif match[1] == 'R'
                    circuitstring = replace(circuitstring,match=>"T")
                elseif match[1] == 'L'
                    circuitstring = replace(circuitstring,match=>"(2im*π*f*"*"T"*")")
                elseif match[1] == 'P'
                    circuitstring = replace(circuitstring, match => "(1/("*"T"*"*(2im*π*f)^"*"N"*"))") # compatible to package "impedance.py"
                elseif match[1] == 'W'
                    circuitstring = replace(circuitstring,match=>"T*(2*π*f)^(-0.5)"*"*(cos((π*0.5)*0.5)-sin((π*0.5)*0.5)im)")
            end
        end
    end
    new_circuit = ""
    counter = 1
    for i in circuitstring
        if i == 'T'
            new_circuit = new_circuit*"T["*string(counter)*"]"
            counter += 1
        elseif i == 'N'
            new_circuit = new_circuit*"T["*string(counter)*"]"
            counter += 1
        else
            new_circuit = new_circuit*i
        end
    end
    circuit_expression = Meta.parse(new_circuit)
    return mk_function([:T,:f],[],circuit_expression)
end

function circtuitfunction_fixed_params(circuitstring,param_indices,param_values)
    temp_circ_func = circuitfunction(circuitstring)
    function fixed_param_function(parameters,frequency)
        for (idx,val) in zip(param_indices,param_values)
            parameters[idx] = val
        end
        return temp_circ_func(parameters,frequency)
    end
    return fixed_param_function
end

function tree_to_function(tree_array)
    dictionary_of_calculations = Dict('R'=>"T",'C'=>"(1/(2im*π*f*"*"T"*"))",'L'=>"(2im*π*f*"*"T"*")",'P'=>"(1/(T[1]*(2im*π*f)^T[2]))") 
    tree = copy(tree_array) 
    # nodecount = length(tree) # curently not in use
    essential_info = [[node.ParentIndex,node.Type,node.Index] for node in tree]
    for i in 1:((length(essential_info)-1)/2)
        parent1, type1, index1 = pop!(essential_info)
        parent2, type2, index2 = pop!(essential_info) # parent2 currently not in use
        type1 = haskey(dictionary_of_calculations,type1) ? replace(dictionary_of_calculations[type1],"T"=>"T["*string(index1)*']') : type1
        type2 = haskey(dictionary_of_calculations,type2) ? replace(dictionary_of_calculations[type2],"T"=>"T["*string(index2)*']') : type2
        operation = essential_info[parent1][2]
        if operation == '+'
            essential_info[parent1][2] = type1*'+'*type2
        elseif operation == '-'
            essential_info[parent1][2] = "((("*type1*")^-1)+("*type2*")^-1)^-1"
        end
    end
    circuitExpression = Meta.parse(essential_info[1][2])
return mk_function([:T,:f],[],circuitExpression)
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