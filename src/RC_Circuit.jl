abstract type Circuit end

mutable struct LRC_Circuit <: Circuit
    karva::String
    circuit::String
    parameters::Array{Any,1}
    parameter_indices::Array{Int,1}
end

function Karva_to_ET(karva)
    level = 1
    nonterminal = true
    idx = 1
    param_indexes = 1:length(karva)
    total_arity = arity(karva[1])
    levels = []
    index_levels = []
    push!(levels,karva[1])
    push!(index_levels,param_indexes[1])
    while total_arity > 0
        push!(levels,karva[idx+1:idx+total_arity])
        push!(index_levels,param_indexes[idx+1:idx+total_arity])
        level += 1
        idx = idx + total_arity
        total_arity  = sum([arity(k) for k in levels[level]])
    end
    ET = "l1"
    ET_params = "l1"
    for (lev,param) in zip(levels,index_levels)
        for (i,(o,pa)) in enumerate(zip(lev,param))
            if is_operation(o)
                ET = replace(ET,"l"*string(i) => o*"[l,l]")
                ET_params = replace(ET_params,"l"*string(i) => o*"[l,l]")
            else
                ET = replace(ET,"l"*string(i) => o)
                ET_params = replace(ET_params,"l"*string(i) => string(pa))
            end
        end
        ET = number_ls(ET)
        ET_params = number_ls(ET_params)
    end
    parameter_indices = parse.(Int,split(replace(ET_params,r"[\]\[+-]"=>""),','))
    return ET, parameter_indices
end

function Generate_karva(operators,terminals,head)
    all_elements = operators*terminals
    karva = rand(operators)*String(rand(all_elements,head-1))*String(rand(terminals,head+1))
    return karva
end

Generate_circuit_karva(head) = Generate_karva("+-","LCRP",head)

function extract_deepest_part(Tree)
    depth = ET_depth(Tree)
    deepest_start = findall(x->x=='[',Tree)[depth]-1
    deepest_end = findall(x->x==']',Tree[deepest_start:end])[1]
    return Tree[deepest_start:deepest_start+deepest_end-1]
end

function deepest_part_to_expression(part)
    operation = part[1]
    seperator = findfirst(x->x==',',part)
    argument_1 = part[3:seperator-1]
    argument_2 = part[seperator+1:end-1]
    if operation == '+'
        return argument_1*'-'*argument_2
    elseif operation == '-'
        return '('*argument_1*'.'*argument_2*')'
    end
end


function ET_to_expression(Tree)
    depth = ET_depth(Tree)
    # for i in 1:depth
    while occursin('[',Tree)
        deepest_part = extract_deepest_part(Tree)
        deepest_part_expression = deepest_part_to_expression(deepest_part)
        Tree = replace(Tree,deepest_part=>deepest_part_expression)
    end
    Tree = replace(Tree,'('=>'[')
    Tree = replace(Tree,')'=>']')
    Tree = replace(Tree,'.'=>',')
    return Tree
end


function ET_depth(ET)
    max_brac = 0
    brac = 0
    for chr in ET
        if chr == '['
            brac += 1
        elseif chr == ']'
            brac -= 1
        end
        if brac > max_brac
            max_brac = brac
        end
    end
    return max_brac
end

function arity(chr)
    if chr in "/*-+"
        return 2
    else
        return 0
    end
end

function number_ls(ET)
    counter = 1
    new_ET = ""
    for i in ET
        if i == 'l'
            new_ET*="l"*string(counter)
            counter += 1
        else
            new_ET*=i
        end
    end
    return new_ET
end

is_operation(chr) = occursin(chr,"-+")

function number_ET(ET)
    numbered_ET = ""
    counter = 1
    for i in ET
        if occursin(i,"LRCP")
            numbered_ET *= i*string(counter)
            counter += 1
        else
            numbered_ET *= i
        end
    end
    return numbered_ET
end


function denumber(expression)
    new_expression = ""
    for i in expression
        if occursin(i,"LRCP+-[],")
            new_expression *= i
        end
    end
    return new_expression
end

function extract_elements(circuit)
    new_expression = ""
    for i in circuit
        if occursin(i,"LRCP")
            new_expression *= i
        end
    end
    return new_expression
end

function circuit_depth(circuit)
    max_depth = 0
    bracket_counter = 0
    for i in circuit
        if i == '['
            bracket_counter += 1
            if bracket_counter > max_depth
                max_depth = bracket_counter
            end
        elseif i == ']'
            bracket_counter -= 1
        end
    end
    return max_depth
end

macro genfun(expr,args...); :(($(args...),)->$expr) end

genfun(expr,args::Union{Vector,Tuple}) = eval(:(($(args...),)->$expr))
genfun(expr,args::Symbol...) = genfun(expr,args)


function Generate_circuit(max_depth,min_elements,terminals)
    karva = "" ; ET = "" ; parameter_inds = ""
    n_elements = 0
    while n_elements < min_elements
        karva = Generate_karva("+-",terminals,max_depth)
        ET,parameter_inds = Karva_to_ET(karva)
        n_elements = count(x->occursin(x,"PLRC"),ET)
    end
    parameters = karva_parameters_inclusive(karva)
    circuit = ET_to_expression(number_ET(ET))
    return LRC_Circuit(karva,circuit,parameters,parameter_inds)
end


function Generate_population(N=20,max_depth=8,min_elements=3,terminals = "LRCP")
    return [Generate_circuit(max_depth,min_elements,terminals) for i in 1:N]
end

function to_flat_array(params)
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

function deflatten_parameters(circuit,parameters)
    elements = get_elements(circuit)
    index_counter = 1
    new_parameters = Array{Any}(undef, length(elements))
    for (n,e) in enumerate(elements)
        if e != 'P'
            new_parameters[n] = parameters[index_counter]
            index_counter += 1
        else
            new_parameters[n] = [parameters[index_counter],parameters[index_counter+1]]
            index_counter += 2
        end
    end
    return new_parameters
end

function get_elements(circuit)
    elements = replace(circuit,r"[-,\[\]0-9]"=>"")
    return elements
end

function Nyquist(circuit,frequencies;title=false)
    CircuitFunc = CircuitFunction_inclusive_Final(circuit.circuit)
    parameters = to_flat_array(circuit.parameters[circuit.parameter_indices])
    model_output = [Base.invokelatest(CircuitFunc,parameters,fr) for fr in frequencies]
    reals = real(model_output)
    imags = imag(model_output)
    if title
        scatter(reals,-imags,title = circuit.circuit)
    else
        scatter(reals,-imags)
    end
end
