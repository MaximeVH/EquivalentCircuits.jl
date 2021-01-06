"""
    InitializeParameters(Circuit,ranges=Dict('R'=>1000,'C'=>0.001))
Provides a random set of parameters for a given circuit. Upper bounds of the
parameters can be imposed with the range argument.

"""

function InitializeParameters_inclusive(Circuit,ranges=Dict('R'=>1000,'C'=>0.001,'L'=>1,'P'=>[100,1]))
    Elements = [Circuit[r] for r in findall(x->x=='R'||x=='C'||x=='L'||x=='P', Circuit)]
    Element_values = Array{Any}(undef, length(Elements))
    for (n,e) in enumerate(Elements)
            if e == 'P'
                Element_values[n] = [ranges[e][1]*rand(),ranges[e][2]*rand()]
            else
                Element_values[n] = ranges[e]*rand()
            end
    end
    return Element_values
end

function karva_parameters_inclusive(karva)
    parameters = Array{Any}(undef, length(karva))
    ranges = Dict('R'=>1000,'C'=> 0.001,'L'=> 1,'+'=> 0,'-'=> 0,'P'=>[100,1])
    for (e,i) in enumerate(karva)
            if e == 'P'
                Element_values[n] = [ranges[e][1]*rand(),ranges[e][2]*rand()]
            else
            parameters[e] = rand()*ranges[i]
        end
    end
    return parameters
end

function InitializeParameters(karva)
    ranges=Dict('a'=>1000,'b'=>0.001,'c'=>1,'d'=>[100,1])
    params = [] #preallocate this.
    for k in karva
        if k in ['a','b','c']
            push!(params,rand()*ranges[k])
        elseif k == 'd'
            push!(params,[rand()*ranges['d'][1],rand()*ranges['d'][2]])
        elseif k in ['e','h'] #RC
            push!(params,[rand()*ranges['a'],rand()*ranges['b']])
        elseif k in ['f','i'] #RL
            push!(params,[rand()*ranges['a'],rand()*ranges['c']])
        elseif k in ['k','m'] #CL
            push!(params,[rand()*ranges['b'],rand()*ranges['c']])
        elseif k in ['g','j'] #RP 
            push!(params,[rand()*ranges['a'],[rand()*ranges['d'][1],rand()*ranges['d'][2]]])
        elseif k in ['l','n'] #CP 
            push!(params,[rand()*ranges['b'],[rand()*ranges['d'][1],rand()*ranges['d'][2]]])
        elseif k == 'p' #LP 
            push!(params,[rand()*ranges['c'],[rand()*ranges['d'][1],rand()*ranges['d'][2]]])
        elseif k in ['+','-']
            push!(params,0)
        end
    end
    return params
end