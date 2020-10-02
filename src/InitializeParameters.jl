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
