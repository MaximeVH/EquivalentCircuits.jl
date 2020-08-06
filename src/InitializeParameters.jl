"""
    InitializeParameters(Circuit,ranges=Dict('R'=>1000,'C'=>0.001))
Provides a random set of parameters for a given circuit. Upper bounds of the
parameters can be imposed with the range argument.

"""
function InitializeParameters(Circuit,ranges=Dict('R'=>1000,'C'=>0.001,'L'=>1))
    Elements = [Circuit[r] for r in findall(x->x=='R'||x=='C'||x=='L', Circuit)]
    Element_values = Array{Float64}(undef, length(Elements))
    for (n,e) in enumerate(Elements)
            Element_values[n] = ranges[e]*rand()
    end
    return Element_values
end
