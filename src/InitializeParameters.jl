function InitializeParameters(Circuit,ranges=Dict('R'=>1000,'C'=>0.001))
    Elements = [Circuit[r] for r in findall(x->x=='R'||x=='C', Circuit)]
    Element_values = Array{Float64}(undef, length(Elements))
    for (n,e) in enumerate(Elements)
            Element_values[n] = ranges[e]*rand()
    end
    return Element_values
end
