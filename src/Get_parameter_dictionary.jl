function Get_parameter_dictionary(Circuit,ranges = 'n')
    Circuit *= '-'
    # Future improvement: additional argument (with reasonable defaults)
    # to define ranges for the parameters of the different circuit elements.
                    # resistances = findall(x->x=='R', Circuit)
                    # capacitances = findall(x->x=='C', Circuit)
                    # elements = findall(x->x=='R'||x=='C', Circuit)
    # Other elements to be included here later on.
    Elements = [Circuit[r:r+2] for r in findall(x->x=='R'||x=='C', Circuit)]
    for (i,j) in enumerate(Elements)
        if 47 < Int(j[end]) < 57
            continue
        else
            Elements[i] = j[1:end-1]
        end
    end
    Element_values = Dict()
    for e in Elements
        if ranges == 'n'
        Element_values[e] = rand()
        else
            Element_values[e] = ranges[e[1]]*rand()
    end
    # This will have to be adjusted according to the reasonable range for parameter values, as well as
    # the inclusion of the specification of the element type.
end
    return Element_values
end
