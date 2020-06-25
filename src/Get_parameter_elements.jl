function Get_parameter_elements(Circuit)
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
return Elements
end
