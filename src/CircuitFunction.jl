function CircuitFunction(Circuit)
    for (f,t) in zip(["-","[",",","]"],["+","((",")^-1+(",")^-1)^-1"])
    Circuit=replace(Circuit,f=>t)
    end
    Cs = eachmatch(r"(C)([0-9])+",Circuit)
        for c in Cs
            Circuit = replace(Circuit,c.match=>"(1/(2im*π*f*"*c.match*"))")
        end
    double_matches = eachmatch(r"(C|R)([0-9]){2}",Circuit)
    for E in double_matches
        Circuit = replace(Circuit,E.match=>"P["*E.match[2:end]*"]")
    end
    single_matches = eachmatch(r"(C|R)([0-9]){1}",Circuit)
    for E in single_matches
        Circuit = replace(Circuit,E.match=>"P["*E.match[2:end]*"]")
    end

    Circuit_expression = Meta.parse(Circuit)

    return genfun(Circuit_expression,[:P,:f])
end
