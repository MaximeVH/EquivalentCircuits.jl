function Circuit_expression(Circuit)
    # Circuit = deepcopy(CircuitString)
    Circuit2=replace(Circuit,"-"=>"+")
    Circuit3=replace(Circuit2,"["=>"((")
    Circuit4=replace(Circuit3,","=>")^-1+(")
    Circuit5=replace(Circuit4,"]"=>")^-1)^-1")
    # Cs = eachmatch(r"(C)([0-9]{2})+",Circuit5)
    #     for c in Cs
    #         Circuit5 = replace(Circuit5,c.match=>"(1/(2im*π*f*"*c.match*"))")
    #     end
    Cs = eachmatch(r"(C)([0-9])+",Circuit5)
        for c in Cs
            Circuit5 = replace(Circuit5,c.match=>"(1/(2im*π*f*"*c.match*"))")
        end
    return Circuit5
end
