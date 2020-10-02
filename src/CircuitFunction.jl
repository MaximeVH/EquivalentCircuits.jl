"""
    CircuitFunction(Circuit)
Converts the string representation of a Circuit to a function with as inputs the circuit's
parameter values and experimental input frequencies. The resulting outputted function
calculates the impedance value that the circuit would produce with a given set of parameters
and for a given frequency of the input voltage.

"""
function CircuitFunction_inclusive_Final(Circuit)
    #Basic adjustments to get an expression.
    for (f,t) in zip(["-","[",",","]"],["+","((",")^-1+(",")^-1)^-1"])
        Circuit=replace(Circuit,f=>t)
    end
    for I in 2:-1:1
        Es = eachmatch(Regex("([CLRP])([0-9]){$(I)}"),Circuit)
            for e in Es
                match = e.match
                if match[1] == 'C'
                 Circuit = replace(Circuit,match=>"(1/(2im*π*f*"*"T"*"))")
                elseif match[1] == 'R'
                 Circuit = replace(Circuit,match=>"T")
                elseif match[1] == 'L'
                 Circuit = replace(Circuit,match=>"(2im*π*f*"*"T"*")")
             elseif match[1] == 'P'
                 Circuit = replace(Circuit,match=>"T*(2*π*f)^(-N)"*"*(cos((π*N)*0.5)-sin((π*T)*0.5)im)")
            end
        end
end
new_circuit = ""
counter = 1
for i in Circuit
    if i == 'T'
        new_circuit = new_circuit*"T["*string(counter)*"]"
        counter += 1
    elseif i == 'N'
        new_circuit = new_circuit*"T["*string(counter)*"]"
    else
        new_circuit = new_circuit*i
    end
end

    Circuit_expression = Meta.parse(new_circuit)
    return genfun(Circuit_expression,[:T,:f])
end
