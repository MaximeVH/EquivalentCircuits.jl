"""
    CircuitFunction(Circuit)
Converts the string representation of a Circuit to a function with as inputs the circuit's
parameter values and experimental input frequencies. The resulting outputted function
calculates the impedance value that the circuit would produce with a given set of parameters
and for a given frequency of the input voltage.

"""
function CircuitFunction(Circuit)
    #Basic adjustments to get an expression.
    for (f,t) in zip(["-","[",",","]"],["+","((",")^-1+(",")^-1)^-1"])
        Circuit=replace(Circuit,f=>t)
    end
    for I in 2:-1:1
        Es = eachmatch(Regex("([CLR])([0-9]){$(I)}"),Circuit)
            for e in Es
                match = e.match
                if match[1] == 'C'
                 Circuit = replace(Circuit,match=>"(1/(2im*π*f*"*"P["*match[2:end]*"]"*"))")
                elseif match[1] == 'R'
                 Circuit = replace(Circuit,match=>"P["*match[2:end]*"]")
                elseif match[1] == 'L'
                 Circuit = replace(Circuit,match=>"(2im*π*f*"*"P["*match[2:end]*"]"*")")
            end
        end
end
    Circuit_expression = Meta.parse(Circuit)
    return genfun(Circuit_expression,[:P,:f])
end
