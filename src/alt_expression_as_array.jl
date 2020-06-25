function alt_expression_as_array(Circuit)
    # Circuit = replace_caps(Circuit_expression2(Circuit
    Circuit = Circuit_expression(Circuit)
    single_matches = eachmatch(r"(C|R)([0-9]){1}",Circuit)
    double_matches = eachmatch(r"(C|R)([0-9]){2}",Circuit)
    for E in double_matches
        Circuit = replace(Circuit,E.match=>"P["*E.match[2:end]*"]")
    end
    for E in single_matches
        Circuit = replace(Circuit,E.match=>"P["*E.match[2:end]*"]")
    end
    # println(Circuit)
    return Meta.parse(Circuit)
end
