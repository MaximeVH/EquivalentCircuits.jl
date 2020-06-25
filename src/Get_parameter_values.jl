function Get_parameter_values(Circuit,ranges)
    param_dict = Get_parameter_dictionary(Circuit,ranges)
    elements = Get_parameter_elements(Circuit)
    values = [param_dict[e] for e in elements]
    return values#,elements
end
