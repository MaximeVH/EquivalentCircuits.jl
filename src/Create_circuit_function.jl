function Create_Circuit_function(Circuit)
    expression = alt_expression_as_array(Circuit)
    c_func = genfun(expression,[:P,:f])
    return c_func
end

# f(x) = 5x + 2
#
# f(5)
#
# @eval(f(5))

# function Evaluate_Circuit_function(Circuit)
