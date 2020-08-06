using ExprRules

abstract type CircuitObject end

mutable struct RC_Circuit <: CircuitObject
    expression::RuleNode
    string::String
    parameters::Array{Float64,1}
end

CircuitGrammar = @grammar begin
D = "-"*B | ""
E = "C"*D | "R"*D
P = "["*B*","*B*"]"*D
B = E | P
end


function CircuitExpression(type,minimum_length = 2)
    generated = false
    Expression = ""
    while generated == false
        try
        Expression = rand(RuleNode,CircuitGrammar,type)
        catch ArgumentError
            continue
        end
        if length(Expression) ≥ minimum_length*3
            generated = true
        end
        # println(Expression)
    end
    return Expression
end

function CircuitExpression()
    generated = false
    Expression = ""
    while generated == false
        try
        Expression = rand(RuleNode,CircuitGrammar,:B)
        catch ArgumentError
            continue
        end
        if length(Expression) ≥ 6
            generated = true
        end
        # println(Expression)
    end
    return Expression
end

function NumberCircuit(Circuit)
    elements = eachmatch(r"([LRC])",Circuit)
    counter = 1
    Numbered_circuit = ""
    prev_ = 1
    start = 0
    for e in elements
        start = e.offset[1]
        Numbered_circuit *= Circuit[prev_:start]*string(counter)
        counter += 1
        prev_ = start + 1
    end
    if start < length(Circuit)
    Numbered_circuit *= Circuit[prev_:end]
    end
    return Numbered_circuit
end

function CircuitString(CircuitExpr,numbered = true)
    exec = get_executable(CircuitExpr,CircuitGrammar)
    string = eval(exec)
    if numbered
        string = NumberCircuit(string)
    end
    return string
end

function CircuitExpression(type,minimum_length = 2)
    generated = false
    Expression = ""
    while generated == false
        try
        Expression = rand(RuleNode,CircuitGrammar,type)
        catch ArgumentError
            continue
        end
        if length(Expression) ≥ minimum_length*3
            generated = true
        end
        # println(Expression)
    end
    return Expression
end

function CircuitExpression()
    generated = false
    Expression = ""
    while generated == false
        try
        Expression = rand(RuleNode,CircuitGrammar,:B)
        catch ArgumentError
            continue
        end
        if length(Expression) ≥ 6
            generated = true
        end
        # println(Expression)
    end
    return Expression
end


function InitializeCircuit(minimum_length=2)
    expression  = CircuitExpression(:B,minimum_length)
    string = CircuitString(expression)
    parameters = InitializeParameters(string)
    return RC_Circuit(expression,string,parameters)
end

function InitializeCircuit(expression::RuleNode)
    string = CircuitString(expression)
    parameters = InitializeParameters(string)
    return RC_Circuit(expression,string,parameters)
end

function InitializePopulation(N=100,minimum_length=3)
    Population = [InitializeCircuit(minimum_length) for i in 1:N]
    return Population
end
