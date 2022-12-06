# --- relevant adresses ----------------------------------------------------------------------------------
# || docu: https://impedancepy.readthedocs.io/en/latest/getting-started.html#step-1-installation
# || GitHub-Repository: https://github.com/ECSHackWeek/impedance.py
# --- installation: --------------------------------------------------------------------------------------
# || julia> using Conda; Conda.add("impedance")
# --- installation via pip:
# --------------------------------------------------------------------------------------------------------------------------
using PyCall, Printf, PlotlyJS, RobustModels, EquivalentCircuits, JLD2
circuits = PyCall.pyimport("impedance.models.circuits")
np = PyCall.pyimport("numpy")
pyplt = PyCall.pyimport("matplotlib.pyplot")

#> --- define mutuable struct for equivalent circuits
mutable struct _MyEquivCircStruct
    CircuitStr :: String
    CircuitParams_initial :: NamedTuple
    Quality_initial :: Float64
    CircuitParams_optimized :: Vector{Float64}
    Quality_optimized :: Float64
end

function _MyLibIndxCircuitStr(_in_CircStr::String)
    _indxCircuitStr = 0
    for i_ = eachindex(equiv_circ_vec)
        # println(typeof(i_))
        if _in_CircStr == equiv_circ_vec[i_].CircuitStr
            _indxCircuitStr = i_
        end
    end
    return _indxCircuitStr
end

function _MyLibFindfirstempty(_in::Vector{_MyEquivCircStruct})
    _indx_firstempty = nothing
    for i_ = eachindex(_in)
        if isassigned(equiv_circ_vec, i_) && isempty(_in[i_].CircuitStr)
            _indx_firstempty = i_
        end
    end
    return _indx_firstempty
end


function _MyLibSortEuivCircuit(_vec_in::Vector{_MyEquivCircStruct})
    # --- sort results:
    _vec_Q = Vector{Float64}(undef, n_best)
    for i_ = 1:n_best
        _vec_Q[i_] = _vec_in[i_].Quality_initial
    end
    _sort_idx = sortperm(_vec_Q)
    return _vec_in[_sort_idx]        
end

function _MylibExpCircuitStr(_str_in::String)
    left_round_brackets =  replace(_str_in, "[" => "(")
    round_brackets =  replace(left_round_brackets, "]" => ")")
    return replace(round_brackets, "(" => "p(")
end

# --- simulate impedance, based on given equivalent-circuit model via own function:
function _simulate_impedance(_in_frequ_vec::Vector{<:Number}, _circuit_str::AbstractString, _parameters_in::NamedTuple)
    _circfunc = EquivalentCircuits.circuitfunction(_circuit_str)
    return EquivalentCircuits.simulateimpedance_noiseless(_circfunc, _parameters_in, _in_frequ_vec)
end

function _build_NamedTuple(_circuit_model_preset::String, _param_vec::Vector{<:Number})
    _r_search_pat = r"[RLC]{1}[1-9]+"
    _vec_param_symbol = Vector{Symbol}(undef, 0)
    idx_start = 1
    for i_ = 1:length(_circuit_model_preset) - 1
        # global idx_start
        match_result    = match(_r_search_pat, _circuit_model_preset, idx_start)
        if isnothing(match_result)
            break
        else
            push!(_vec_param_symbol, Symbol(match_result.match))
            _, idx_start = findfirst(match_result.match, _circuit_model_preset) 
            idx_start += 1
            # println("match: ", match_result.match, ", idx_start: ", idx_start, ", str[idx_start:end]: ", circuit_model_preset[idx_start:end])
        end
    end
    # ---
    return NamedTuple.(zip.(Ref(_vec_param_symbol), zip(_param_vec...)))
end


# --- measured data: -------------------------------------------------------------------------------------------------------
# origin of data: https://github.com/MaximeVH/EquivalentCircuits.jl/blob/master/example_measurements.csv
# --- frequencies:
frequ_data = [0.1, 0.20691380811147891, 0.42813323987193935, 0.8858667904100828, 1.8329807108324356, 3.792690190732248, 7.847599703514613, 16.237767391887218, 
33.59818286283781, 69.51927961775601, 143.84498882876616, 297.63514416313194, 615.8482110660267, 1274.2749857031336, 2636.650898730358, 5455.5947811685255, 11288.378916846883, 
23357.214690901215, 48329.30238571752, 100000.0]

# --- sample impedance data:
Z_data = ComplexF64[5919.90084073586 - 15.794826681804063im, 5919.575521325405 - 32.677443741883025im, 5918.183674897797 - 67.57666460870544im, 
5912.242152808868 - 139.49441081175667im, 5887.119965779756 - 285.73699600024963im, 5785.038233646888 - 566.878749499386im, 5428.935296370544 - 997.1881947423717im, 
4640.2144606930815 - 1257.8277219098052im, 3871.8361085209845 - 978.9656717819273im, 3537.682636142598 - 564.9627167404748im, 3442.9419240480606 - 315.3996363823805im, 
3418.140460121871 - 219.68986957046025im, 3405.513645508888 - 242.57272592660013im, 3373.904450003017 - 396.0671717029891im, 3249.673719526995 - 742.0301829777005im, 
2808.423185495856 - 1305.924162464249im, 1779.4087896271944 - 1698.9660879948128im, 701.9588433822435 - 1361.4674351816855im, 208.28978681589047 - 777.6453690080142im, 
65.93273498232111 - 392.50667235657im]

# --- generate frequency vector with n_elements with the same range as given in the measurement:
n_elements = 100
frequ_vec = exp10.(LinRange(log10(frequ_data[1]), log10(frequ_data[end]), n_elements))

# --- more general parameters: ---------------------------------------------------------------------------------------------
b_run_EC_evolution = true
nf_equ_circ_evo = raw"C:\tmp\data_log\equvialent_circuit_evolution_results.jld"
n_best          = 15
n_test          = 3
Q_th            = 250

# --- find suitable equivalent circuits via "EquivalentCircuits": -------------------------------------------------------------------------------
# **************************************************************************************************************************
# --- function "circuitevolution()" to find a suitable equivalent circuit model and its parameters:
# **************************************************************************************************************************
# # Arguments
# - `generations::Integer=10`: the maximum number of iterations of the evolutionary algorithm.
# - `population_size::Integer=30`: the number of individuals in the population during each iteration.
# - `terminals::String="RCLP"`: the circuit components that are to be included in the circuit identification.
# - `head::Integer=8`: a hyperparameter than controls the maximum considered complexity of the circuits.
# - `cutoff::Float64=0.8`: a hyperparameter that controls the circuit complexity by removing redundant components.
# - `initial_population::Array{Circuit,1}=nothing`:the option to provide an initial population of circuits
#   (obtained by using the loadpopulation function) with which the algorithm starts.
# -------------------------------------------------------------------------------------------------------------------------
terminals_      = "RC"
head_           = 6

if isfile(nf_equ_circ_evo)
    equiv_circ_vec = load(nf_equ_circ_evo, "equiv_circ_vec")
else
    equiv_circ_vec = Vector{_MyEquivCircStruct}(undef, n_best)
    for i_ = 1:n_best
        equiv_circ_vec[i_] = _MyEquivCircStruct("", (R = 1.0, B = 0.002, R2 = 2.43),   +Inf, [], +Inf)
    end
end

# ---
if b_run_EC_evolution || ~isfile(nf_equ_circ_evo)
    for i_ = 1:n_test
        global equiv_circ_vec
        println("#: ", i_)
        println("  DBG#: 1") # sometimes I have observed chrashes.
        i_equiv_circ_evo  = EquivalentCircuits.circuitevolution(Z_data, frequ_data, terminals=terminals_, generations=25, population_size=100, head=head_)
        println("  DBG#: 2")
        # --- simulate Impedance:
        _circfunc_evo    = EquivalentCircuits.circuitfunction(i_equiv_circ_evo.Circuit)
        println("  DBG#: 3")
        # --- Calc quality as the mean of the distances between measured and simulated impedance:
        _impedance_evo_data_pts = EquivalentCircuits.simulateimpedance_noiseless(_circfunc_evo, i_equiv_circ_evo.Parameters, frequ_data)
        println("  DBG#: 4")
        i_Q              = RobustModels.mean(abs.(Z_data - _impedance_evo_data_pts))
        # ---
        idx_CiruitStr   = _MyLibIndxCircuitStr(i_equiv_circ_evo.Circuit)
        idx_firstEmpty  = _MyLibFindfirstempty(equiv_circ_vec)
        println("  DBG, 5, idx_CiruitStr: ", idx_CiruitStr, ", idx_firstEmpty: ", idx_firstEmpty)
        if  idx_CiruitStr == 0 && ~isnothing(idx_firstEmpty) && i_Q < Q_th  
            equiv_circ_vec[idx_firstEmpty] = _MyEquivCircStruct(i_equiv_circ_evo.Circuit, i_equiv_circ_evo.Parameters, i_Q, [], +Inf)
        elseif (idx_CiruitStr != 0) && i_Q < Q_th
            println(i_equiv_circ_evo.Circuit)
            if i_Q < equiv_circ_vec[idx_CiruitStr].Quality_initial 
                equiv_circ_vec[idx_CiruitStr] = _MyEquivCircStruct(i_equiv_circ_evo.Circuit, i_equiv_circ_evo.Parameters, i_Q, [], +Inf)
            end
            # equiv_circ_vec = _MyLibSortEuivCircuit(equiv_circ_vec)
            jldsave(nf_equ_circ_evo, true; equiv_circ_vec)
        else
            println("B1: ", idx_CiruitStr == 0, ", B2: ",  i_Q < Q_th, ", B3: ", idx_CiruitStr != 0, ", i_Q: ", i_Q)
            if i_ == n_test
                jldsave(nf_equ_circ_evo, true; equiv_circ_vec)
            end
        end
        println("  DBG#: 6")
    end
end
# --- load / re-load results and sort:
if isfile(nf_equ_circ_evo)
    equiv_circ_vec = load(nf_equ_circ_evo, "equiv_circ_vec")
else
    error("\"", nf_equ_circ_evo, "\" does not exist!")
end
equiv_circ_vec = _MyLibSortEuivCircuit(equiv_circ_vec)

# --- plot results:
for i_ = 1:n_best
    println(@sprintf("%30s, \tQ: %9.6f", equiv_circ_vec[i_].CircuitStr, equiv_circ_vec[i_].Quality_initial) )
end
println('-'^100, "\n")

# --- optimization:
i_better_th = 0
for i_ = 1:n_best
    global i_better_th, equiv_circ_vec
    if equiv_circ_vec[i_].Quality_initial < Q_th || equiv_circ_vec[i_].Quality_optimized < Q_th
        i_better_th += 1
        i_circuit_str = _MylibExpCircuitStr(equiv_circ_vec[i_].CircuitStr)
        if equiv_circ_vec[i_].Quality_initial < equiv_circ_vec[i_].Quality_optimized
            i_initial_guess_ = collect(equiv_circ_vec[i_].CircuitParams_initial)
        else
            i_initial_guess_ = equiv_circ_vec[i_].CircuitParams_optimized
        end
        i_customCircuit = circuits.CustomCircuit(initial_guess= i_initial_guess_, circuit= i_circuit_str)
        i_customCircuit.fit(frequ_data, Z_data)
        i_simulated_on_data_pts = i_customCircuit.predict(frequ_data)
        i_Q_validation = RobustModels.mean(abs.(Z_data - i_simulated_on_data_pts))
        equiv_circ_vec[i_] = _MyEquivCircStruct(equiv_circ_vec[i_].CircuitStr, equiv_circ_vec[i_].CircuitParams_initial, equiv_circ_vec[i_].Quality_initial, 
                                i_customCircuit.parameters_, i_Q_validation)
    end
end


# --- 
_vec_Q = Vector{Float64}(undef, n_best)
for i_ = 1:n_best
    _vec_Q[i_] = equiv_circ_vec[i_].Quality_optimized
end
_sort_idx = sortperm(_vec_Q)

# --- print results:
for i_ = 1:n_best
    i_idx = _sort_idx[i_]
    println(@sprintf("%-30s, Q_init: %9.6g, Q_opt: %9.6g", equiv_circ_vec[i_idx].CircuitStr, equiv_circ_vec[i_idx].Quality_initial, equiv_circ_vec[i_idx].Quality_optimized) )
end

# --- get best result:
idx_best = _sort_idx[1]
strg_Circuit_best = equiv_circ_vec[idx_best].CircuitStr
equiv_par_nt = _build_NamedTuple(strg_Circuit_best, equiv_circ_vec[idx_best].CircuitParams_optimized)
Z_simulated_data = _simulate_impedance(frequ_data, strg_Circuit_best, equiv_par_nt)
Z_simulated_vec  = _simulate_impedance(frequ_vec, strg_Circuit_best, equiv_par_nt)


# --- plot nyquist: --------------------------------------------------------------------------------------------------------
function plot_nyquist(_Circuit_str::String, _Z_simulated_data::Vector{ComplexF64}, _Z_simulated_vec::Vector{ComplexF64})
    s_pts_info = Vector{String}(undef, 0)
    for i_ndx in eachindex(frequ_data)
        push!(s_pts_info, @sprintf("#:%i, f=%.2fHz", i_ndx, frequ_data[i_ndx]))
    end
    # ---
    trace_impedance = PlotlyJS.scatter(; x= real(Z_data), y=  imag(Z_data), name = "measured", text = s_pts_info, mode = "markers")
    trace_simulated_on_data_pts  = PlotlyJS.scatter(; x= real(_Z_simulated_data), y= imag(_Z_simulated_data), name = "simulated",
                                    text = s_pts_info, mode = "markers")
    trace_simulated_on_frequ_vec = PlotlyJS.scatter(; x= real(_Z_simulated_vec),    y= imag(_Z_simulated_vec),    name = "fine")
    # --- annotations:
    _annotations = [PlotlyJS.attr(; 
    xref        ="paper",           yref    = "paper",
    x           = 0.1,              y       = 0.9, 
    text        = _Circuit_str,
    xanchor     = "left",
    yanchor     = "top",
    font_family = "PT Sans Narrow, monospace",
    showarrow   = false)]
    # ---
    plt_layout = PlotlyJS.Layout(;
        title_text          = "Sample <-> Optimization",
        xaxis_title_text    = "z<sub>Real</sub>",
        xaxis_dtick         = 1000,
        yaxis_title_text    = "z<sub>Imag</sub>",
        yaxis_dtick         = 1000,
        # --- Fixed Ratio Axes Configuration
        yaxis_scaleanchor   = "x",
        yaxis_scaleratio    = 1,
        annotations         = _annotations,
    )
    return PlotlyJS.Plot([trace_impedance, trace_simulated_on_data_pts, trace_simulated_on_frequ_vec], plt_layout)
end

plot_nyquist(strg_Circuit_best, Z_simulated_data, Z_simulated_vec)
