# --------------------------------------------------------------------------------------------------------------------------
# --- load packages after having changed to the activated environment:
using PyCall, Printf, PlotlyJS, RobustModels, JLD2
using EquivalentCircuits
import Pkg, XLSX, FileIO
circuits = PyCall.pyimport("impedance.models.circuits")

# --- transfer EquivalentCircuits_jl-syntax ==> impdance_py-syntax
function _MylibExpCircuitStrToImpPy(_str_in::String)
    left_round_brackets =  replace(_str_in, "[" => "(")
    round_brackets =  replace(left_round_brackets, "]" => ")")
    char_P_replaced = replace(round_brackets, "P" => "CPE") # P = constant phase element
    return replace(char_P_replaced, "(" => "p(")
end

# --- frequencies:
frequ_data_all = [0.1, 0.20691380811147891, 0.42813323987193935, 0.8858667904100828, 1.8329807108324356, 3.792690190732248, 7.847599703514613, 16.237767391887218,
33.59818286283781, 69.51927961775601, 143.84498882876616, 297.63514416313194, 615.8482110660267, 1274.2749857031336, 2636.650898730358, 5455.5947811685255, 11288.378916846883,
23357.214690901215, 48329.30238571752, 100000.0]

# --- parameters: ----------------------------------------------------------------------------------------------------
case_nr = 4
# ---
if case_nr == 1
    s_label_EquivCirc     = "SimpleRC"
    circuit_model_preset  = "[C1,R2-[R3,C4]]"
    circuit_params_preset = (C1 = 0.025036871360843482, R2 = 396.73873944116787, R3 = 2178.061127814435, C4 = 1.1589755057609664e-5)
elseif case_nr == 2
    s_label_EquivCirc     = "Matthias"
    circuit_model_preset  = "[C1-P2,[R3,C4-L5-P6]]"
    circuit_params_preset = (C1 = 1.761822019169847, P2w = 3.7508155012864477, P2n = 0.24165725697257198, R3 = 0.01897566655203251, C4 = 0.12994316115177285, L5 = 6.365197784775357e-7, P6w = 	17.49835814269397, 	P6n = 0.11636029150118028)
elseif case_nr == 3
    s_label_EquivCirc     = "SimpleCPE"
    circuit_model_preset  = "[C1,P2]"
    circuit_params_preset = (C1 = 1.761822019169847, P2w = 3.7508155012864477, P2n = 0.24165725697257198)
elseif case_nr == 4
    s_label_EquivCirc     = "SimpleCPE"
    circuit_model_preset  = "R1-[C2,P3]"
    circuit_params_preset = (R1 = 321.124, C2 = 1.761, P3w = 3.75, P3n = 0.24)
elseif case_nr == 5
    s_label_EquivCirc     = "SimpleCL"
    circuit_model_preset  = "[C1,L2]"
    circuit_params_preset = (C1 = 1.761822019169847, L2 = 1.7508155012864477)
elseif case_nr == 6
    s_label_EquivCirc     = "ECircJL_Expl"
    circuit_model_preset  = "R1-[C2,R3-[R4,C5]]"
    circuit_params_preset = (R1 = 19.999999999999996, C2 = 4.0e-9, R3 = 3400.0, R4 = 2500.0, C5 = 4.0e-6)
else
    error("Case not defined!")
end

# --- simulate via ImpPy (impedance.py) based on initual guess / set of equivalent circuit parameters:
circuitStr_ImpPy = _MylibExpCircuitStrToImpPy(circuit_model_preset)
param_ImpPy = collect(circuit_params_preset)
equiv_circuit_PyObj = circuits.CustomCircuit(initial_guess= param_ImpPy, circuit= circuitStr_ImpPy, name = s_label_EquivCirc)
Z_simulated_ImpPy = equiv_circuit_PyObj.predict(frequ_data_all, use_initial = true)


# --- simulate via ECircJL (EquivalentCircuits.jl): ---------------------------------------------------
circfunction = EquivalentCircuits.circuitfunction(circuit_model_preset)
# circfunction = EquivalentCircuits.circuitfunction(circuit_model_preset)

Z_simulated_ECircJL = EquivalentCircuits.simulateimpedance_noiseless(circfunction, circuit_params_preset, frequ_data_all);

Q_ImpPy_ECircJL = RobustModels.mean(abs.(Z_simulated_ImpPy - Z_simulated_ECircJL))
println(" Δ(Q_ImpPy - Q_ECircJL): ", Q_ImpPy_ECircJL)

# --- plot nyquist: --------------------------------------------------------------------------------------------------------
function plot_nyquist(_Circuit_str::String, _frequ_data::Vector{<:Number}, _Z_ImpPy::Vector{ComplexF64}, _Z_ECircJL::Vector{ComplexF64})
    s_pts_info = Vector{String}(undef, 0)
    for i_ndx in eachindex(_frequ_data)
        push!(s_pts_info, @sprintf("#:%i, f=%.2fHz", i_ndx, _frequ_data[i_ndx]))
    end
    # ---
    trace_ImpPy   = PlotlyJS.scatter(; x= real(_Z_ImpPy),       y= imag(_Z_ImpPy),       name = "ImpPy",
                                    text = s_pts_info, mode = "markers")
    trace_ECircJL = PlotlyJS.scatter(; x= real(_Z_ECircJL),     y= imag(_Z_ECircJL),     name = "ECircJL")
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
        title_text          = "Comparison Impedance Simulation: ECircJL <-> ImpPy",
        xaxis_title_text    = "z<sub>Real</sub>",
        xaxis_dtick         = 1000,
        yaxis_title_text    = "z<sub>Imag</sub>",
        yaxis_dtick         = 1000,
        # --- Fixed Ratio Axes Configuration
        yaxis_scaleanchor   = "x",
        yaxis_scaleratio    = 1,
        annotations         = _annotations,
    )
    return PlotlyJS.Plot([trace_ImpPy, trace_ECircJL], plt_layout)
end

hdl_plt = plot_nyquist(circuit_model_preset, frequ_data_all, Z_simulated_ImpPy, Z_simulated_ECircJL)
PlotlyJS.display(hdl_plt)

