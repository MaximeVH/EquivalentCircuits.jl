using EquivalentCircuits, PlotlyJS

# --- generate frequency vector with n_elements with the same range_data as given in the measurement:
n_elements          = 100
frequ_range         = [1e-1, 1e+4]
frequ_vec           = exp10.(LinRange(log10(frequ_range[1]), log10(frequ_range[2]), n_elements))
b_switch_imag_sign  = false
case_nr = 3

# ---
if case_nr == 1
    # test warburg element:
    circuit_model_preset = "[R1, W2]"
    circuit_params_preset = (R1=1000.0, W2 = 10.0e-6, )
elseif case_nr == 2
    circuit_model_preset = "[R1, P2]"
    circuit_params_preset = (R1=1000.0, Pw2 = 10.0e-6, Pn2 = 10.0e-4, )
elseif case_nr == 3
    circuit_model_preset = "[R1, C2]"
    circuit_params_preset = (R1=1000.0, C2 = 10.0e-6, )
else
    error("Case not defined!")
end

# --- simulate impedance based on suitable preset of a circuit model and its corresponding parameter-set: 
circfunc_preset = EquivalentCircuits.circuitfunction(circuit_model_preset)
impedance_preset = EquivalentCircuits.simulateimpedance_noiseless(circfunc_preset, circuit_params_preset, frequ_vec)

real_part = real(impedance_preset)
# imag_part = - imag(impedance_preset) 
if b_switch_imag_sign
    imag_part = imag(impedance_preset)
else
    imag_part = - imag(impedance_preset)
end

# --- plot: -------------------------------------------------------------------
line_data = AbstractTrace[]
push!(line_data, scatter(;x = real_part, y = imag_part, name = circuit_model_preset, 
    mode = "markers", marker_symbol = "circle-open-dot", marker_size = 10, marker_line_width = 2))   
PlotlyJS.Plot(line_data, PlotlyJS.Layout(; showlegend = true, ) )

