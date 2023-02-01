using  EquivalentCircuits, PlotlyJS
using Optim, BlackBoxOptim

# --- sample frequencies:
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

# --- Gamry's fitted circuit parameters:

circ_strg_ref = "R1-L2-[P3,R4]-[P5,R6]-[P7,R8]"
const R1_ref    = 0.007031                                          # Gamry: HFR
const L2_ref    = 0.00000004257                                     # Gamry: Lstray
# ---
const P3w_ref   = 149.9;            const P3n_ref   = 0.9763        # Gamry: Q_3, a_3
const R4_ref    = 0.00132                                           # Gamry: R_3
# ---
const P5w_ref   = 1.948;            const P5n_ref   = 0.7817        # Gamry: Q_2, a_2
const R6_ref    = 0.009341                                          # Gamry: R2
# ---
const P7w_ref   = 0.2224;           const P7n_ref   = 9.97E-01      # Gamry: Q_1, a_1
const R8_ref    = 0.001118                                          # Gamry: R_1
gamry_params = [R1_ref,L2_ref,P3w_ref,P3n_ref,R4_ref,P5w_ref,P5n_ref,R6_ref,P7w_ref,P7n_ref,R8_ref]
gamry_trace = simulateimpedance_noiseless(circuitfunction(circ_strg_ref),gamry_params,frequ_data)

# Tested optimisation methods from BlackBoxOptim package.

bbopt_methods =  [:separable_nes,:xnes,:dxnes,:adaptive_de_rand_1_bin,:adaptive_de_rand_1_bin_radiuslimited,
:de_rand_1_bin,:de_rand_1_bin_radiuslimited,:de_rand_2_bin,:de_rand_2_bin_radiuslimited,:generating_set_search,
:probabilistic_descent,:resampling_memetic_search,:resampling_inheritance_memetic_search,
:simultaneous_perturbation_stochastic_approximation]

bbopt_method_names =  ["separable_nes","xnes", "dxnes","adaptive_de_rand_1_bin","adaptive_de_rand_1_bin_radiuslimited",
"de_rand_1_bin","de_rand_1_bin_radiuslimited","de_rand_2_bin","de_rand_2_bin_radiuslimited","generating_set_search",
"probabilistic_descent","resampling_memetic_search","resampling_inheritance_memetic_search",
"simultaneous_perturbation_stochastic_approximation"]

# Function used to optimise the circuit parameters, using a method from BlackBoxOptim, followed by Nelder--Mead simplex fine-tuning.
function parameopt(circuitstring::String,measurements,frequencies,method)
    elements = foldl(replace,["["=>"","]"=>"","-"=>"",","=>""],init = denumber_circuit(circuitstring))
    initial_parameters = flatten(karva_parameters(elements));
    circfunc = circuitfunction(circuitstring)
    objective = objectivefunction(circfunc,measurements,frequencies) 
    lower = zeros(length(initial_parameters))
    upper = get_parameter_upper_bound(circuitstring)

    ### First step ###
    SR = Array{Tuple{Float64,Float64},1}(undef,length(initial_parameters))
    for (e,(l,u)) in enumerate(zip(lower,upper))
        SR[e] = (l,u)
    end
    
    res = bboptimize(objective; SearchRange = SR, Method = method,TraceMode = :silent); #MaxSteps=70000,
    initial_parameters = best_candidate(res);
    fitness_1 = best_fitness(res);
    ### Second step ###
    inner_optimizer = NelderMead()
    results = optimize(objective, lower, upper, initial_parameters, Fminbox(inner_optimizer), Optim.Options(time_limit = 50.0)); #20.0
    fitness_2 = results.minimum
    best = results.minimizer
    parameters = fitness_2 < fitness_1 ? best : initial_parameters

    return parameters
end

# Optional: The parameter bounds for this particular application can be made more specific than the bounds that are generally used in the package,
# which is application-agnostic.
function get_parameter_upper_bound(readablecircuit::String)
    elements = foldl(replace,["["=>"","]"=>"","-"=>"",","=>""],init = denumber_circuit(readablecircuit))
    ranges = Dict('R'=>400,'C'=>0.01,'L'=>5,'P'=>[400,1],'W'=>400,'+'=>0,'-'=>0) 
    return flatten([ranges[e] for e in elements])
end

# Optimise the circuit parameters and simulate the resulting impedance spectrum.
function optimise_and_simulate(circuitstring,measurements,frequencies,optim_method)
    parameters = parameopt(circuitstring,measurements,frequencies,optim_method)
    trace = simulateimpedance_noiseless(circuitfunction(circuitstring),parameters,frequencies)
    return trace
end

#Simulate impedance spectra for circuit parameters using the considered optimisation methods.
function generate_traces(circuitstring,gamry_data,measurements,frequencies_,opt_methods,method_names)
    names_ = Vector{String}(undef,length(opt_methods)+2)
    traces_ = Vector{Vector{ComplexF64}}(undef,length(opt_methods)+2)
    traces_[1] = Z_data
    names_[1] = "Measurements"
    traces_[2] = gamry_data
    names_[2] = "Gamry_fit"
    for i in 3:length(traces_)
        traces_[i] = optimise_and_simulate(circuitstring,measurements,frequencies_,opt_methods[i-2])
        names_[i] = method_names[i-2]
    end
return traces_, names_
end

# The fitting errors calculated using the modulus weighted objective function,
# you can adjust the function to see other fitting quality metrics (e.g. removal of the denominator gives the MSE).
function trace_quality(measurements,trace)
    return mean((abs.(measurements - trace).^2)./(abs.(measurements).^2 .+ abs.(trace).^2))
end

# Nyquist plots of arbitraty number of impedance spectra, compared to measurements and Gamry's fit.

function plot_Nyquist(traces_::Vector{Vector{ComplexF64}},names::Vector{String}, _frequ_data::Vector{Float64}) #_Z_Gamry::Vector{ComplexF64}, _Z_GenAlg1::Vector{ComplexF64}, _Z_GenAlg2::Vector{ComplexF64},
    _dtick = 2.0e-3
    s_pts_info = []
    for i_ndx in eachindex(_frequ_data)
        push!(s_pts_info, @sprintf("#:%i, f=%.2fHz", i_ndx, _frequ_data[i_ndx]))
    end
    # --- traces:
    T = Array{GenericTrace{Dict{Symbol, Any}}}(undef,length(traces_))
    for i in 1:length(traces_)
        if i == 1
        T[i] = PlotlyJS.scatter(; x= real(traces_[i]),    y=  imag(traces_[i]),    name = names[i], text = s_pts_info, mode = "markers")
        else
        T[i] = PlotlyJS.scatter(; x= real(traces_[i]),    y=  imag(traces_[i]),    name = names[i])
        end
    end
    # --- layout:
    plt_layout = PlotlyJS.Layout(;
        title_text          = "Comparison of Fit Results",
        xaxis_title_text    = "z<sub>Real</sub> / Ω",
        xaxis_dtick         = _dtick,
        yaxis_title_text    = "z<sub>Imag</sub> / Ω",
        yaxis_dtick         = _dtick,
        yaxis_scaleanchor   = "x",
        yaxis_scaleratio    = 1,
    )
    return PlotlyJS.Plot(T, plt_layout)
end

#Plot the fittin errors to compare the different optimisation methods.
function Plot_fitting_errors(method_names,fit_errors)
    bar_plot = PlotlyJS.bar(x=method_names, y=fit_errors, text=y, textposition="auto")
    plt_layout = PlotlyJS.Layout(;
        title_text          = "Comparison of fiting methods",
        yaxis_title_text    = "Fitting error",
        yaxis_scaleanchor   = "x",
        yaxis_scaleratio    = 1,
    )
    return PlotlyJS.Plot(bar_plot,plt_layout)
end

#Evaluate the optimisation methods and make plots:

@time begin
    traces_,names_ = generate_traces(circ_strg_ref, gamry_trace, Z_data, frequ_data, bbopt_methods, bbopt_method_names)
end
Fitting_errors = [trace_quality(Z_data,traces_[i]) for i in 1:length(traces_)]
errors_plt = Plot_fitting_errors(names_,Fitting_errors)
nyqs_plot = plot_Nyquist(traces_,names_,frequ_data)