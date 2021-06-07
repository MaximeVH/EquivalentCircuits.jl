# EquivalentCircuits.jl

[![Build Status](https://travis-ci.com/MaximeVH/EquivalentCircuits.jl.svg?branch=master)](https://travis-ci.com/MaximeVH/EquivalentCircuits.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/MaximeVH/EquivalentCircuits.jl?svg=true)](https://ci.appveyor.com/project/MaximeVH/EquivalentCircuits-jl)
[![Coverage](https://coveralls.io/repos/github/MaximeVH/EquivalentCircuits.jl/badge.svg?branch=master)](https://coveralls.io/github/MaximeVH/EquivalentCircuits.jl?branch=master)


# EquivalentCircuits.jl

This Julia module allows users to analyse their **electrochemical impedance spectroscopy** data using **equivalent electrical circuits**. EquivalentCircuits.jl can be used to either find the optimal paramters of a given equivalent electrical circuit , or to get recommendations for an appropriate equivalent electrical circuit configuration. The latter is based on a **gene expression programming** approach.

## Installation
The package can be installed using the package manager.
```julialang
] add EquivalentCircuits
```

## Usage
### Circuit notation
Equivalent electrical circuit models are composed of electrical elements, connected in series or in parallel. The four fundamental elements that are most commonly encountered in equivalent electrical circuits, are resistors, capacitors, inductors and constant phase elements. These four elements are represented by the capital letters R, C, L and P, respectively. serially connected elements have dashes between them, wereas parallely connected elements are placed between square brackets and separated by a comma. Finally all the elements in a circuit are numbered. Using these notation rules, the circuit `R1-[C2,R3-[C4,R5]]` corresponds to: ![](example_circuit.png)

### Parameter fitting
When an appropriate circuit model is available, the parameters can be fitted to experimental data using the `parameteroptimisation(circuit,data)` function. the `circuit`argument is the equivalent circuit, provided as a string with the circuit notation displayed above. The electrochemical impedance measurement data can be provided as a CSV file with three columns: imaginary impedance, real impedance and frequency (see example_measurements.csv).

### Circuit fitting
When only the electochemical impedance measurements are available, equivalent electrical circuit recommendations can be obtained using the `circuitevolution(data;kwargs)` function. The data can once again be provided as CSV file (or alternatively as a DataFrame). A variety of keyword arguments can be adjusted to adjust the gene expression programming circuit identification procedure.
