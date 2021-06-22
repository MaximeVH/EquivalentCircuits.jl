# EquivalentCircuits.jl

This Julia module allows users to analyse their **electrochemical impedance spectroscopy** data using **equivalent electrical circuits**. EquivalentCircuits.jl can be used to either fit the parameters of a given equivalent electrical circuit , or to get recommendations for an appropriate equivalent electrical circuit configuration. The latter is based on a [gene expression programming](https://en.wikipedia.org/wiki/Gene_expression_programming) approach.

## Installation
The package can be installed using the package manager.
```julialang
] add EquivalentCircuits
```

## Usage
### Circuit notation
Equivalent electrical circuit models are composed of electrical elements, connected in series or in parallel. The four fundamental elements that are most commonly encountered in equivalent electrical circuits, are resistors, capacitors, inductors and constant phase elements. These four elements are represented by the capital letters R, C, L and P, respectively. serially connected elements have dashes between them, wereas parallely connected elements are placed between square brackets and separated by a comma. Finally all the elements in a circuit are numbered. Using these notation rules, the circuit `R1-[C2,R3-[C4,R5]]` corresponds to:

![](example_circuit.png)

When using this package, the circuit should input as a String:  
```julia
using EquivalentCircuits

circuit = "R1-[C2,R3-[C4,R5]]"
```

### Parameter fitting
When an appropriate circuit model is available, the parameters can be fitted to experimental data using the `parameteroptimisation` function which accepts two arguments:
- `circuit` : the equivalent circuit, provided as a string with the circuit notation displayed above.
- `data` : the filepath of the electrochemical impedance measurement data.
The data should be provided as a CSV file with three columns: imaginary impedance, real impedance and frequency (see example_measurements.csv).

Lets first take a look at what the çontents of the [example_measurements.csv](https://github.com/MaximeVH/EquivalentCircuits.jl/blob/master/example_measurements.csv) file look like:

```julia
using CSV, DataFrames

#Load the measurement data.

data = "example_measurements.csv"; #This should be the filepath of the example_measurements.csv file.

df = CSV.read("example_measurements.csv",DataFrame,header = false);

#Rename the columns for illustration purposes.

rename_dict = Dict("Column1"=>"Reals","Column2"=>"Imags","Column3"=>"Frequencies");

rename!(df, rename_dict);

println(df)

```

Next we can fit the parameters of our example ciruit to the example measurement data as follows:
```julia
circuitparams = parameteroptimisation(circuit,data)
```

### Circuit fitting
When only the electochemical impedance measurements are available, equivalent electrical circuit recommendations can be obtained using the `circuitevolution(data;kwargs)` function. The data can once again be provided as CSV file. A variety of keyword arguments can be adjusted to fine-tune the gene expression programming circuit identification procedure.The possible keyword agruments to tune the cirucit identification are:

- `generations` : the maximum number of algorithm iterations.
- `population_size` : the number of individuals in the population during each iteration.
- `terminals` : the circuit components that are to be included in the circuit identification.
- `cutoff`: a hyperparameter that controls the circuit complexity by removing redundant components. Lower values lead to more simple circuits, however too low values will lead to circuits that no longer fit the measurements well.
- `head` : a hyperparameter than controls the maximum considered complexity of the circuits.
- `initial_population` : the option to provide an initial population of circuits with which the algorithm starts.

The defaults values are as follows:


| Argument      | Default value |
| ----------- | ----------- |
| generations      | 10       |
| population_size   | 30        |
| terminals   | "RCLP"        |
| head   | 8        |
| cutoff  | 0.80     |
| initial_population  | nothing      |

As an example, by running the code below you can see if a circuit can be found, consisting of only resistors and capacitors, that is capable of fitting the example measurement data.

```julia
circuitevolution("example_measurements.csv",terminals="RC")
```

Next, the file [Circuitlibrary.csv](https://github.com/MaximeVH/EquivalentCircuits.jl/blob/master/Circuitlibrary.csv) contains a collection of various circuit topologies. We can allow the algorithm to start from this circuit collection as initial population as follows:

```julia
# Load the population from the CSV file, using the loadpopulation function.
# The input of the loadpopulation should be the filepath of Circuitpopulation.csv.
circuit_library = loadpopulation("Circuitlibrary.csv"); #The input should be the filepath of the Circuitlibrary.csv file.

# Now find a circuit that fits the data, starting from the initial population of circuits
circuitevolution("example_measurements.csv",initial_population = circuit_library)

```
