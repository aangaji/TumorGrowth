module TumorGrowth

print("Loading Packages... ")
import DataFrames: DataFrame, DataFrame!, DataFrameRow, nrow
using CSV
using Plots: palette, distinguishable_colors
using Makie
using LinearAlgebra
using Statistics
using StaticArrays
using StatsBase: fit, Histogram
using Random
using Distributions
using AbstractPlotting: automatic
using ProgressMeter
println("Done!")

(
    "simulation.jl",
    "plotting.jl",
    "pushing.jl",
    "pushing_animate.jl",
    "sampling.jl",
    "analysis.jl",
    "time_series.jl"
) .|> include

_version = "1.0_sim-rework"

export data_import

# function to import a saved tumor from a .csv file
function data_import(path::String; delim="\t")
    data = CSV.File(path, delim=delim) |> DataFrame!
    fields = names(data)
    types = typeof.([first(data)[field] for field in fields])

    for field in fields[isequal.(types, String)]
        data[!,field] = data[!,field] .|> Meta.parse .|> eval
    end
    return data
end

end  # module TumorGrowth
