module TumorGrowth

begin
print("Loading Packages... ")
import DataFrames: DataFrame!, DataFrame
using CSV
using Plots: palette, distinguishable_colors
using Makie
using LinearAlgebra
using StatsBase: fit, Histogram
using Random
using AbstractPlotting: automatic
using ProgressMeter
println("Done!")
end

(
    "pushing_simulation.jl",
    "plotting.jl",
    "pushing.jl",
    "pushing_animate.jl",
    "sampling.jl",
    "analysis.jl",
    "time_series.jl"
) .|> include


export data_import

#function to import a saved tumor from a .csv file
function data_import(path::String)
    data = path |> CSV.File |> DataFrame!
    for field in data |> names .|> Symbol
        try data[!, field] = data[!, field] .|> Meta.parse .|> eval catch e end
    end
    return data
end

end  # module TumorGrowth
