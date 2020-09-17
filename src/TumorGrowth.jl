module TumorGrowth

begin
print("Loading Packages... ")
import DataFrames: DataFrame, DataFrame!, DataFrameRow
using CSV
using Plots: palette, distinguishable_colors
using Makie
using LinearAlgebra
using StaticArrays
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
    data = CSV.File(path, delim="\t") |> DataFrame!
    fields = names(data)
    types = typeof.([first(data)[field] for field in fields])

    for field in fields[isequal.(types, String)]
        data[!,field] = data[!,field] .|> Meta.parse .|> eval
    end
    return data
end

end  # module TumorGrowth
