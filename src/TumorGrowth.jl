module TumorGrowth

print("Loading Packages... ")
import DataFrames: DataFrame
using DataFrames: DataFrameRow, nrow
using CSV
using Plots: palette, distinguishable_colors
using GLMakie
using LinearAlgebra
using Statistics
using StaticArrays
using StatsBase: fit, Histogram
using Random
using Distributions
using ProgressMeter
println("Done!")

(
    "simulation.jl",
    "plotting.jl",
    "pushing.jl",
    "pushing_animate.jl",
    "sampling.jl",
    "analysis.jl",
    "time_series.jl",
    "nonspatial.jl"
) .|> include

export data_import

"""
    data_import(path::String; delim="\t") :: DataFrame

Convenience function which loads a DataFrame from a '.csv' file and correctly
parses columns with elementtype 'Vector'.
"""
function data_import(path::String; delim="\t")
    data = DataFrame(CSV.File(path; delim=delim); copycols=false)
    arraycols = collect(first(data)) .|> val -> isequal(typeof(val), String) && all(occursin.(("[","]"), val))
    data[!,arraycols] .= data[!,arraycols] .|> Meta.parse .|> eval
    for col in findall(.!(arraycols))
            data[!, col] .= data[!, col] |> Vector
    end
    return data
end

end  # module TumorGrowth
