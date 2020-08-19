module TumorGrowth

begin
print("Loading Packages... ")
import DataFrames: DataFrame!, DataFrame
using CSV
using Plots
using Makie: Scene, meshscatter, meshscatter!, text!, save
using LinearAlgebra
using Random
#using Juno
using ProgressMeter
println("Done!")
end

(
    "pushing_simulation.jl",
    "plotting.jl",
    "sampling.jl",
    "pushing_prototype.jl",
    "analysis.jl"
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
