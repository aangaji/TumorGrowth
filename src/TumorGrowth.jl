module TumorGrowth

begin
print("Loading Packages... ")
import DataFrames: DataFrame, DataFrameRow
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
    data = DataFrame!(CSV.File(path))
    data.position = data.position .|> Meta.parse .|> eval
    data.p_birth = data.p_birth .|> Meta.parse .|> eval
    data.mutations = data.mutations .|> Meta.parse .|> eval
    return data
end

end  # module TumorGrowth
