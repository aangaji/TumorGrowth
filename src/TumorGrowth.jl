module TumorGrowth

begin
print("Loading Packages... ")
using DataFrames
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
    "2d_pushing.jl",
    "2d_restricted.jl",
    "3d_pushing.jl",
    "3d_restricted.jl",
    "plotting.jl",
    "sampling.jl",
    "pushing_prototype.jl",
    "analysis.jl"
) .|> include


const ε = 300
const r0 = 1
const kb = 1
const T = 310.15 #37 degree Celsius

#Lenneard-Jones like potential to describe the repulsive force by pushing away when there is an overlap
U(r1, r2) = 4 * ε * ((2 * r0 / norm(r1 - r2))^(12) - (2 * r0 / norm(r1 - r2))^6)


export data_import

#function to import a saved tumor from a .csv file
function data_import(path::String)
    data = DataFrame!(CSV.File(path))
    data.position = data.position .|> Meta.parse .|> eval
    data.mutations = data.mutations .|> Meta.parse .|> eval
    return data
end

end  # module TumorGrowth
