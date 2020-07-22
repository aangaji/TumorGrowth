module TumorGrowth

begin
print("Loading Packages... ")
using DataFrames
using CSV
using Plots
using Makie: meshscatter, text!, save
using LinearAlgebra
using Random
println("Done!")
end




(
    "2d_pushing.jl",
    "2d_restricted.jl",
    "3d_pushing.jl",
    "3d_restricted.jl",
    "plotting_functions.jl",
    "sampling_functions.jl",
    "pushing_prototype.jl"
) .|> include


const ε = 300            ### THESE SHOULD PROBABLY BE CONSTANTS
const r0 = 1
const kb = 1
const T = 310.15 #37 degree Celsius

#Lenneard-Jones like potential to describe the repulsive force by pushing away when there is an overlap
U(r1, r2) = 4 * ε * ((2 * r0 / norm(r1 - r2))^(12) - (2 * r0 / norm(r1 - r2))^6)



#####################################################
################### DATA Handling ###################
#####################################################



#function to import a saved tumor correctly (especially Arrays aren't imported correctly (here positions and mutations))
#works for 2D AND 3D
function data_import(path::String)
    #Import data and convert positions and mutations back to array type
    data = CSV.read(path)
    data[!, :position] = Meta.parse.(data[!, :position])
    B = Array{Array}(undef, length(data[!, :index]))
    for i = 1:length(data[!, :index])
        B[i] = data[!, :position][i].args
    end
    data[2] = B

    data[!, :mutations] = Meta.parse.(data[!, :mutations])
    C = Array{Array}(undef, length(data[!, :index]))
    for i = 1:length(data[!, :index])
        C[i] = deleteat!(data[!, :mutations][i].args, 1)
    end

    data[!, :mutations] = C

    return data
end

function data_import_transpose(path::String)
    #Import data and convert positions and mutations back to array type
    data = CSV.read(path)
    data = Meta.parse.(data)

    B = Array{Array}(undef, length(data[2, :]))         ### TYPE OF ARRAY
    for i = 1:length(data[2, :])
        B[i] = data[2, :][i].args
    end
    data[2, :] = B

    C = Array{Array}(undef, length(data[2, :]))
    for i = 1:length(data[2, :])
        C[i] = deleteat!(data[6, :][i].args, 1)
    end

    data[6, :] = C

    return data
end


end  # module TumorGrowth
