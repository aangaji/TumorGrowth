export cross_section, radial_sample

function cross_section(tumor; x=nothing, y=nothing, z=nothing, width = 2., reduce_dim=false)
    for (d,coord) in enumerate((x,y,z))
        if !isnothing(coord)
            rows = -width .< getindex.(tumor.position,d).-coord .< width
            if reduce_dim
                plane = tumor[rows, :]
                plane.position = plane.position .|> p->reduce_dimension(p, d)
                return plane
            else
                return tumor[rows, :]
            end
        end
    end
end

# this function returns a copy without entries at dimension d
function reduce_dimension(position, d)
    dim = position |> length
    select = mod.((d-1:2:d+1).-1, dim).+1 |> unique
    return getindex(position, select)
end

radial_sample(tumor; r=2., width=2.) = tumor[-width .< [sqrt(sum(p.^2)) for p in tumor.position] .- r .< width, :]
