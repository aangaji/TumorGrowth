export cross_section_3D_simulation, cross_section_3D_simulation_transpose

#determine cells in the cross sections along the x-axis and the y-axis
#(attention: take all cells which lay partly in this plane -> distance from center of the cell to the plane is less than the radius)
function cross_section_3D_simulation(data::DataFrame)

    cross_section_x = DataFrame()
    cross_section_y = DataFrame()
    cross_section_z = DataFrame()

    for i in 1:3
        for j in 1:length(data.position)
            if abs(data.position[j][i]) <= r0
                if i == 1
                    push!(cross_section_x, data[j,:])
                end
                if i == 2
                    push!(cross_section_y, data[j,:])
                end
                if i == 3
                    push!(cross_section_z, data[j,:])
                end
            end
        end
    end

    # calculate cross sections as real 2D to use the already implemented functions (no special functions for 3D)

    for i in eachrow(cross_section_x)
        i.position = i.position[2:3]
    end

    for i in eachrow(cross_section_y)
        i.position = i.position[1:2:3]
    end

    for i in eachrow(cross_section_z)
        i.position = i.position[1:2]
    end

    return cross_section_x, cross_section_y, cross_section_z
end
