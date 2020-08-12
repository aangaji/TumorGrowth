export plotting_2d, plotting_2d!, plotting_2d_colored_mutations, plotting, plotting!, plotting_colored_mutations

#function to create circles in a 2D plot (here: cells with radius r)
function circle(pos_x, pos_y)
    theta = LinRange(0, 2*pi, 20)
    return pos_x .+ sin.(theta), pos_y .+ cos.(theta)
end

#function to plot a tumor with each mutation colored differently
function plotting_2d!(p, data::DataFrame; path="", annotate=false, color=nothing)

    pos = data.position
    color_scheme = palette(:tab10)
    choose_color = false
    isnothing(color) && (choose_color=true)
    for (i,mut) in enumerate(data.mutations)
        color_choice = 1
        if ! isempty(mut)
            color_choice = mod(maximum(mut), 10)+1 #choose a different color for each mutation out of a colorset with 10 colors
            annotate && annotate!(pos[i]..., text(mut[1], 5, :black)) #annotate the number of the first mutation to each cell
        end
        choose_color && (color = color_scheme[color_choice])
        plot!(p, circle(pos[i]...), seriestype = [:shape,], lw = 0.5, c = color, linecolor = :black, lab = :none, fillalpha = 0.5, aspect_ratio = 1) #plot the cell with its specific color choice
    end

    first = findfirst(isequal(1), data.id)
    isnothing(first) || plot!(p, circle(pos[first]...), seriestype = [:shape,], lw = 0.5, c = :black, linecolor = :black, lab = "first", fillalpha = 1.0, aspect_ratio = 1) #plot the start cell (lowest id) with black filling
    isempty(path) || savefig("$(path)")
    return p
end

plotting_2d(data::DataFrame; path="", annotate=false) = plotting_2d!(plot(), data; path=path, annotate=annotate)

#function to plot the tumor similar to Ling et al. (2015), color code in experiments
function plotting_2d_colored_mutations(data::DataFrame; path="")

    color_scheme = palette(:tab10)
    p = plot()

    for (i,color) in enumerate(colors_by_muts(data.mutations; color_set_size=10))
        plot!(p, circle(data.position[i]...), seriestype = [:shape,], lw = 0.5, c = color_scheme[color], linecolor = :black, lab=:none, fillalpha = 1.0, aspect_ratio = 1)
    end

    first = findfirst(isequal(1), data.id)
    isnothing(first) || plot!(p, circle(data.position[first]...), seriestype = [:shape,], lw = 0.5, c = :black, linecolor = :black, lab = "first", fillalpha = 1.0, aspect_ratio = 1)

    isempty(path) || savefig("$(path)_colored", p)
    return p
end

#function to plot a tumor with color determined by the latest mutation of the cell
function plotting!(scene, data::DataFrame; path="", color=nothing, dim = length(data.position[1]))

    p = (getindex.(data.position, i) for i=1:dim)
    isnothing(color) && (color = [isempty(mut) ? 1 : mod(mut[end], 68)+1 for mut in data.mutations])
    meshscatter!(scene, p..., markersize = 1.0, color = color, colormap=distinguishable_colors(68), scale_plot = false)
    isempty(path) || save("$(path)", scene)
    return scene
end
plotting(data::DataFrame; path="", color=nothing, dim=length(data.position[1])) = plotting!(Scene(), data::DataFrame; path=path, color=color, dim=dim)


#function to plot the tumor similar to Ling et al. (2015), color code in experiments
function plotting_colored_mutations(data::DataFrame; path="", dim=length(data.position[1]))

    p = ( getindex.(data.position, i) for i=1:dim )

    color_choice = colors_by_muts(data.mutations; color_set_size=10)

    scene = meshscatter(p..., markersize = 1., color = color_choice, colormap=:tab10, scale_plot = false)

    isempty(path) || save("$(path)_colored", scene)
    return scene
end

function colors_by_muts(mutations; color_set_size=1) :: Vector{Int64}

    color_choice = fill(color_set_size, length(mutations))

    colored_mutations = Vector{Vector{Int64}}(undef,0)

    for (i, mut_list) in enumerate(mutations)
        color_found = false
        if ! isempty(mut_list) #if the cell has mutations
            if isempty(colored_mutations) #if this is the first mutated cell
                push!(colored_mutations, mut_list) #push the mutations in colored mutations
                color_choice[i] = 1
            else
                for (j, cmut) in enumerate(colored_mutations) #search for the specific color (color j) and plot the cell with this color (multiple times if the succesive mutations are already known)
                    if mut_list[1] in cmut
                        for mut in mut_list
                            color_choice[i] = mod(j,color_set_size)+1
                            if !(mut in cmut)
                                push!(colored_mutations[j], mut)
                            end
                        end
                        color_found = true
                    end
                end
                if !color_found #if the first mutation of a cell wasn't found already, choose a new color by pushing it into the array colored_mutations
                    color_choice[i] = mod(length(colored_mutations)+1,color_set_size)+1
                    push!(colored_mutations, mut_list)
                end
            end
        end
    end
    return color_choice
end
