export plotting, plotting!, plotting_colored_mutations, plotting_colored_mutations!

#function to create circles in a 2D plot (here: cells with radius r)
function circle(pos_x, pos_y)
    theta = LinRange(0, 2*pi, 20)
    return pos_x .+ sin.(theta), pos_y .+ cos.(theta)
end

#function to plot a tumor with each mutation colored differently
function plotting_2d!(p, data::DataFrame; path="", annotate=false, color=nothing)
    isempty(data) && return
    pos = data.position
    color_scheme = palette(:tab10)
    choose_color = false
    isnothing(color) && (choose_color=true)
    for (i,mut) in enumerate(data.mutations)
        color_choice = 1
        if ! isempty(mut)
            color_choice = mod(maximum(mut), 10)+1 #choose a different color for each mutation out of a colorset with 10 colors
            annotate && annotate!(pos[i]..., Plots.text(mut[1], 5, :black)) #annotate the number of the first mutation to each cell
        end
        choose_color && (color = color_scheme[color_choice])
        Plots.plot!(p, circle(pos[i]...), seriestype = [:shape,], lw = 0.5, c = color, linecolor = :black, lab = :none, fillalpha = 0.5, aspect_ratio = 1) #plot the cell with its specific color choice
    end

    first = findfirst(isequal(1), data.index)
    isnothing(first) || Plots.plot!(p, circle(pos[first]...), seriestype = [:shape,], lw = 0.5, c = :black, linecolor = :black, lab = "first", fillalpha = 1.0, aspect_ratio = 1) #plot the start cell (lowest index) with black filling
    isempty(path) || Plots.savefig("$(path)")
    return p
end

plotting_2d(data::DataFrame; path="", annotate=false) = plotting_2d!(Plots.plot(), data; path=path, annotate=annotate)

#function to plot the tumor similar to Ling et al. (2015), color code in experiments
function plotting_2d_colored_mutations(data::DataFrame; path="")
    isempty(data) && return
    color_scheme = palette(:tab10)
    p = Plots.plot()

    for (i,color) in enumerate(colors_by_muts(data.mutations; color_set_size=10))
        plot!(p, circle(data.position[i]...), seriestype = [:shape,], lw = 0.5, c = color_scheme[color], linecolor = :black, lab=:none, fillalpha = 1.0, aspect_ratio = 1)
    end

    first = findfirst(isequal(1), data.index)
    isnothing(first) || Plots.plot!(p, circle(data.position[first]...), seriestype = [:shape,], lw = 0.5, c = :black, linecolor = :black, lab = "first", fillalpha = 1.0, aspect_ratio = 1)

    isempty(path) || Plots.savefig("$(path)_colored", p)
    return p
end

#function to plot a tumor with color determined by the latest mutation of the cell
function plotting!(scene, data::DataFrame; path="", color=nothing, shading = false)
    isempty(data) && return
    p = [Point(pos...) for pos in data.position]
    isnothing(color) && (color = [isempty(mut) ? 1 : mod(mut[end], 68)+1 for mut in data.mutations])
    meshscatter!(scene, p, markersize = 1.0, color = color, colormap=Plots.distinguishable_colors(68), scale_plot = false, shading=shading)
    isempty(path) || save("$(path)", scene)
end
function plotting(data::DataFrame; path="", color=nothing, shading=false, inline=false)
    AbstractPlotting.inline!(inline)
    scene = Scene()
    plotting!(scene, data::DataFrame; path=path, color=color, shading=shading)
    return scene
end


#function to plot the tumor similar to Ling et al. (2015), color code in experiments

function plotting_colored_mutations!(scene, data::DataFrame; colormap=:tab10, color_set_size = 10, path="", shading = false, limits = AbstractPlotting.automatic)
    isempty(data) && return
    p = [Point(pos...) for pos in data.position]

    color_choice = colors_by_muts(data.mutations; color_set_size=color_set_size)

    meshscatter!(scene, p, markersize = 1., color = color_choice, colormap=colormap, scale_plot = false, shading=shading, limits = limits)

    isempty(path) || save("$(path)_colored", scene)
end

function plotting_colored_mutations(data::DataFrame; colormap=:tab10, color_set_size = 10, path="", shading = false, limits = AbstractPlotting.automatic, inline=false)
    AbstractPlotting.inline!(inline)
    scene = Scene()
    plotting_colored_mutations!(scene, data; colormap=colormap, color_set_size = color_set_size, path=path, shading = shading, limits = limits)
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
