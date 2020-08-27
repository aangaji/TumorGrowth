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
    isempty(path) || save(path, scene)
end
function plotting(data::DataFrame; path="", color=nothing, shading=false, inline=false)
    AbstractPlotting.inline!(inline)
    scene = Scene()
    plotting!(scene, data::DataFrame; path=path, color=color, shading=shading)
    return scene
end


#function to plot the tumor similar to Ling et al. (2015), color code in experiments

function plotting_colored_mutations!(scene, tumor::DataFrame; colorpalette = Plots.palette(:tab10), path="", shading = false, limits = AbstractPlotting.automatic, sub=1)
    isempty(tumor) && return

    p = [Point(pos...) for pos in tumor.position]

    subclones = clones(tumor; sub=sub)
    for (i,c) in enumerate(subclones)
         meshscatter!(scene, [Point(pos...) for pos in c.position], markersize = 1., color = colorpalette[(i-1)%length(colorpalette)+1], scale_plot = false, shading=shading, limits = limits)
    end
    isempty(path) || save(path, scene)
end

function plotting_colored_mutations(tumor::DataFrame; colorpalette=Plots.palette(:tab10), path="", shading = false, limits = AbstractPlotting.automatic, inline=false)
    AbstractPlotting.inline!(inline)
    scene = Scene()
    plotting_colored_mutations!(scene, tumor; colorpalette = colorpalette, path=path, shading = shading, limits = limits)
    return scene
end
