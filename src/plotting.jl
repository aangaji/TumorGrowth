export plotting, plotting!, plotting_colored_mutations, plotting_colored_mutations!

#function to plot a tumor with color determined by the latest mutation of the cell
function plotting!(scene, data::DataFrame; path="", color=nothing, shading = false)
    isempty(data) && return
    p = [Point(pos...) for pos in data.position]
    isnothing(color) && (color = [isempty(mut) ? 1 : mod(mut[end], 68)+1 for mut in data.mutations])
    meshscatter!(scene, p, markersize = 1.0, color = color, colormap=distinguishable_colors(68), scale_plot = false, shading=shading)
    isempty(path) || save(path, scene)
end
function plotting(data::DataFrame; path="", color=nothing, shading=false, inline=false)
    AbstractPlotting.inline!(inline)
    scene = Scene()
    plotting!(scene, data::DataFrame; path=path, color=color, shading=shading)
    return scene
end


#function to plot the tumor similar to Ling et al. (2015), color code in experiments

function plotting_colored_mutations!(scene, tumor::DataFrame; colorpalette = palette(:tab20), path="", shading = false, limits = automatic, sub=1)
    isempty(tumor) && return

    p = [Point(pos...) for pos in tumor.position]

    subclones = clones(tumor; sub=sub)
    for (i,c) in enumerate(subclones)
         meshscatter!(scene, [Point(pos...) for pos in c.position], markersize = 1., color = colorpalette[(i-1)%length(colorpalette)+1], scale_plot = false, shading=shading, limits = limits)
    end
    isempty(path) || save(path, scene)
end

function plotting_colored_mutations(tumor::DataFrame; colorpalette=palette(:tab10), path="", shading = false, limits = automatic, inline=false)
    AbstractPlotting.inline!(inline)
    scene = Scene()
    plotting_colored_mutations!(scene, tumor; colorpalette = colorpalette, path=path, shading = shading, limits = limits)
    return scene
end
