export plotting, plotting!, plotting_colored_mutations, plotting_colored_mutations!

function plotting!(scene, tumor;
        path="", limits = automatic,
        markersize = 1.0, color=nothing, colormap=distinguishable_colors(68), shading = false
        )

    isempty(tumor) && return
    p = (tumor isa DataFrameRow) ? [Point(tumor.position...)] : [Point(pos...) for pos in tumor.position]

    isnothing(color) && (color = [isempty(mut) ? 1 : mod(mut[end], 68)+1 for mut in tumor.mutations])
    meshscatter!(scene, p, markersize = markersize, color = color, colormap=colormap, scale_plot = false, shading=shading, limits = limits)
    isempty(path) || save(path, scene)
end
function plotting(tumor;
        size=(500,500), inline=false, path="", limits = automatic,
        markersize = 1.0, color=nothing, colormap=distinguishable_colors(68), shading=false
        )

    AbstractPlotting.inline!(inline)
    scene = Scene()
    resize!(scene, size)
    plotting!(scene, tumor; path=path, markersize = markersize, color=color, colormap=colormap, shading=shading, limits = limits)
    return scene
end

function plotting_colored_mutations!(scene, tumor;
        path="", limits = automatic,
        markersize = 1.0, colorpalette = palette(:tab20), shading = false, sub=1
        )

    isempty(tumor) && return
    for (i,c) in enumerate( clones(tumor; sub=sub) )
         meshscatter!(scene, [Point(pos...) for pos in c.position], markersize = markersize, color = colorpalette[(i-1)%length(colorpalette)+1], scale_plot = false, shading=shading, limits = limits)
    end
    isempty(path) || save(path, scene)
end

function plotting_colored_mutations(tumor;
        size=(500,500), limits = automatic, inline=false,
        markersize = 1.0, colorpalette=palette(:tab10), path="", shading = false, sub=1
        )
    AbstractPlotting.inline!(inline)
    scene = Scene()
    resize!(scene, size)
    plotting_colored_mutations!(scene, tumor; markersize = markersize, colorpalette = colorpalette, path=path, shading = shading, limits = limits, sub=sub)
    return scene
end
