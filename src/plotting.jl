export plotting, plotting!, plotting_colored_mutations, plotting_colored_mutations!, colors_by_mutations

function plotting!(scene, tumor;
        path="", limits = Makie.automatic,
        markersize = 1.0, color=nothing, colormap=distinguishable_colors(68), shading = false
        )

    isempty(tumor) && return
    p = (tumor isa DataFrameRow) ? [Point(tumor.position...)] : [Point(pos...) for pos in tumor.position]

    isnothing(color) && (color = [isempty(mut) ? 1 : mod(mut[end], length(colormap))+1 for mut in tumor.mutations])
    meshscatter!(scene, p , markersize = markersize, color = color, colormap=colormap, scale_plot = false, shading=shading, limits = limits)
    isempty(path) || save(path, scene)
    return scene
end
function plotting(tumor;
        size=(500,500), show_axis=false, inline=false, plotargs...)

    Makie.inline!(inline)
    scene = Scene(resolution=size, show_axis=show_axis)
    plotting!(scene, tumor; plotargs...)
end

function colors_by_mutations(tumor; colorpalette = palette(:tab20), sub = 1, autodepth = false, show_warning = true)
    points = [Point(pos...) for pos in tumor.position]
    colors = fill(colorpalette[1], length(points))
    subclones = clones(tumor; sub=sub, autodepth=autodepth, show_warning = show_warning)
    for (i,cl) in enumerate(subclones)
        colors[findall(in(cl.index), tumor.index)] .= colorpalette[(i-1)%length(colorpalette)+1]
    end
    return points, colors
end

function plotting_colored_mutations!(scene, tumor;
        path="", limits =  Makie.automatic,
        markersize = 1.0, colorpalette = palette(:tab20), shading = false,
        sub=1, autodepth = false, show_warning = true
        )

    isempty(tumor) && return
    points, colors = colors_by_mutations(tumor; colorpalette = colorpalette, sub = sub, autodepth = autodepth, show_warning = show_warning)

    newscene = meshscatter!(scene, points, markersize = markersize, color = colors, scale_plot = false, shading=shading, limits = limits)

    isempty(path) || save(path, scene)
    return scene
end

function plotting_colored_mutations(tumor;
        size=(500,500), show_axis=false, inline=false, plotargs...)
    Makie.inline!(inline)
    scene = Scene(resolution=size, show_axis=show_axis)
    plotting_colored_mutations!(scene, tumor; plotargs...)
end
