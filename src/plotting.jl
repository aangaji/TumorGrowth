export plotting, plotting!, plotting_colored_mutations, plotting_colored_mutations!, colors_by_mutations
export plotting_sampletumor_pies!, plotting_sampletumor_pies

function plotting!(scene, tumor;
        path="", markersize = 1.0, color=nothing, colormap=distinguishable_colors(68), alpha = 1.,
        plotargs...
        )

    isempty(tumor) && return
    p = (tumor isa DataFrameRow) ? [Point(tumor.position...)] : [Point(pos...) for pos in tumor.position]

    colors = if isnothing(color)
        [isempty(mut) ? 1 : mod(mut[end], length(colormap)) + 1 for mut in tumor.mutations]
    elseif !(color isa Vector)
        (color isa Number ? colormap[color] : color, alpha)
    elseif length(color) == nrow(tumor)
        [(c isa Number ? colormap[c] : c, alpha) for c in color]
    else
        error("Options for color: nothing, single value or value for every cell.")
    end
    
    meshscatter!(scene, p; markersize = markersize, color = colors)
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
        sub=1, autodepth = false, show_warning = true, alpha = 1., plotargs...
        )

    isempty(tumor) && return
    points, colors = colors_by_mutations(tumor; colorpalette = colorpalette, sub = sub, autodepth = autodepth, show_warning = show_warning)

    newscene = meshscatter!(scene, points, markersize = markersize, color = [(c,alpha) for c in colors], scale_plot = false, shading=shading, limits = limits, plotargs...)

    isempty(path) || save(path, scene)
    return scene
end

function plotting_colored_mutations(tumor;
        size=(500,500), show_axis=false, inline=false, plotargs...)
    Makie.inline!(inline)
    scene = Scene(resolution=size, show_axis=show_axis)
    plotting_colored_mutations!(scene, tumor; plotargs...)
end

function plotting_sampletumor_pies!(scene, sampletumor;
        path="", colorpalette = palette(:tab20), overdraw=false,
        strokecolor=nothing, strokewidth = 1, sample_r = 0.,
        mutations = nothing, first = isnothing(mutations) ? 10 : length(mutations)
        )

        top_muts = isnothing(mutations) ? sort( sampletumor_mfreqs(sampletumor), :frequency, rev=true).mutation[1:first,:] : mutations
        m_colors = Dict( top_muts .=> colorpalette[1:first])

        for s in eachrow(sampletumor)
                mask = map( m-> m in top_muts, s.mutations)
                center, freqs, colors = s.position, s.frequencies[mask], [m_colors[m] for m in s.mutations[mask]]
                if isempty(freqs)
                        freqs, colors = [1.],[:lightgrey]
                end
                freqs ./= sum(freqs)
                startphis = cumsum(freqs).-freqs
                sections = map(zip(startphis,freqs)) do (fs,f)
                        r = iszero(sample_r) ? s.sample_r : sample_r
                        circ = Point2f0[ center .+ r .* (cos(2π*p), sin(2π*p)) for p in (0.:0.002:f).+fs]
                        !isone(length(freqs)) && push!(circ,Point2f0(center))
                        circ
                end
                for (sec,c) in zip(sections,colors)
                        poly!(scene, sec , color = c, strokewidth = strokewidth,
                        strokecolor = isnothing(strokecolor) ? c : strokecolor
                         )
                end
        end
        isempty(path) || save(path, scene)
        return scene
end

function plotting_sampletumor_pies(sampletumor;
        size=(500,500), show_axis=false, inline=false, plotargs...)
        Makie.inline!(inline)
        scene = Scene(resolution=size, show_axis=show_axis)
        plotting_sampletumor_pies!(scene, sampletumor; plotargs...)
end
