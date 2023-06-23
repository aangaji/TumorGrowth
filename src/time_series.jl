export tumor_stepper, tumor_stepper!, record_growth

function tumor_stepper!(simresult, trange;
        seed = nothing)
    isnothing(seed) || Random.seed!(seed)
    id, mut_id, t, tumor, mutations = simresult |> r -> (r[:index], r[:mutation], r[:time], r[:tumor], r[:mutations])

    tumor_dict = Dict(getfield.(tumor,:index) .=> tumor)

    prog = Progress(length(trange), dt=0.01)
    time_series = map(trange) do tp
        up = birth_death_pushing!(tumor_dict, mutations, tp; cur_id=id, cur_mutation=mut_id, t=t, showprogress=false)
        id, mut_id, t = up[:index], up[:mutation], up[:time]

        next!(prog)

        DataFrame(values(tumor_dict))
    end
    return time_series
end
function tumor_stepper(trange;
        seed = abs(rand(Int)), dim, simparams...)
    simresult = birth_death_pushing(first(trange); dim = dim, simparams..., seed=seed)
    tumor_stepper!(simresult, trange)
end

function record_growth(time_series; path="temp.gif", frames=1,
        shading = false, markersize=1., size = (500,500),
        show_axis=false, colorpalette = palette(:tab20),
        points_colors = snap -> colors_by_mutations(snap; colorpalette=colorpalette, show_warning=false)
        )
    dim = length(first(first(time_series).position))
    width = maximum(abs.(vcat(last(time_series).position...)))
    origin = fill(width, dim)
    limits = dim == 2 ? FRect2D(-origin, 2*origin) : FRect3D(-origin, 2*origin)

    pointsNode = Node([Point(zeros(dim)...)])
    colorsNode = Node([colorpalette[1]])
    scene = Scene(resolution=size, show_axis=show_axis)
    meshscatter!(scene, pointsNode, color = colorsNode, scale_plot=false, limits = limits, shading = shading, markersize = markersize)
    record(scene, path) do io
        @showprogress for snapshot in time_series
            pointsNode[], colorsNode[] = points_colors(snapshot)
            for _=1:frames recordframe!(io) end
        end
    end
end
