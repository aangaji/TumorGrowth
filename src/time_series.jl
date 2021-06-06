export tumor_stepper, tumor_stepper!, record_growth

function tumor_stepper!(simresult, trange;
        seed = abs(rand(Int)),simparams...)
    Random.seed!(seed)
    id, mut_id, t, tumor = simresult |> r -> (r[:index], r[:mutation], r[:time], r[:tumor])
    # dim = length(first(tumor).position)
    time_series = Vector{DataFrame}()
    prog = Progress(length(trange), dt=0.01)
    for tp in trange
        up = birth_death_pushing!(tumor, tp; simparams..., cur_id=id, cur_mutation=mut_id, t=t, showprogress=false)
        push!(time_series, deepcopy(DataFrame(tumor)) )
        id, mut_id, t = up[:index], up[:mutation], up[:time]

        next!(prog)
    end
    return time_series
end
function tumor_stepper(trange;
        seed = abs(rand(Int)), dim, simparams...)
    simresult = birth_death_pushing(first(trange); dim = dim, simparams..., seed=seed)
    state = simresult |> r -> (r[:index], r[:mutation], r[:time])
    tumor_stepper!(simresult, trange; simparams...)
end

function record_growth(time_series; path="temp.gif", frames=1,
        shading = false, size = (500,500),
        points_colors = colors_by_mutations
        )
    dim = length(first(first(time_series).position))
    width = maximum(abs.(vcat(last(time_series).position...)))
    origin = fill(width, dim)
    limits = dim == 2 ? FRect2D(-origin, 2*origin) : FRect3D(-origin, 2*origin)

    pointsNode = Node([Point(zeros(dim)...)])
    colorsNode = Node([palette(:tab20)[1]])
    scene = Scene(resolution=size, show_axis=false)
    meshscatter!(scene, pointsNode, color = colorsNode, scale_plot=false, limits = limits, shading = shading, markersize = 1.)
    record(scene, path) do io
        for snapshot in time_series
            pointsNode[], colorsNode[] = points_colors(snapshot)
            for _=1:frames recordframe!(io) end
        end
    end
end
