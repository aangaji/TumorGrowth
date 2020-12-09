export tumor_stepper, tumor_stepper!, record_growth

function tumor_stepper!(params, tumor, trange; b, d, μ, ρc=ρc)
    id, mut_id, t = params
    dim = length(first(tumor).position)
    time_series = Vector{DataFrame}()
    for tp in trange
        up = birth_death_pushing!(tumor, tp; b=b, d=d, μ=μ, ρc=ρc, cur_id=id, cur_mutation=mut_id, t=t, dim=dim)
        push!(time_series, deepcopy(DataFrame(tumor)) )
        id, mut_id, t = up[:index], up[:mutation], up[:time]
    end
    return time_series
end
function tumor_stepper(trange; b, d, μ, ρc, dim, seed = nothing)
    simresult = birth_death_pushing(first(trange); b=b, d=d, μ=μ, ρc=ρc, dim=dim, seed=seed)
    params = simresult |> r -> (r[:index], r[:mutation], r[:time])
    tumor_stepper!(params, simresult[:tumor], trange; b=b, d=d, μ=μ, ρc=ρc)
end

function record_growth(time_series; path="test.gif", frames=1,
        shading = false,
        points_colors = colors_by_mutations
        )
    dim = length(first(first(time_series).position))
    width = maximum(abs.(vcat(last(time_series).position...)))
    origin = fill(width, dim)
    limits = dim == 2 ? FRect2D(-origin, 2*origin) : FRect3D(-origin, 2*origin)

    pointsNode = Node([Point(zeros(dim)...)])
    colorsNode = Node([palette(:tab20)[1]])
    scene = meshscatter(pointsNode, color = colorsNode, scale_plot=false, limits = limits, shading = shading, markersize = 1.)
    record(scene, path) do io
        for snapshot in time_series
            pointsNode[], colorsNode[] = points_colors(snapshot)
            for _=1:frames recordframe!(io) end
        end
    end
end
