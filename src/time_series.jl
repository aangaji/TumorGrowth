export tumor_stepper, tumor_stepper!, record_growth

function tumor_stepper!(params, tumor, trange; b, d, mu)
    id, mut_id, t = params
    dim = length(first(tumor).position)
    time_series = Vector{DataFrame}()
    for tp in trange
        id, mut_id, t = birth_death_pushing!(tumor, tp; b=b, d=d, mu=mu, cur_id=id, cur_mutation=mut_id, t=t, dim=dim)
        push!(time_series, deepcopy(DataFrame(tumor)) )
    end
    return time_series
end
function tumor_stepper(trange; b, d, mu, dim, seed = nothing)
    (params, tumor) = birth_death_pushing(first(trange); b=b, d=d, mu=mu, dim=dim, seed=seed)
    tumor_stepper!(params, tumor, trange; b=b, d=d, mu=mu)
end


function record_growth(time_series; path="test.gif", plot_func! = plotting_colored_mutations!, frames=3)
    origin = maximum( last(time_series).position |> p-> ( getindex.(p,i) for i=1:length(first(p)) ) .|> x-> maximum(abs.(x)) ) - 1
    limits = FRect(-origin, -origin, 2*origin, 2*origin)

    scene = meshscatter([Point(origin...)], scale_plot=false, limits = limits)
    record(scene, path) do io
        for snapshot in time_series
            while length(scene) > 1
                delete!(scene, scene[end])
            end
            plot_func!(scene, snapshot)
            for _=1:frames recordframe!(io) end
        end
    end
end
