export tumor_stepper, tumor_stepper!, record_growth

function tumor_stepper!(params, tumor, mut_events, trange; b, d, mu)
    id, mut_id, t = params
    dim = length(first(tumor).position)
    time_series = Vector{DataFrame}()
    for tp in trange
        id, mut_id, t = birth_death_pushing!(tumor, mut_events, tp; b=b, d=d, mu=mu, cur_id=id, cur_mutation=mut_id, t=t, dim=dim)
        push!(time_series, deepcopy(DataFrame(tumor)) )
    end
    return time_series
end
function tumor_stepper(trange; b, d, mu, dim, seed = nothing)
    (params, tumor, mut_events) = birth_death_pushing(trange[1]; b=b, d=d, mu=mu, dim=dim, seed=seed)
    tumor_stepper!(params, tumor, mut_events, trange; b=b, d=d, mu=mu)
end


function record_growth(time_series; path="test.gif", manipulate! = nothing)
    origin = maximum( last(time_series).position |> p-> ( getindex.(p,i) for i=1:length(first(p)) ) .|> x-> maximum(abs.(x)) ) - 1
    limits = FRect(-origin, -origin, 2*origin, 2*origin)

    scene = plotting_colored_mutations(first(time_series), limits=limits)
    record(scene, path) do io
        for snapshot in time_series
            while length(scene) > 1
                delete!(scene, scene[end])
            end
            plotting_colored_mutations!(scene, snapshot)
            isnothing(manipulate!) || manipulate!(scene, snapshot)
            recordframe!(io)
        end
    end
end
