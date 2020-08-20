export birth_death_pushing, birth_death_pushing!, DataFrame!, DataFrame

mutable struct Cell
    index :: Int64
    position :: Vector{Float64}
    parent :: Int64
    mutations :: Vector{Int64}
    t_birth :: Float64
    p_birth :: Vector{Float64}
end

struct Mutation_event
    mutation :: Int64
    parent :: Int64
    t_birth :: Float64
    p_birth :: Vector{Float64}
end

DataFrame(tumor::Vector{T}) where T = DataFrame( [getfield.(tumor, field) for field in fieldnames(T)], collect(fieldnames(T)) )

function birth!(tumor::Vector{Cell}, parent, cur_id, cur_mutation, mu, cellbox, t, dim)

    δ = rand(dim) .- 0.5
    pos = parent.position + 2 .* δ/norm(δ)
    push!(cellbox, pos2box.(pos))

    new = Cell(cur_id, copy(pos), parent.index, copy(parent.mutations), t, copy(pos))

    push!(tumor, new)

    if rand() < mu
        push!(new.mutations, cur_mutation+1)
        return Mutation_event(cur_mutation+1, cur_id, t, pos)
    end
end

function birth_death_pushing!( tumor::Vector{Cell}, mutation_events::Vector{Mutation_event}, tumor_size; b, d, mu, cur_id = 1, cur_mutation = 0, t = 0.0, dim=length(tumor[1].position), seed=nothing)
    isnothing(seed) || Random.seed!(seed)

    N = length(tumor)

    cellbox = [pos2box.(p) for p in getfield.(tumor,:position)]

    b_prob = 1/(1 + d/b)

    #Juno.progress() do id
    prog = ProgressUnknown("Progress: Tumor size ")
    while N < tumor_size
            N==0 && error("Tumor died")
            parent_row = rand(1:N)
            parent = tumor[parent_row]

            t += randexp()/((b+d)*N)

            if rand() <= b_prob
                cur_id += 1
                N += 1
                mutation_event = birth!(tumor, parent, cur_id, cur_mutation, mu, cellbox, t, dim)

                if !isnothing(mutation_event)
                    cur_mutation += 1
                    push!(mutation_events, mutation_event)
                end
                pushing!(tumor, N, cellbox)
            else
                N -= 1
                deleteat!.( (tumor, cellbox), parent_row)
            end
        #print("Progress: $(size(tumor,1))/$tumor_size | Recent mutation: $cur_mutation \r")
        #@info "Tumor size" progress=N/tumor_size _id=id
        ProgressMeter.update!(prog, N)
    end
    return (cur_id, cur_mutation, t)
end

function birth_death_pushing( tumor_size; b, d, mu, cur_mutation = 0, dim , seed=nothing)

    tumor = [Cell(1, zeros(Float64,dim), 0, collect(1:cur_mutation), 0.0, zeros(Float64, dim))]

    mutation_events = Vector{Mutation_event}()

    cur_id, cur_mutation, t = birth_death_pushing!(tumor, mutation_events, tumor_size; b=b, d=d, mu=mu, cur_mutation=cur_mutation, dim=dim, seed=seed)

    return (cur_id, cur_mutation, t), tumor, mutation_events
end
