export birth_death_pushing, birth_death_pushing!, DataFrame

mutable struct Cell
    index :: Int64
    position :: Vector{Float64}
    parent :: Int64
    mutations :: SVector{len, Int64} where len
    t_birth :: Float64
    p_birth :: Vector{Float64}
end

Cell(;index, position, mutations, parent, b, t_birth) = Cell(index, position, mutations, parent, b, t_birth)

struct Mutation_event
    mutation :: Int64
    parent :: Int64
    t_birth :: Float64
    p_birth :: Vector{Float64}
end

function DataFrame(tumor::Vector{T}) where T
    df = DataFrame()
    for field in fieldnames(T)
        ftype = fieldtype(Cell, field)
        df[!,field] = if ftype <: SVector
            collect.(reinterpret.(eltype(ftype),getfield.(tumor, field)))
        else
            getfield.(tumor, field)
        end
    end
    return df
end

function birth!(tumor::Vector{Cell}, parent, cur_id, cur_mutation, mu, cellbox, t, dimv::Val{dim}) where dim

    δ = rand(dim) .- 0.5
    pos = parent.position + 2.1 .* δ/norm(δ)
    push!(cellbox, pos2box(pos, dimv))

    if rand() < mu
        new = Cell(cur_id, copy(pos), parent.index, push(parent.mutations, cur_mutation+1), t, copy(pos))
        push!(tumor, new)
        return true
    else
        new = Cell(cur_id, copy(pos), parent.index, parent.mutations, t, copy(pos))
        push!(tumor, new)
        return false
    end
    # return Mutation_event(cur_mutation+1, cur_id, t, pos)
end

function birth_death_pushing!( tumor::Vector{Cell}, mutation_events::Vector{Mutation_event}, until; b, d, mu, cur_id = 1, cur_mutation = 0, t = 0.0, dim=length(tumor[1].position), seed=nothing)
    dimv = Val(dim)

    isnothing(seed) || Random.seed!(seed)

    N = length(tumor)

    cellbox = [pos2box(p, dimv) for p in getfield.(tumor,:position)]

    #Juno.progress() do id
    prog = ProgressUnknown("Progress: Tumor size ")
    while loop_condition(N,t, until)

        N==0 && error("Tumor died")
        parent_row = rand(1:N)
        parent = tumor[parent_row]

        p = rand()*(b+d) - b

        t += randexp()/((b+d)*N)

        if p < 0.
            cur_id += 1
            N += 1
            mutation_event = birth!(tumor, parent, cur_id, cur_mutation, mu, cellbox, t, dimv)

            mutation_event && (cur_mutation += 1)
            #push!(mutation_events, mutation_event)

            pushing!(tumor, N, cellbox, dimv)
        else
            N -= 1
            deleteat!.( (tumor, cellbox), parent_row)
        end
        #print("Progress: $(size(tumor,1))/$tumor_size | Recent mutation: $cur_mutation \r")
        #@info "Tumor size" progress=N/tumor_size _id=id
        ProgressMeter.update!(prog, N )
    end
    return (cur_id, cur_mutation, t)
end

# loop conditons
loop_condition(N,t, Nfinal::Int) = N < Nfinal
loop_condition(N,t, tfinal::Float64) = t < tfinal


function birth_death_pushing( until; b, d, mu, cur_mutation = 0, dim , seed=nothing)
    p₀ = zeros(Float64,dim)
    tumor = [Cell(1, p₀, 0, SVector{cur_mutation, Int64}(1:cur_mutation...), 0.0, p₀)]

    mutation_events = Vector{Mutation_event}()

    cur_id, cur_mutation, t = birth_death_pushing!(tumor, mutation_events, until; b=b, d=d, mu=mu, cur_mutation=cur_mutation, dim=dim, seed=seed)

    return (cur_id, cur_mutation, t), tumor, mutation_events
end
