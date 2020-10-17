export birth_death_pushing, birth_death_pushing!, DataFrame

mutable struct Cell
    index :: Int64
    position :: Vector{Float64}
    mutations :: SVector{len, Int64} where len
    parent :: Int64
    b :: Float64
    t_birth :: Float64
	p_birth :: SVector{dim,Float64} where dim
end

Cell(;index, position, mutations, parent, b, t_birth, p_birth) = Cell(index, position, mutations, parent, b, t_birth, p_birth)

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

w(r; σ=1.)  = exp(-r^2/(2*σ^2))/√(2π*σ^2)

b_curve(ρ; bup = log(2), ρc) = max(0., bup*(1. - ρ/ρc))

function update_birthrate!(cell::Cell, tumor; bup=log(2), ρc = 20.)
	isempty(tumor) && return cell.b = bup
    ws = getfield.(tumor, :position) .|> p -> w(norm(cell.position - p); σ=4.)
    return cell.b = b_curve(sum(ws); ρc = ρc)
end

function birth!(tumor::Vector{Cell}, parent, cur_id, cur_mutation, mu, cellbox, t, dimv::Val{dim}) where dim

    δ = rand(dim) .- 0.5
    pos = parent.position + 2.1 .* δ/norm(δ)
    push!(cellbox, pos2box(pos, dimv))

    if rand() < mu
        if rand(Bool)
        	new = Cell(cur_id, copy(pos), push(parent.mutations, cur_mutation+1), parent.index, parent.b, t, SVector{dim,Float64}(pos))
		else
			new = Cell(cur_id, copy(pos), parent.mutations, parent.index, parent.b, t, SVector{dim,Float64}(pos))
			parent.mutations = push(parent.mutations, cur_mutation+1)
		end
        push!(tumor, new)
        return true
    else
        new = Cell(cur_id, copy(pos), parent.mutations, parent.index, parent.b, t, SVector{dim,Float64}(pos))
        push!(tumor, new)
        return false
    end
    # return Mutation_event(cur_mutation+1, cur_id, t, pos)
end

function birth_death_pushing!( tumor::Vector{Cell}, until; b, d, mu, ρc=Inf, cur_id = 1, cur_mutation = 0, t = 0.0, dim=length(tumor[1].position), seed=nothing)
    dimv = Val(dim)

    isnothing(seed) || Random.seed!(seed)

    N = length(tumor)

    cellbox = [pos2box(p, dimv) for p in getfield.(tumor,:position)]

    #Juno.progress() do id
    prog = ProgressUnknown("Progress: Tumor size ")
    while loop_condition(N,t, until)

        N==0 && error("Tumor died")
        row = rand(1:N)
        parent = tumor[row]

        p = rand()*(b+d) - update_birthrate!(parent, view(tumor, find_neighbors(cellbox, row; s=4)); bup=b, ρc = ρc)

        t += randexp()/(d*N+ sum(getfield.(tumor,:b)))

        if p < 0.
            cur_id += 1
            N += 1
            mutation_event = birth!(tumor, parent, cur_id, cur_mutation, mu, cellbox, t, dimv)

            mutation_event && (cur_mutation += 1)

            pushing!(tumor, N, cellbox, dimv)
        elseif p < d
            N -= 1
            deleteat!.( (tumor, cellbox), row)
        end
        #print("Progress: $(size(tumor,1))/$tumor_size | Recent mutation: $cur_mutation \r")
        #@info "Tumor size" progress=N/tumor_size _id=id
        ProgressMeter.update!(prog, N )
    end
	for (i, cell) in enumerate(tumor)
	    update_birthrate!(cell, view(tumor, find_neighbors(cellbox, i; s=4)); bup=b, ρc = ρc)
	end
    return (cur_id, cur_mutation, t)
end

# loop conditons
loop_condition(N,t, Nfinal::Int) = N < Nfinal
loop_condition(N,t, tfinal::Float64) = t < tfinal


function birth_death_pushing( until; b, d, mu, ρc=Inf, cur_mutation = 0, dim , seed=nothing)
    p₀ = zeros(Float64,dim)
    tumor = [Cell(index=1, position=p₀, parent=0, mutations=SVector{cur_mutation, Int64}(1:cur_mutation...), b=b, t_birth=0.0, p_birth=SVector{dim,Float64}(p₀))]

    cur_id, cur_mutation, t = birth_death_pushing!(tumor, until; b=b, d=d, mu=mu, ρc=ρc, cur_mutation=cur_mutation, dim=dim, seed=seed)

    return (cur_id, cur_mutation, t), tumor
end
