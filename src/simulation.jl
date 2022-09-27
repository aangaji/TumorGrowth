export birth_death_pushing, birth_death_pushing!, DataFrame

########---FORMAT---########

"""
    Cell(index, position, mutations, parent, b, t_birth, p_birth)

Mutable structure used in the simulation algorithm. function. The output of the simulation
method `birth_death_pushing` is a `Vector{Cell}` and parameters.
"""
mutable struct Cell
    index :: Int64
    position :: Vector{Float64}
    mutations :: Vector{Int64}
    parent :: Int64
    b :: Float64
    t_birth :: Float64
	p_birth :: SVector{dim,Float64} where dim
end
"""
    Cell(;index, position, mutations, parent, b, t_birth, p_birth)

Constructor of type `Cell` with kwrgs for fields.
"""
Cell(;index, position, mutations, parent, b, t_birth, p_birth) = Cell(index, position, mutations, parent, b, t_birth, p_birth)

struct Mutation
	origin :: Int64
	ancestor :: Int64
	t_birth :: Float64
	N_birth :: Int64
	p_birth :: SVector{dim,Float64} where dim
end
Mutation(; origin, ancestor, t_birth, N_birth, p_birth) :: Mutation = Mutation(origin, ancestor, t_birth, N_birth, p_birth)

"""
    DataFrame(tumor::Vector{T})

Generates a DataFrame out of a Vector of type T (Cell). Columnnames of the output
are fieldnames of the type T (Cell).
"""
function DataFrame(tumor::Vector{T}) where T
    df = DataFrame()
    for field in fieldnames(T)
        ftype = fieldtype(T, field)
        df[!,field] = if ftype <: SVector
            collect.(reinterpret.(eltype(ftype),getfield.(tumor, field)))
        else
            getfield.(tumor, field)
        end
    end
    return df
end

############################
######---BIRTH RATE---######

w(r::Float64; σ=1.)  = exp(-r^2/(2*σ^2))/√(2π*σ^2)

"""
    b_linear()
Call before running the simulation to set the birth rate profile.
Linear is set by default. A custom profile can also be set by overwriting TumorGrowth.b_curve(ρ; bup, ρc).
"""
function b_linear()
	@eval b_curve(ρ::Float64; bup::Float64, ρc::Float64) = max(0., bup*(1. - ρ/ρc))
end

"""
    b_hill(n=2)
Call before running the simulation to set the birth rate profile to a hill curve o f order `n`.
A custom profile can also be set by overwriting TumorGrowth.b_curve(ρ; bup, ρc).
"""
function b_hill(n=2)
	@eval b_curve(ρ::Float64; bup::Float64, ρc::Float64) = bup/(1+exp(ρ-ρc)^($n))
end
b_linear()

"""
    update_birthrate!(cell::Cell, tumor; bup::Float64, ρ = Inf)

Determines a cells birth rate based on its surrounding from the set birth rate profile `b_curve`. For each cell in the variable `tumor` a gaussian weight (hardcoded width `σ=4`) is calculated on its distance. The sum of weights is the density that gives the cells birth rate at set denstiy threshold `ρ` and maximal birth rate `bup`.
"""
function update_birthrate!(cell::Cell, tumor; bup::Float64, ρ = Inf)
	isempty(tumor) && return cell.b = bup
    ws = getfield.(tumor, :position) .|> p -> w(norm(cell.position - p); σ=4.)
    return cell.b = b_curve(sum(ws); bup = bup, ρc = ρ)
end

#############################
#######---SIMULATION---######

"""
    birth!(tumor::Vector{Cell}, parent, cur_id, cur_mutation, μ, cellbox, t, dimv::Val{dim}) where dim
Birth event handler in simulation algorithm.
At division the birth! function is called and parent or daughter may or may not gain a new mutation.
"""
function birth!(tumor::Vector{Cell}, mutations::Vector{Mutation}, parent, cur_id, cur_mutation, μ, cellbox, t, dimv::Val{dim}) where {dim}

    δ = rand(dim) .- 0.5
    pos = parent.position + 2.1 .* δ / norm(δ)
    push!(cellbox, pos2box(pos, dimv))

    new = Cell(cur_id, copy(pos), copy(parent.mutations), parent.index, parent.b, t, SVector{dim,Float64}(pos))
    push!(tumor, new)

    m1, m2 = rand(Poisson(μ / 2), 2)
    ancestor_1 = isempty(parent.mutations) ? 0 : last(parent.mutations)
    ancestor_2 = isempty(new.mutations) ? 0 : last(new.mutations)
    for m in 1:m1
        push!(mutations,
            Mutation(parent.index, ancestor_1,
                t, length(tumor), SVector{dim,Float64}(parent.position)
            )
        )
        push!(parent.mutations, cur_mutation + m)
    end
    for m in 1:m2
        push!(mutations,
            Mutation(new.index, ancestor_2,
                t, length(tumor), SVector{dim,Float64}(pos)
            )
        )
        push!(new.mutations, cur_mutation + m1 + m)
    end

    return m1 + m2
end

"""
    birth_death_pushing!( tumor::Vector{Cell}, until;
    b, d, μ, ρ=Inf, dim=length(tumor[1].position), seed=abs(rand(Int)),
    cur_id = 1, cur_mutation = 0, t = 0.0, showprogress=true)

Evolves an existing collection of cells `Vector{Cell}`.
mutates `tumor` variable and returns current `index`, `mutation` and time`t` as a `Dict`.

Rejection algorithm for density dependent tumor growth:
At each iteration a parent is uniformly chosen and its birth rate updated. It then either devides, dies (removed) or nothing happens. Time is incremented appropriatly.
For division the birth! function is called and parent or daughter may or may not gain a new mutation. The pushing! algorithm then removes overlaps.

Call `birth_death_pushing` to run a simulation from a single cell
returns:    `Dict{:tumor, :index, :mutation, :time}`
arg:        `until`       loopcondition, if Float final time, if Int final size
kwargs:     `b, d, μ, ρ=Inf, cur_mutation = 0, dim, seed=nothing`

To further evolve a tumor call `birth_death_pushing!` on existing `Cell` array and set kwargs `cur_id` (next cell index), `cur_mutation` (next mutation), `t` (time) appropriatly.
"""
function birth_death_pushing!( tumor::Vector{Cell}, mutations::Vector{Mutation}, until;
	b, d, μ, ρ=Inf, dim=length(tumor[1].position), seed=nothing,
	cur_id = 1, cur_mutation = 0, t = 0.0, showprogress=true)
    dimv = Val(dim)

    isnothing(seed) || Random.seed!(seed)

    N = length(tumor)

    cellbox = [pos2box(p, dimv) for p in getfield.(tumor,:position)]

    prog = showprogress ? ProgressUnknown("Progress: Tumor size ") : nothing
    while loop_condition(N,t, until)

        N==0 && error("Tumor died")
        row = rand(1:N)
        parent = tumor[row]

        p = rand()*(b+d) - update_birthrate!(parent, view(tumor, find_neighbors(cellbox, row; s=4)); bup=b, ρ = ρ)

        t += randexp()/((d+b)*N)

        if p < 0.
            cur_id += 1
            N += 1
            cur_mutation += birth!(tumor, mutations, parent, cur_id, cur_mutation, μ, cellbox, t, dimv)

            pushing!(tumor, N, cellbox, dimv)
        elseif p < d
            N -= 1
            deleteat!.( (tumor, cellbox), row)
        end
        showprogress && ProgressMeter.update!(prog, N )
    end
	for (i, cell) in enumerate(tumor)
	    update_birthrate!(cell, view(tumor, find_neighbors(cellbox, i; s=4)); bup=b, ρ = ρ)
	end
    return Dict{Symbol, Any}(:index => cur_id, :mutation => cur_mutation, :time => t)
end
birth_death_pushing!( until; tumor::Vector{Cell}, mutations::Vector{Mutation}, simparams...) =  birth_death_pushing!( tumor, mutations, until; simparams...)

loop_condition(N::Int,t::Float64, Nfinal::Int) = N < Nfinal
loop_condition(N::Int,t::Float64, tfinal::Float64) = t < tfinal


"""
    function birth_death_pushing( until; b, d, μ, ρ=Inf, cur_mutation = 0, dim , seed=abs(rand(Int)), showprogress=true)

Rejection algorithm for density dependent tumor growth:
At each iteration a parent is uniformly chosen and its birth rate updated. It then either devides, dies (removed) or nothing happens. Time is incremented appropriatly.
For division the birth! function is called and parent or daughter may or may not gain a new mutation. The pushing! algorithm then removes overlaps.

Call `birth_death_pushing` to run a simulation from a single cell
returns:    `Dict{:tumor, :index, :mutation, :time}`
arg:        `until`       loopcondition, if Float final time, if Int final size
kwargs:     `b, d, μ, ρ=Inf, cur_mutation = 0, dim, seed=nothing`

The output `:tumor` is a `Vector{Cell}` which can be passed to `DataFrame` for analysis or plotting.

To further evolve a tumor call `birth_death_pushing!` on existing `Cell` vector and set kwargs `cur_id` (next cell index), `cur_mutation` (next mutation), `t` (time) appropriatly.
"""
function birth_death_pushing( until; b, cur_mutation = 0, dim, simparams...)
    p₀ = zeros(Float64,dim)
    tumor = [Cell(; index=1, position=p₀, parent=0, mutations=collect(1:cur_mutation), b=b, t_birth=0.0, p_birth=SVector{dim,Float64}(p₀))]

	mutations = [Mutation(; origin=1, ancestor=0, t_birth=0., N_birth=1, p_birth=SVector{dim,Float64}(p₀)) for m=1:cur_mutation]

    output = birth_death_pushing!(tumor, mutations, until; b=b, cur_mutation=cur_mutation, dim=dim, simparams...)
	output[:tumor] = tumor
	output[:mutations] = mutations

	return output
end
