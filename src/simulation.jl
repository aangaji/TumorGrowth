export birth_death_pushing, birth_death_pushing!, DataFrame

########---FORMAT---########

"""
    Cell(index, position, mutations, parent, b, d, mu, rho, t_birth)

Mutable structure used in the simulation algorithm. function. The output of the simulation
method `birth_death_pushing` is a `Vector{Cell}` and parameters.
"""
mutable struct Cell
    index :: Int64
    position :: Vector{Float64}
    mutations :: Vector{Int64}
    parent :: Int64
    b :: Float64
    d :: Float64
    mu :: Float64
    rho :: Float64
    t_birth :: Float64
end
"""
    Cell(;index, position, mutations, parent, b, d, mu, rho, t_birth)

Constructor of type `Cell` with kwrgs for fields.
"""
Cell(;index, position, mutations, parent, b, d, mu, rho, t_birth) = Cell(index, position, mutations, parent, b, d, mu, rho , t_birth)

struct Mutation
	origin :: Int64
	ancestor :: Int64
	t_birth :: Float64
	N_birth :: Int64
end
Mutation(; origin, ancestor, t_birth, N_birth) :: Mutation = Mutation(origin, ancestor, t_birth, N_birth)

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

w(r::Float64; sigma=1.)  = exp(-r^2/(2*sigma^2))/√(2π*sigma^2)

b_linear(rho::Float64; b_max::Float64, rho_c::Float64) = max(0., b_max*(1. - rho/rho_c))

"""
    get_birthrate(cell::Cell, tumor; b_of_rho)

Determines a cells birth rate based on its surrounding from the birth rate profile `b_of_rho`. For each cell in the variable `tumor` a gaussian weight (hardcoded width `sigma=4`) is calculated on its distance. The sum of weights is the density that gives the cells birth rate at its density threshold `rho` and maximal birth rate `b`.
"""
function get_birthrate(cell::Cell, neighbors; b_of_rho )
	isempty(neighbors) && return cell.b
    ws = getfield.(neighbors, :position) .|> p -> w(norm(cell.position - p); sigma=4.)
    return b_of_rho(sum(ws); b_max = cell.b, rho_c = cell.rho)
end

#############################
#######---SIMULATION---######

"""
    birth!(tumor::Vector{Cell}, parent, cur_id, cur_mutation, t)
Birth event handler in simulation algorithm.
At division the birth! function is called and parent or daughter may or may not gain a new mutation.
"""
function birth!(tumor::Vector{Cell}, mutations::Vector{Mutation}, parent, cur_id, cur_mutation, t)

    new = deepcopy(parent)
    new.index = cur_id
    new.t_birth = t
    push!(tumor, new)

    m1, m2 = rand(Poisson(parent.mu / 2), 2)
    ancestor_1 = isempty(parent.mutations) ? 0 : last(parent.mutations)
    ancestor_2 = isempty(new.mutations) ? 0 : last(new.mutations)
    for m in 1:m1
        push!(mutations,
            Mutation(parent.index, ancestor_1, t, length(tumor)
            )
        )
        push!(parent.mutations, cur_mutation + m)
    end
    for m in 1:m2
        push!(mutations,
            Mutation(new.index, ancestor_2, t, length(tumor)
            )
        )
        push!(new.mutations, cur_mutation + m1 + m)
    end

    return m1 + m2
end

"""
    birth_death_pushing!( tumor::Vector{Cell}, until;
    dim=length(tumor[1].position), seed=abs(rand(Int)),
    cur_id = 1, cur_mutation = 0, t = 0.0, showprogress=true)

Evolves an existing collection of cells `Vector{Cell}`.
mutates `tumor` variable and returns current `index`, `mutation` and time`t` as a `Dict`.

Rejection algorithm for density dependent tumor growth:
At each iteration a parent is uniformly chosen and its birth rate updated. It then either devides, dies (removed) or nothing happens. Time is incremented appropriatly.
For division the birth! function is called and parent or daughter may or may not gain a new mutation. The pushing! algorithm then removes overlaps.

Call `birth_death_pushing` to run a simulation from a single cell
returns:    `Dict{:tumor, :index, :mutation, :time}`
arg:        `until`       loopcondition, if Float final time, if Int final size
kwargs:     `cur_mutation = 0, dim, seed=nothing`

To further evolve a tumor call `birth_death_pushing!` on existing `Cell` array and set kwargs `cur_id` (next cell index), `cur_mutation` (next mutation), `t` (time) appropriatly.
"""
function birth_death_pushing!( tumor::Vector{Cell}, mutations::Vector{Mutation}, until; b_of_rho = b_linear,
	dim=length(tumor[1].position), seed=nothing, cur_id = 1, cur_mutation = 0, t = 0.0, showprogress=true)
    dimv = Val(dim)

    b_max = maximum(getfield.(tumor, :b))
    d_max = maximum(getfield.(tumor, :d))

    isnothing(seed) || Random.seed!(seed)

    N = length(tumor)

    cellbox = [pos2box(p, dimv) for p in getfield.(tumor,:position)]

    prog = showprogress ? ProgressUnknown("Progress: Tumor size ") : nothing
    while loop_condition(N,t, until)

        N==0 && error("Tumor died")
        row = rand(1:N)
        parent = tumor[row]

        p = rand()*(b_max+d_max) - get_birthrate(parent, view(tumor, find_neighbors(cellbox, row; s=4)); b_of_rho=b_of_rho)

        t += randexp()/((b_max+d_max)*N)

        if p < 0.
            cur_id += 1
            N += 1
            cur_mutation += birth!(tumor, mutations, parent, cur_id, cur_mutation, t)
            
            del = rand(dim) .- 0.5
            last(tumor).position = parent.position + 2.1 .* del / norm(del)
            push!(cellbox, pos2box(last(tumor).position, dimv))

            pushing!(tumor, N, cellbox, dimv)
        elseif p < parent.d
            N -= 1
            deleteat!.( (tumor, cellbox), row)
        end
        showprogress && ProgressMeter.update!(prog, N )
    end

    return Dict{Symbol, Any}(:index => cur_id, :mutation => cur_mutation, :time => t)
end
birth_death_pushing!( until; tumor::Vector{Cell}, mutations::Vector{Mutation}, simparams...) =  birth_death_pushing!( tumor, mutations, until; simparams...)

loop_condition(N::Int,t::Float64, Nfinal::Int) = N < Nfinal
loop_condition(N::Int,t::Float64, tfinal::Float64) = t < tfinal


"""
    function birth_death_pushing( until; b, d, mu, rho=Inf, cur_mutation = 0, dim , seed=abs(rand(Int)), showprogress=true)

Rejection algorithm for density dependent tumor growth:
At each iteration a parent is uniformly chosen and its birth rate updated. It then either devides, dies (removed) or nothing happens. Time is incremented appropriatly.
For division the birth! function is called and parent or daughter may or may not gain a new mutation. The pushing! algorithm then removes overlaps.

Call `birth_death_pushing` to run a simulation from a single cell
returns:    `Dict{:tumor, :index, :mutation, :time}`
arg:        `until`       loopcondition, if Float final time, if Int final size
kwargs:     `b, d, mu, rho=Inf, cur_mutation = 0, dim, seed=nothing`

The output `:tumor` is a `Vector{Cell}` which can be passed to `DataFrame` for analysis or plotting.

To further evolve a tumor call `birth_death_pushing!` on existing `Cell` vector and set kwargs `cur_id` (next cell index), `cur_mutation` (next mutation), `t` (time) appropriatly.
"""
function birth_death_pushing( until; b, d, mu, rho=Inf, cur_mutation = 0, dim, simparams...)

    tumor = [Cell(; index=1, position=zeros(Float64,dim), parent=0, mutations=collect(1:cur_mutation), b=b, d=d, mu=mu, rho=rho, t_birth=0.)]

	mutations = [Mutation(; origin=1, ancestor=0, t_birth=0., N_birth=1) for _=1:cur_mutation]

    output = birth_death_pushing!(tumor, mutations, until; cur_mutation=cur_mutation, dim=dim, simparams...)
	output[:tumor] = tumor
	output[:mutations] = mutations

	return output
end
