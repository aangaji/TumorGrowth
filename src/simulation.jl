export birth_death_pushing, birth_death_pushing!, DataFrame

########---FORMAT---########

#=
	The simulation uses and returns an array of the mutable Cell type. To obtain a DataFrame with the same fields for analysis or plotting call the overloaded DataFrame function on the array.
=#
mutable struct Cell
    index :: Int64
    position :: Vector{Float64}
    mutations :: Vector{Int64}
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

w(r; σ=1.)  = exp(-r^2/(2*σ^2))/√(2π*σ^2)

#=
	Call one of the following before running the simulation to set the birth rate profile.
	Linear is set by default. A custom profile can also be set by overwriting TumorGrowth.b_curve(ρ; bup, ρc).

	A cell's birth rate is updated by computing the local density as a weighted (Gaussian(0, 4.)) sum of distances to a given neighborhood and passing it to the birth rate profile.
=#

function b_linear()
	@eval b_curve(ρ; bup = log(2), ρc) = max(0., bup*(1. - ρ/ρc))
end
function b_hill(n=2)
	@eval b_curve(ρ; bup = log(2), ρc) = bup/(1+exp(ρ-ρc)^($n))
end
b_linear()

function update_birthrate!(cell::Cell, tumor; bup=log(2), ρc = Inf)
	isempty(tumor) && return cell.b = bup
    ws = getfield.(tumor, :position) .|> p -> w(norm(cell.position - p); σ=4.)
    return cell.b = b_curve(sum(ws); bup = bup, ρc = ρc)
end

#############################
#######---SIMULATION---######

#=
	Rejection algorithm for density dependent tumor growth:
	At each iteration a parent is uniformly chosen and its birth rate updated. It then either devides, dies (removed) or nothing happens. Time is incremented appropriatly.
	For division the birth! function is called and parent or daughter may or may not gain a new mutation. The pushing! algorithm then removes overlaps.

	Call birth_death_pushing to run a simulation from a single cell
	returns:	Dict{:tumor, :index, :mutation, :time}
	arg:		until       loopcondition, if Float final time, if Int final size
	kwargs:	  	b, d, μ, ρc=Inf, cur_mutation = 0, dim, seed=nothing

	To further evolve a tumor call birth_death_pushing! on existing Cell array and set kwargs cur_id (next cell index), cur_mutation (next mutation), t (time) appropriatly.
=#

function birth!(tumor::Vector{Cell}, parent, cur_id, cur_mutation, μ, cellbox, t, dimv::Val{dim}) where dim

    δ = rand(dim) .- 0.5
    pos = parent.position + 2.1 .* δ/norm(δ)
    push!(cellbox, pos2box(pos, dimv))

	new = Cell(cur_id, copy(pos), copy(parent.mutations), parent.index, parent.b, t, SVector{dim,Float64}(pos))
	push!(tumor, new)

	m = 0
	if rand() < μ/2
		m += 1
		push!(parent.mutations, cur_mutation+m)
	end
	if rand() < μ/2
		m += 1
		push!(new.mutations, cur_mutation+m)
	end
	return m
end

function birth_death_pushing!( tumor::Vector{Cell}, times, Ns, until; b, d, μ, ρc=Inf, cur_id = 1, cur_mutation = 0, t = 0.0, dim=length(tumor[1].position), seed=nothing)
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

        t += randexp()/((d+b)*N)

        if p < 0.
            cur_id += 1
            N += 1
            cur_mutation += birth!(tumor, parent, cur_id, cur_mutation, μ, cellbox, t, dimv)

            pushing!(tumor, N, cellbox, dimv)

        elseif p < d
            N -= 1
            deleteat!.( (tumor, cellbox), row)
        end
		push!(times, t)
		push!(Ns, N)

        #print("Progress: $(size(tumor,1))/$tumor_size | Recent mutation: $cur_mutation \r")
        #@info "Tumor size" progress=N/tumor_size _id=id
        ProgressMeter.update!(prog, N )
    end
	for (i, cell) in enumerate(tumor)
	    update_birthrate!(cell, view(tumor, find_neighbors(cellbox, i; s=4)); bup=b, ρc = ρc)
	end
    return Dict{Symbol, Any}(:index => cur_id, :mutation => cur_mutation, :times => times, :sizes => Ns)
end

loop_condition(N,t, Nfinal::Int) = N < Nfinal
loop_condition(N,t, tfinal::Float64) = t < tfinal


function birth_death_pushing( until; b, d, μ, ρc=Inf, cur_mutation = 0, dim , seed=nothing)
    p₀ = zeros(Float64,dim)
    tumor = [Cell(index=1, position=p₀, parent=0, mutations=collect(1:cur_mutation), b=b, t_birth=0.0, p_birth=SVector{dim,Float64}(p₀))]

	times = [0.0]
	Ns = [1]

    output = birth_death_pushing!(tumor, times, Ns, until; b=b, d=d, μ=μ, ρc=ρc, cur_mutation=cur_mutation, dim=dim, seed=seed)
	output[:tumor] = tumor

	return output
end
