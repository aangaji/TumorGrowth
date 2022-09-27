export nonspatial, nonspatial!

function birth!(tumor::Vector{Cell}, mutations::Vector{Mutation}, parent, cur_id, cur_mutation, μ, t)

    pos = Float64[]

	new = Cell(cur_id, pos, copy(parent.mutations), parent.index, parent.b, t, SVector{0,Float64}())
	push!(tumor, new)

	m1, m2 = rand(Poisson(μ/2), 2)
    ancestor_1 = isempty(parent.mutations) ? 0 : last(parent.mutations)
    ancestor_2 = isempty(new.mutations) ? 0 : last(new.mutations)
    for m in 1:m1
        push!(mutations, 
            Mutation(parent.index, ancestor_1,
                t, length(tumor), SVector{dim,Float64}(parent.position)
                )  
            )
        push!(parent.mutations, cur_mutation+m)
    end
    for m in 1:m2
        push!(mutations, 
            Mutation(new.index, ancestor_2,
                t, length(tumor), SVector{dim,Float64}(pos)
                )
            )
        push!(new.mutations, cur_mutation+m1+m)
    end

	return m1+m2
end

function nonspatial!( tumor::Vector{Cell}, mutations::Vector{Mutation}, until; b, d, μ, cur_id = 1, cur_mutation = 0, t = 0.0, seed=nothing, showprogress=true)

    isnothing(seed) || Random.seed!(seed)

    N = length(tumor)

    prog = showprogress ? ProgressUnknown("Progress: Tumor size ") : nothing
    while loop_condition(N,t, until)

        N==0 && error("Tumor died")
        row = rand(1:N)
        parent = tumor[row]

        p = rand()*(b+d) - b
        t += randexp()/((d+b)*N)

        if p < 0.
            cur_id += 1
            N += 1
            cur_mutation += birth!(tumor, mutations, parent, cur_id, cur_mutation, μ, t)
        elseif p < d
            N -= 1
            deleteat!(tumor, row)
        end
        showprogress && ProgressMeter.update!(prog, N )
    end
    return Dict{Symbol, Any}(:index => cur_id, :mutation => cur_mutation, :time => t)
end
nonspatial!( until; tumor::Vector{Cell}, mutations::Vector{Mutation}, simparams...) =  nonspatial!( tumor, mutations, until; simparams...)


function nonspatial( until; cur_mutation = 0, b, simparams...)
    p₀ = Float64[]
    tumor = [Cell(index=1, position=p₀, parent=0, mutations=collect(1:cur_mutation), b=b, t_birth=0.0, p_birth=SVector{0,Float64}(p₀))]

	mutations = [Mutation(; origin=1, ancestor=0, t_birth=0., N_birth=1, p_birth=SVector{0,Float64}(p₀)) for m=1:cur_mutation]

    output = nonspatial!(tumor, mutations, until; cur_mutation=cur_mutation, b=b, simparams...)
	output[:tumor] = tumor
	output[:mutations] = mutations

	return output
end
