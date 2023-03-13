export nonspatial, nonspatial!

function nonspatial!(tumor::Vector{Cell}, mutations::Vector{Mutation}, until;
    cur_id=1, cur_mutation=0, t=0.0, seed=nothing, showprogress=true)

    isnothing(seed) || Random.seed!(seed)

    b_max = maximum(getfield.(tumor, :b))
    d_max = maximum(getfield.(tumor, :d))

    N = length(tumor)

    prog = showprogress ? ProgressUnknown("Progress: Tumor size ") : nothing
    while loop_condition(N, t, until)

        N == 0 && error("Tumor died")
        row = rand(1:N)
        parent = tumor[row]

        p = rand() * (b_max+d_max) - parent.b
        t += randexp() / ((b_max+d_max) * N)

        if p < 0.0
            cur_id += 1
            N += 1
            cur_mutation += birth!(tumor, mutations, parent, cur_id, cur_mutation, t)
        elseif p < parent.d
            N -= 1
            deleteat!(tumor, row)
        end
        showprogress && ProgressMeter.update!(prog, N)
    end
    return Dict{Symbol,Any}(:index => cur_id, :mutation => cur_mutation, :time => t)
end
nonspatial!( until; tumor::Vector{Cell}, mutations::Vector{Mutation}, simparams...) =  nonspatial!( tumor, mutations, until; simparams...)


function nonspatial( until; cur_mutation = 0, b, d, mu, simparams...)

    tumor = [Cell(index=1, position=Float64[], parent=0, mutations=collect(1:cur_mutation), b=b, d=d, mu=mu, rho=Inf, t_birth=0.0)]

	mutations = [Mutation(; origin=1, ancestor=0, t_birth=0., N_birth=1) for _=1:cur_mutation]

    output = nonspatial!(tumor, mutations, until; cur_mutation=cur_mutation, simparams...)
	output[:tumor] = tumor
	output[:mutations] = mutations

	return output
end
