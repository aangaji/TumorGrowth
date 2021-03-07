export nonspatial, nonspatial!

function birth!(tumor::Vector{Cell}, parent, cur_id, cur_mutation, μ, t, dimv::Val{dim}) where dim

    pos = zeros(Float64, dim)

	new = Cell(cur_id, pos, copy(parent.mutations), parent.index, parent.b, t, SVector{dim,Float64}(pos))
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

function nonspatial!( tumor::Vector{Cell}, until; b, d, μ, cur_id = 1, cur_mutation = 0, t = 0.0, dim=length(tumor[1].position), seed=nothing)
    dimv = Val(dim)

    isnothing(seed) || Random.seed!(seed)

    N = length(tumor)

    prog = ProgressUnknown("Progress: Tumor size ")
    while loop_condition(N,t, until)

        N==0 && error("Tumor died")
        row = rand(1:N)
        parent = tumor[row]

        p = rand()*(b+d) - b
        t += randexp()/((d+b)*N)

        if p < 0.
            cur_id += 1
            N += 1
            cur_mutation += birth!(tumor, parent, cur_id, cur_mutation, μ, t, dimv)
        elseif p < d
            N -= 1
            deleteat!(tumor, row)
        end
        ProgressMeter.update!(prog, N )
    end
    return Dict{Symbol, Any}(:index => cur_id, :mutation => cur_mutation, :time => t)
end


function nonspatial( until; b, d, μ, cur_mutation = 0, dim , seed=nothing)
    p₀ = zeros(Float64,dim)
    tumor = [Cell(index=1, position=p₀, parent=0, mutations=collect(1:cur_mutation), b=b, t_birth=0.0, p_birth=SVector{dim,Float64}(p₀))]

    output = nonspatial!(tumor, until; b=b, d=d, μ=μ, cur_mutation=cur_mutation, dim=dim, seed=seed)
	output[:tumor] = tumor

	return output
end
