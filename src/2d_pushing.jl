export birth_death_pushing

#function for creating a new cell and moving all other cells accordingly to avoid overlapping
function birth_pushing!(tumor::DataFrame, cell::DataFrameRow, index, current_mut::Int64, mutation_rate::Float64, cellbox)

    δ = rand(2) .- 0.5
    position_new_cell = cell[:position] + 2 .* δ/norm(δ)
    push!(cellbox, pos2box.(position_new_cell))

    new_b = cell[:birth_rate]
    new_d = cell[:death_rate]
    new_mutations = copy(cell[:mutations])

    push!(tumor, [index, position_new_cell, new_b, new_d, cell[:index], new_mutations])

    if rand() < mutation_rate
        push!(tumor.mutations[end], current_mut+1)
        return current_mut, position_new_cell
    end
end

function birth_death_pushing(; tumor_size::Int64, b::Float64, d::Float64, mu::Float64)
    cells = DataFrame([[1], [[0.0,0.0]], [b], [d], [0], [Int[]]], [:index, :position, :birth_rate, :death_rate, :parental, :mutations]) #DataFrame to save all relevant information about the single cells
    last_index = 1
    cur_mutation = 0

    mutation_events = DataFrame([Int[],Vector{Vector{Float64}}(undef,0)], [:mutation_number, :appearing_position])

    N = size(cells, 1)
    cellbox = [pos2box.(p) for p in cells.position]

    #Juno.progress() do id
    prog = ProgressUnknown("Progress: Tumor size ")
    while N < tumor_size
            N==0 && return "Tumor died"
            row_number_chosen_cell = rand(1:N) #chose randomly the row number (NOT index!) of the parental cell
            chosen_cell = cells[row_number_chosen_cell, :] #get the specific cell from the DataFrame

            b_prob = 1/(1 + chosen_cell.death_rate/chosen_cell.birth_rate) #probability from the rate
            if rand() <= b_prob
                last_index = last_index .+ 1 #add a new index
                N += 1
                mutation_event = birth_pushing!(cells, chosen_cell, last_index, cur_mutation, mu, cellbox)

                if !isnothing(mutation_event)
                    cur_mutation += 1
                    push!(mutation_events, mutation_event)
                end
                pushing!(cells, N, cellbox)
            else
                N -= 1
                delete!(cells, row_number_chosen_cell)
            end
        #print("Progress: $(size(cells,1))/$tumor_size | Recent mutation: $cur_mutation \r")
        #@info "Tumor size" progress=N/tumor_size _id=id
        ProgressMeter.update!(prog, N)
    end
    return cells, mutation_events
end
