export birth_death_space_restricted


#function to give birth to a new cell in the tumor (called with the tumor, the parental cell and the index of mutations)
function birth_space_restricted(all_cells::DataFrame, cell::DataFrameRow, current_mut::Int64, mutation_rate::Float64, steps::Int64)

    #determine position

    pos_change = false
    position_new_cell = [0.0, 0.0]
    count = 0
    count_max = 1#round(Int64, 0.1*steps)

    while (pos_change == 0) & (count <= count_max)

        overlap = length(all_cells[!, :index]) #"negative" counter of the overlapping cells

        δ = 2 .* rand(2) .- 1 #normalized random vector
        δ = δ ./ norm(δ)

        pos = cell[:position]
        position_new_cell = pos + 2 .* δ #position of the new cell at the surface of the parental cell randomly

        c_j = 0 #Index of the cell in the next loop

        for j in all_cells[!, :position]

            c_j += 1

            energy = U(j, position_new_cell)
            if energy < 0
                overlap -= 1
            else
                if rand() < exp(-energy/(kb*T)) #Monte Carlo like behavior (also allow small overlap of two cells -> squeezing in nature)
                    overlap -= 1
                end
            end
        end

        if overlap == 0
            pos_change = true #if there is no overlap of the new cell with any other: Accept the position
        else
            count += 1 #else: try again with a new position around the same parental cell
        end
    end

    if count == count_max+1         ### JUST !pos_change ?
        return nothing, steps+1 #If the positioning was tried several times (count_max) and it didn't work: Abort!      ### WHY steps+1 ?
    end

    #determine new properties (mutation, birth and death rate etc)

    new_b = cell[:birth_rate]
    new_d = cell[:death_rate]
    new_mutations = copy(cell[:mutations])
    mutation_event = [] #Array to save the mutation event positions -> Otherwise probably effects of algorithm on analyzed data!

    if rand() < mutation_rate
        current_mut += 1

        mutation_event = [current_mut, position_new_cell]

        push!(new_mutations, current_mut)
    end

    return [position_new_cell, new_b, new_d, cell[:index], new_mutations], current_mut, mutation_event, count
end


function birth_death_space_restricted(; tumor_size::Int, b::Float64, d::Float64, mu::Float64)
    cells = DataFrame([[1], [[0.0,0.0]], [b], [d], [0], [Int[]]], [:index, :position, :birth_rate, :death_rate, :parental, :mutations])
    last_index = 0
    cur_mutation = 0
    count_distribution = []

    mutation_events = DataFrame([[],[]], [:mutation_number, :appearing_position]) #Array to save the mutation events (number and position at that time)

    while length(cells[:,1]) < tumor_size

        isempty(cells) && return "Tumor died"

        row_number_chosen_cell = rand(1:size(cells, 1))
        chosen_cell = cells[row_number_chosen_cell, :]
        b_rate = chosen_cell[:birth_rate]
        d_rate = chosen_cell[:death_rate]
        b_prob = b_rate/(b_rate+d_rate)

        r = rand()
        if r <= b_prob #if birth is chosen
            last_index = last_index .+ 1
            new_born_cell = birth_space_restricted(cells, chosen_cell, cur_mutation, mu, tumor_size) #give birth to a new cell
            if new_born_cell[1] != nothing #if new cell could be positioned
                cur_mutation = new_born_cell[2]
                push!(cells, vcat(last_index .+ 1, new_born_cell[1])) #push in the cells the index and the informations from the new cell

                mutation_event = new_born_cell[3]
                push!(count_distribution, new_born_cell[4])
                if mutation_event != [] #If there was a new mutation event, push it into the array
                    push!(mutation_events, mutation_event)
                end

            else
                push!(count_distribution, new_born_cell[2])
            end
        else
            delete!(cells, row_number_chosen_cell)
        end

        print("Progress: $(size(cells,1))/$tumor_size | Recent mutation: $cur_mutation \r")
    end

    return cells, mutation_events, count_distribution

end
