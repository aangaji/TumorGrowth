export plotting_2d, plotting_2d!, plotting_2d_colored_mutations, plotting, plotting!, plotting_colored_mutations

#function to create circles in a 2D plot (here: cells with radius r)
function circle(pos_x, pos_y, r)
    theta = LinRange(0, 2*pi, 20)
    return pos_x .+ r*sin.(theta), pos_y .+ r*cos.(theta)
end

#function to plot a tumor with each mutation colored differently
function plotting_2d!(p, data::DataFrame; path="", annotate=false, color=nothing)

    pos = data.position
    color_scheme = palette(:tab10)
    choose_color = false
    isnothing(color) && (choose_color=true)
    for (i,mut) in enumerate(data.mutations)
        color_choice = 1
        if ! isempty(mut)
            color_choice = mod(maximum(mut), 10)+1 #choose a different color for each mutation out of a colorset with 10 colors
            annotate && annotate!(pos[i]..., text(mut[1], 5, :black)) #annotate the number of the first mutation to each cell
        end
        choose_color && (color = color_scheme[color_choice])
        plot!(p, circle(pos[i]..., r0), seriestype = [:shape,], lw = 0.5, c = color, linecolor = :black, lab = :none, fillalpha = 0.5, aspect_ratio = 1) #plot the cell with its specific color choice
    end
    first = findfirst(isequal(1), data.index)
    if !isnothing(first)
        plot!(p, circle(pos[first]..., r0), seriestype = [:shape,], lw = 0.5, c = :black, linecolor = :black, lab = "first", fillalpha = 1.0, aspect_ratio = 1) #plot the start cell (lowest index) with black filling
    end
    if !isempty(path) savefig("$(path).pdf") end
    return p
end

plotting_2d(data::DataFrame; path="", annotate=false) = plotting_2d!(plot(), data; path=path, annotate=annotate)

#function to plot the tumor similar to Ling et al. (2015), color code in experiments
function plotting_2d_colored_mutations(data::DataFrame; path="")

    pos = data.position

    color_scheme = palette(:tab10)
    p = plot()

    #initialize an array to distinguish colors of mutations
    colored_mutations = Vector{Vector{Int64}}(undef,0)

    for (i, mut_list) in enumerate(data.mutations)
        color_found = false
        if ! isempty(mut_list) #if the cell has mutations
            if isempty(colored_mutations) #if this is the first mutated cell
                push!(colored_mutations, mut_list) #push the mutations in colored mutations
            else
                for (j,c_mut) in enumerate(colored_mutations) #search for the specific color (color j) and plot the cell with this color (multiple times if the succesive mutations are already known)
                    if mut_list[1] in c_mut
                        for mut in mut_list
                            plot!(p, circle(pos[i]..., r0), seriestype = [:shape,], lw = 0.5, c = color_scheme[mod(j,10)+1], linecolor = :black, lab = :none, fillalpha = 0.2, aspect_ratio = 1)
                            if !(mut in c_mut)
                                push!(colored_mutations[j], mut)
                            end
                        end
                        color_found = true
                    end
                end
                if !color_found #if the first mutation of a cell wasn't found already, choose a new color by pushing it into the srray colored_mutations
                    for mut in mut_list       ### WHAT DOES THIS LOOP DO?
                        plot!(p, circle(pos[i]..., r0), seriestype = [:shape,], lw = 0.5, c = color_scheme[mod(length(colored_mutations)+1,10)+1], linecolor = :black, lab=:none, fillalpha = 0.2, aspect_ratio = 1)
                    end
                    push!(colored_mutations, mut_list)
                end
            end
        else
            plot!(p, circle(pos[i]..., r0), seriestype = [:shape,], lw = 0.5, c = color_scheme[1], linecolor = :black, lab=:none, fillalpha = 0.0, aspect_ratio = 1)
        end
    end

    first = findfirst(isequal(1), data.index)
    if !isnothing(first)
        plot!(p, circle(pos[first]..., r0), seriestype = [:shape,], lw = 0.5, c = :black, linecolor = :black, lab = "first", fillalpha = 1.0, aspect_ratio = 1)
    end
    if !isempty(path) savefig("$(path)_colored.pdf") end
    return p
end


#function to plot a tumor with color determined by the latest mutation of the cell
function plotting!(scene, data::DataFrame; path="", color=nothing)

    p = (getindex.(data.position, i) for i=1:length(data.position[1]))

    if isnothing(color)
        color = ones(Int64, size(data,1))

        for (i,mut) in enumerate(data.mutations)
            if ! isempty(mut) color[i] = mod(mut[end], 68)+1 end #choose a different color for each mutation out of a colorset with 35 colors
        end
    end
    meshscatter!(scene, p..., markersize = 1.0, color = color, colormap=distinguishable_colors(68), scale_plot = false)
    if !isempty(path) save("$(path).jpg") end
    return scene
end
plotting(data::DataFrame; path="", color=nothing) = plotting!(Scene(), data::DataFrame; path=path, color=color)


#function to plot the tumor similar to Ling et al. (2015), color code in experiments
function plotting_colored_mutations(data::DataFrame; path="")

    p = ( getindex.(data.position, i) for i=1:length(data.position[1]) )

    color_scheme = palette(:tab10)

    color_choice = Vector{Int64}(undef, size(data,1))

    colored_mutations = Vector{Vector{Int64}}(undef,0)

    for (i, mut_list) in enumerate(data.mutations)
        color_found = false
        if ! isempty(mut_list) #if the cell has mutations
            if isempty(colored_mutations) #if this is the first mutated cell
                push!(colored_mutations, mut_list) #push the mutations in colored mutations
            else
                for (j, cmut) in enumerate(colored_mutations) #search for the specific color (color j) and plot the cell with this color (multiple times if the succesive mutations are already known)
                    if mut_list[1] in cmut
                        for mut in mut_list
                            color_choice[i] = mod(j,10)+1
                            if !(mut in cmut)
                                push!(colored_mutations[j], mut)
                            end
                        end
                        color_found = true
                    end
                end
                if !color_found #if the first mutation of a cell wasn't found already, choose a new color by pushing it into the array colored_mutations
                    color_choice[i] = mod(length(colored_mutations)+1,10)+1
                    push!(colored_mutations, mut_list)
                end
            end
        else
            color_choice[i] = 10
        end
    end

    scene = meshscatter(p..., markersize = 1., color = color_choice, colormap=distinguishable_colors(10), scale_plot = false)

    if !isempty(path) save("$(path)_colored.jpg") end
    return scene
end
