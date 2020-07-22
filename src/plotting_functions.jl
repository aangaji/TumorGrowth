export plotting_normal, plotting_colored_mutations, plotting_normal_3D, plotting_pushing_colored_mutations_3D

#function to create circles in a 2D plot (here: cells with radius r)
function circle(pos_x, pos_y, r)
    theta = LinRange(0, 2*pi, 500)
    return pos_x .+ r*sin.(theta), pos_y .+ r*cos.(theta)
end

#function to create spheres in a 3D plot (here: cells with radius r)
function sphere(pos_x, pos_y, pos_z, r)
    theta = LinRange(0, 2*pi, 500)
    phi = LinRange(0, pi, 500)
    return pos_x .+ r*sin.(theta).*cos.(phi), pos_y .+ r*sin.(theta).*sin.(phi), pos_z .+ r*cos.(theta)
end

#function to plot a tumor with each mutation colored differently
function plotting_normal(data::DataFrame; path="", annotate=true)

    pos = data.position

    color_scheme = palette(:tab10)
    Plots.gr()
    p = Plots.plot()

    for (i,mut) in enumerate(data.mutations)
        color_choice = 1
        if ! isempty(mut)
            color_choice = mod(maximum(mut), 10)+1 #choose a different color for each mutation out of a colorset with 10 colors
            annotate && annotate!(pos[i]..., Plots.text(mut[1], 5, :black)) #annotate the number of the first mutation to each cell
        end
        Plots.plot!(p, circle(pos[i]..., r0), seriestype = [:shape,], lw = 0.5, c = color_scheme[color_choice], linecolor = :black, legend = false, fillalpha = 0.5, aspect_ratio = 1) #plot the cell with its specific color choice
        #Plots.scatter!(p, marker=(1, 0.5, color_scheme[color_choice]), legend = false)
    end

    Plots.plot!(p, circle(pos[1]..., r0), seriestype = [:shape,], lw = 0.5, c = :yellow, linecolor = :black, legend = false, fillalpha = 1.0, aspect_ratio = 1) #plot the start cell (lowest index) with yellow filling

    if !isempty(path) savefig("$(path).pdf") end
    return p
end

#function to plot the tumor similar to Ling et al. (2015), color code in experiments
function plotting_colored_mutations(data::DataFrame; path="")

    positions = data[!, :position]
    mutations = data[!, :mutations]

    #Plotting

    x = []          ### Same here
    y = []
    for i in positions
        push!(x, i[1])
        push!(y, i[2])
    end

    color_scheme = palette(:tab10)
    Plots.gr()
    p = Plots.plot()

    #initialize an array to distinguish colors of mutations
    colored_mutations = [[]]            ### NOT Any[Any[]]
    deleteat!(colored_mutations,1)      ### UNNECESSARY IF NOT Any[Any[]]

    for i in 2:length(x)
        color_found = false
        if ! isempty(mutations[i]) #if the cell has mutations
            if isempty(colored_mutations) #if this is the first mutated cell
                push!(colored_mutations, mutations[i]) #push the mutations in colored mutations
            else
                for j in 1:length(colored_mutations) #search for the specific color (color j) and plot the cell with this color (multiple times if the succesive mutations are already known)           ### (j,colored_mut) in enumerate(colored_mutations)
                    if mutations[i][1] in colored_mutations[j]
                        for k in 1:length(mutations[i])         ### mut in mutations[i] ?
                            Plots.plot!(p, circle(x[i], y[i], r0), seriestype = [:shape,], lw = 0.5, c = color_scheme[mod(j,10)+1], linecolor = :black, legend = false, fillalpha = 0.2, aspect_ratio = 1)
                            if !(mutations[i][k] in colored_mutations[j])
                                push!(colored_mutations[j], mutations[i][k])
                            end
                        end
                        color_found = true
                    end
                end
                if color_found == false #if the first mutation of a cell wasn't found already, choose a new color by pushing it into the srray colored_mutations        ### ! color_found
                    for k in 1:length(mutations[i])         ### WHAT DOES THIS LOOP DO?
                        Plots.plot!(p, circle(x[i], y[i], r0), seriestype = [:shape,], lw = 0.5, c = color_scheme[mod(length(colored_mutations)+1,10)+1], linecolor = :black, legend = false, fillalpha = 0.2, aspect_ratio = 1)
                    end
                    push!(colored_mutations, mutations[i])
                end
            end
        else
            Plots.plot!(p, circle(x[i], y[i], r0), seriestype = [:shape,], lw = 0.5, c = color_scheme[1], linecolor = :black, legend = false, fillalpha = 0.0, aspect_ratio = 1)
        end
    end

    Plots.plot!(p, circle(x[1], y[1], r0), seriestype = [:shape,], lw = 0.5, c = :yellow, linecolor = :black, legend = false, fillalpha = 1.0, aspect_ratio = 1)

    if !isempty(path) savefig("$(path)_colored.pdf") end
    return p
end



#function to plot a tumor with color determined by the latest mutation of the cell
function plotting_normal_3D(data::DataFrame; path="")

    p = [getindex.(data.position, i) for i=1:3]

    color_choice = ones(Int64, size(data,1))

    for (i,mut) in enumerate(data.mutations)
        if ! isempty(mut) color_choice[i] = mod(mut[end], 68)+1 end #choose a different color for each mutation out of a colorset with 35 colors
    end

    scene = meshscatter(p..., markersize = 1.0, color = color_choice, colormap=distinguishable_colors(68)) #plot the cell with its specific color choice

    if !isempty(path) save("$(path).jpg") end
    return scene

end


### THIS MIGHT BE BROKEN ??

#function to plot the tumor similar to Ling et al. (2015), color code in experiments
function plotting_pushing_colored_mutations_3D(pushing_data::DataFrame; path="")

    positions = pushing_data[!, :position]
    mutations = pushing_data[!, :mutations]

    x = Float64[r[1] for r in positions]
    y = Float64[r[2] for r in positions]
    z = Float64[r[3] for r in positions]

    color_scheme = palette(:tab10)

    color_choice = Array{Int64}(undef, length(x))
    alphas = zeros(Int64, length(x))

    #initialize an array to distinguish colors of mutations
    colored_mutations = [[]]
    deleteat!(colored_mutations,1)

    for i in 1:length(x)
        color_found = false
        if ! isempty(mutations[i]) #if the cell has mutations
            if isempty(colored_mutations) #if this is the first mutated cell
                push!(colored_mutations, mutations[i]) #push the mutations in colored mutations
            else
                for j in 1:length(colored_mutations) #search for the specific color (color j) and plot the cell with this color (multiple times if the succesive mutations are already known)
                    if mutations[i][1] in colored_mutations[j]
                        for k in 1:length(mutations[i])
                            alphas[i] += 1
                            color_choice[i] = mod(j,10)+1
                            if !(mutations[i][k] in colored_mutations[j])
                                push!(colored_mutations[j], mutations[i][k])
                            end
                        end
                        color_found = true
                    end
                end
                if color_found == false #if the first mutation of a cell wasn't found already, choose a new color by pushing it into the array colored_mutations
                    alphas[i] += length(mutations[i])
                    color_choice[i] = mod(length(colored_mutations)+1,10)+1
                    push!(colored_mutations, mutations[i])
                end
            end
        else
            color_choice[i] = 10
        end
    end

    alpha = alphas/maximum(alphas)

    scene = meshscatter(x,y,z, markersize = 0.5, color = color_choice, colormap=distinguishable_colors(10))

    if !isempty(path) save("$(path)_colored.jpg") end
    return scene
end
