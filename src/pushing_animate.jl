function update_plot!(scene, tumor, root, n, queue; color, sleeptime)
    p = [Point(cell.position...) for cell in tumor]
    while length(scene) > 1
      delete!(scene, scene[end])
    end
    meshscatter!(scene, p, markersize = 1.0, color = :lightgrey, scale_plot = false, shading=false)
    isempty(queue) || meshscatter!(scene, p[queue], markersize = 1.0, color = :grey, scale_plot = false, shading=false)
    isempty(root) || meshscatter!(scene, p[root], markersize = 1.0, color = :blue, scale_plot = false, shading=false)
    isempty(n) || meshscatter!(scene, p[n], markersize = 1.0, color = color, scale_plot = false, shading=false)
    scene[Axis][:ticks][:ranges] = (-12:4:12, -12:4:12)
    scene[Axis][:grid][:linewidth] = (3,3)
    #sleep(sleeptime)
end


function pushing!_animate(tumor, root, cellbox)
    queue = [root]
    N = length(tumor)

    scene= Scene()
    p = [Point(cell.position...) for cell in tumor]
    meshscatter!(scene, p, markersize = 1.0, color = :lightgrey, scale_plot = false, shading=false)

    record(scene, "test.mp4") do io
        while !isempty(queue)
            root = popfirst!(queue)
            r1 = tumor[root].position

            for n in shuffle!( (1:N)[find_neighbors(cellbox, root)] )
                r2 = tumor[n].position
                d = r2.-r1

                update_plot!(scene, tumor, root, n, queue; color=:orange, sleeptime=0.02)
                recordframe!(io)

                if norm(d) < 2.

                    update_plot!(scene, tumor, root, n, queue; color=:red,sleeptime=0.5)
                    recordframe!(io)

                    r2 .+= d*(2.2/norm(d) - 1)#/2
                    #r1 .-= d*(2.2/norm(d) - 1)/2
                    cellbox[n] = pos2box.(r2)
                    push!(queue, n)

                    update_plot!(scene, tumor, root, n, queue; color=:red,sleeptime=0.5)
                    recordframe!(io)
                end
            end
        end
        update_plot!(scene, tumor, root, Int[], queue; color=:red,sleeptime=0.5)
        recordframe!(io)
    end
end

function pushing_recursive!_animate(tumor, root, cellbox)
    scene= Scene()
    p = [Point(cell.position...) for cell in tumor]
    meshscatter!(scene, p, markersize = 1.0, color = :lightgrey, scale_plot = false, shading=false)
    queue = [root]
    record(scene, "test.mp4") do io
        pushing_recursive!_animate_exec(scene, io, tumor, root, queue, cellbox)
        update_plot!(scene, tumor, root, Int[], Int[]; color=:red,sleeptime=0.5)
        recordframe!(io)
    end
end

function pushing_recursive!_animate_exec(scene, io, tumor, root, queue, cellbox)
    r1 = tumor[root].position

    N=length(tumor)
    for n in shuffle!( (1:N)[find_neighbors(cellbox, root)] )

        update_plot!(scene, tumor, root, n, queue; color=:orange,sleeptime=0.02)
        recordframe!(io)

        r2 = tumor[n].position
        d = r2.-r1
        if norm(d) < 2.

            update_plot!(scene, tumor, root, n, queue; color=:red,sleeptime=0.5)
            recordframe!(io)

            push!(queue, n)
            r2 .+= d*(2.2/norm(d) - 1)#/2
            #r1 .-= d*(2.2/norm(d) - 1)/2
            cellbox[n] = pos2box.(r2)

            update_plot!(scene, tumor, root, n, queue; color=:red,sleeptime=0.5)
            recordframe!(io)

            pushing_recursive!_animate_exec(scene, io, tumor, n, queue, cellbox)
            pop!(queue)
        end
    end
end
