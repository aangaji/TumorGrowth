function pushing!_animate(tumor, root, cellbox, dimv::Val)
    queue = [root]
    N = length(tumor)

    points = [Point2f0(cell.position...) for cell in tumor]
    pointsNode = Node(points)
    colors = fill(:lightgrey, length(tumor))
    colorsNode = Node(colors)

    scene = meshscatter(pointsNode, markersize = 1.0, color = colorsNode, scale_plot = false, shading=false)
    xlims!(-30,30)
    ylims!(-30,30)
    scene[Axis][:ticks][:ranges] = (-20:4:20, -20:4:20)
    display(scene)

    #record(scene, "test.mp4") do io
        while !isempty(queue)
            root = popfirst!(queue)
            r1 = tumor[root].position

            setindex!(colors, :blue, root)

            for n in shuffle!( findall(find_neighbors(cellbox, root; s=2)) )
                r2 = tumor[n].position
                d = r2.-r1

                setindex!(colors, :orange, n)
                colorsNode[] = colors
                sleep(0.01)
                # recordframe!(io)

                if norm(d) < 2.

                    setindex!(colors, :red, n)
                    colorsNode[] = colors
                    sleep(0.01)
                    # recordframe!(io)

                    r2 .+= d*(2.2/norm(d) - 1)#/2
                    #r1 .-= d*(2.2/norm(d) - 1)/2
                    cellbox[n] = pos2box(r2, dimv)
                    push!(queue, n)

                    setindex!(points, Point2f0(r2...), n)
                    pointsNode[] = points
                    sleep(0.01)
                    setindex!(colors, :grey, n)
                    # recordframe!(io)
                else
                    setindex!(colors, :lightgrey, n)
                end
            end
            setindex!(colors, :lightgrey, root)
        end
        # recordframe!(io)
    #end
end

function pushing_recursive!_animate(tumor, root, cellbox, dimv::Val)
    points = [Point2f0(cell.position...) for cell in tumor]
    pointsNode = Node(points)
    colors = fill(:lightgrey, length(tumor))
    colorsNode = Node(colors)

    scene = meshscatter(pointsNode, markersize = 1.0, color = colorsNode, scale_plot = false, shading=false)
    xlims!(-30,30)
    ylims!(-30,30)
    scene[Axis][:ticks][:ranges] = (-20:4:20, -20:4:20)
    display(scene)
    queue = Int[]
    record(scene, "test_rec.mp4") do io
        pushing_recursive!_animate_exec(points, pointsNode, colors, colorsNode, io, tumor, root, cellbox, dimv, queue)
        recordframe!(io)
    end
end

function pushing_recursive!_animate_exec(points, pointsNode, colors, colorsNode, io, tumor, root, cellbox, dimv, queue)
    r1 = tumor[root].position

    colors[root] = :blue
    push!(queue, root)

    N=length(tumor)
    for n in shuffle!( findall(find_neighbors(cellbox, root; s=2)) )

        setindex!(colors, :orange, n)
        colorsNode[] = colors
        sleep(0.01)
        recordframe!(io)

        r2 = tumor[n].position
        d = r2.-r1
        if norm(d) < 2.

            setindex!(colors, :red, n)
            colorsNode[] = colors
            sleep(0.01)
            recordframe!(io)

            r2 .+= d*(2.2/norm(d) - 1)#/2
            #r1 .-= d*(2.2/norm(d) - 1)/2
            cellbox[n] = pos2box(r2, dimv)

            setindex!(points, Point2f0(r2...), n)
            pointsNode[] = points
            sleep(0.01)
            setindex!(colors, :grey, root)
            recordframe!(io)

            pushing_recursive!_animate_exec(points, pointsNode, colors, colorsNode, io, tumor, n, cellbox, dimv, queue)
        else
            colors[n] = n in queue ? :grey : :lightgrey
        end
    end
    pop!(queue)
    colors[root] = :lightgrey
end
