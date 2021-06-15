export cross_section, radial_sample, punch, bulk, multi_region_sampling

function cross_section(tumor; x=nothing, y=nothing, z=nothing, width = 2., reduce_dim=false)
    for (d,coord) in enumerate((x,y,z))
        if !isnothing(coord)
            rows = -width .< getindex.(tumor.position,d).-coord .< width
            if reduce_dim
                plane = tumor[rows, :]
                plane.position = plane.position .|> p->reduce_dimension(p, d)
                return plane
            else
                return tumor[rows, :]
            end
        end
    end
end

# this function returns a copy without entries at dimension d
function reduce_dimension(position, d)
    dim = position |> length
    select = mod.((d-1:2:d+1).-1, dim).+1 |> unique
    return getindex(position, select)
end

radial_sample(tumor; r=2., width=2.) = tumor[ tumor.position .|> p-> abs(norm(p)-r)<width, :]

punch(tumor; pos, r=10.) = tumor[ tumor.position .|> p-> norm(p.-pos) < r, :]

bulk(tumor; pos, box) = tumor[ tumor.position .|> p-> all(abs.(p .- pos) .< box./2), :]

struct Square
    anchorpoint
    width
end

Base.in(point, c::Circle) = norm(point .- c.center) < c.r

function triangular_lattice(sq::Square; n=0, a=sqrt(π/n)*radius)
    v, w = [1.0, 0.0], [0.5, -cosd(30)]

    steps_x = ceil(Int, sq.width/a)
    steps_y = ceil(Int, sq.width/(a*cosd(30)))

    lattice = Matrix{Vector{Float64}}(undef, steps_x, steps_y)
    lattice[:, 1] .= [sq.anchorpoint .+ a*i*v for i=1:steps_x]
    for j = 2:steps_y
        @. lattice[:, j] = [p + a*(w - iseven(j)*v) for p in lattice[:, j-1]]
    end
    return lattice
end

function multi_region_sampling(tumor; n=0, a=0., cells_per_sample=0, sample_r=a/2)

    cm = mean(tumor.position)
    r = mean(norm.(p .- cm for p in tumor.position))
    density = size(punch(tumor; pos=cm, r=r), 1)/(π*r ^2)
    R = sqrt(size(tumor,1)/(π*density))

    sample_r = iszero(cells_per_sample) ? a/2 : sqrt(cells_per_sample/(π*density))
    if !iszero(n)
        a = sqrt(π/n)*R
        sample_r = iszero(sample_r) ? a/2 : sample_r
    elseif iszero(a)
        a = 2*sample_r
    end
    iszero(sample_r) && error("specify parameters")
    sample_r > a/2 && error("samples too large for given spacing")

    square = Square(cm .- R.*(1,-1), 2*R)
    lattice = [triangular_lattice(square; n=n, a=a)...]
    filter!(p -> p in Circle(Point{2}(cm), R - a/2), lattice)

    samples = lattice .|> p -> punch(tumor; pos=p, r=sample_r)

    return lattice, samples
end
