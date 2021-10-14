using Pkg; Pkg.activate(pwd()); Pkg.instantiate()
using Revise

@time using TumorGrowth

c = TumorGrowth.HyperCube(zeros(3), 10)
lattice = TumorGrowth.triangular_lattice(c; n=11)

sc = TumorGrowth.GLMakie.meshscatter(
    reshape(TumorGrowth.Point{3}.(lattice), prod(size(lattice))), markersize=0.1,
    )
display(sc)

lattice

lattice[1]
sc = plotting( DataFrame(position = reshape(lattice, prod(size(lattice)))),
    color = :black, markersize=0.1)
TumorGrowth.GLMakie.zlims!(sc, (0,4))

lattice
