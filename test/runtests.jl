using Test, MixedLayerThermoclineDynamics

@time @testset "Grid tests" begin
    include("test_grids.jl")

    nx, ny = 10, 12
    Lx, Ly = 2.0, 2.4
    
    grid1D = MixedLayerThermoclineDynamics.Grid1D(nx, 0, Lx)
    
    @test test_dx(grid1D, Lx/nx)
end