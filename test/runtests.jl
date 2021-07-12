using Test, MixedLayerThermoclineDynamics

@time @testset "Grid tests" begin
    include("test_grids.jl")

    nx, ny = 10, 12
    Lx, Ly = 2.0, 2.4
    
    grid1D = Grid1D(nx, 0, Lx)

    @test test_dx(grid1D, Lx/nx)
    @test test_xF(grid1D, LinRange(0, Lx, nx))
    @test xdomain_length(grid1D, Lx - 0)

    grid2D = Grid2D(nx, ny, 0, Lx, 0, Ly)

    @test test_dx(grid2D, Lx/nx)
    @test test_dy(grid2D, Ly/ny)
    @test test_xF(grid2D, LinRange(grid2D.dx/2, Lx - grid2D.dx/2, nx))
    @test test_xC(grid2D, LinRange(0, Lx - grid2D.dx, nx))
    @test test_xF(grid2D, LinRange(0, Ly - grid2D.dy, ny))
    @test test_yC(grid2D, LinRange(grid2D.dy/2, Ly - grid2D.dy/2, ny))
    @test xdomain_length(grid2D, Lx - 0)
    @test ydomain_length(grid2D, Ly - 0)
end
