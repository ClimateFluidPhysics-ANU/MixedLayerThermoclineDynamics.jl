using Test, MixedLayerThermoclineDynamics

@time @testset "Grid tests" begin
    include("test_grids.jl")

    nx, ny = 10, 12
    Lx, Ly = 2.0, 2.4
    
    grid1D = Grid1D(nx, 0, Lx)

    @test test_dx(grid1D, Lx/nx)
    @test test_xF(grid1D, range(0, stop = Lx - grid1D.dx, length = nx))
    @test xdomain_length(grid1D, Lx - 0)

    grid2D = Grid2D(nx, ny, 0, Lx, 0, Ly)

    @test test_dx(grid2D, Lx/nx)
    @test test_dy(grid2D, Ly/ny)
    @test test_xC(grid2D, range(grid2D.dx/2, stop = Lx - grid2D.dx/2, length = nx))
    @test test_xF(grid2D, range(0, stop = Lx - grid2D.dx, length = nx))
    @test test_yF(grid2D, range(0, stop = Ly - grid2D.dy, length = ny))
    @test test_yC(grid2D, range(grid2D.dy/2, stop = Ly - grid2D.dy/2, length = ny))
    @test xdomain_length(grid2D, Lx - 0)
    @test ydomain_length(grid2D, Ly - 0)
end
