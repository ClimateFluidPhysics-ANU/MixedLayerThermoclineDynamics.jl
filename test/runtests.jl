using Test, MixedLayerThermoclineDynamics

@time @testset "Grid tests" begin
    include("test_grids.jl")

    nx, ny = 10, 12
    Lx, Ly = 2.0, 2.4
    
    dx, dy = Lx/nx, Ly/ny
    
    grid1D = Grid1D(nx, 0, Lx)

    @test test_dx(grid1D, dx)
    @test test_xF(grid1D, range(0, stop = Lx - dx, length = nx))
    @test xdomain_length(grid1D, Lx)

    grid2D = Grid2D(nx, ny, 0, Lx, 0, Ly)

    @test test_dx(grid2D, dx)
    @test test_dy(grid2D, dy)
    @test test_xC(grid2D, range(dx/2, stop = Lx - dx/2, length = nx))
    @test test_xF(grid2D, range(0, stop = Lx - dx, length = nx))
    @test test_yF(grid2D, range(0, stop = Ly - dy, length = ny))
    @test test_yC(grid2D, range(dy/2, stop = Ly - dy/2, length = ny))
    @test xdomain_length(grid2D, Lx)
    @test ydomain_length(grid2D, Ly)
end
