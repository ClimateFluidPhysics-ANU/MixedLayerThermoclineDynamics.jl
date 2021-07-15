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

@time @testset "Field tests" begin
    include("test_fields.jl")

    nx, ny = 10, 12
    Lx, Ly = 2.0, 2.4
    
    dx, dy = Lx/nx, Ly/ny
    
    grid1D = Grid1D(nx, 0, Lx)
    hdata = zeros(nx)
    udata = zeros(nx)
    hdata = sin.(2*π*grid1D.xC/Lx)
    udata = sin.(2*π*grid1D.xF/Lx)

    h1D = Field(Centre, hdata, grid1D)
    u1D = Field(Face, udata, grid1D)
    @test location_centre(h1D)
    @test location_face(u1D)

    @test test_grid(h1D, grid1D)
    @test test_grid(u1D, grid1D)

    @test test_data(h1D, hdata)
    @test test_data(u1D, udata)
end


