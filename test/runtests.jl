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
    grid2D = Grid2D(nx, ny, 0, Lx, 0, Ly)
    
    # 1D Fields
    hdata = @. sin(2π * grid1D.xC / Lx)
    udata = @. cos(2π * grid1D.xF / Lx)

    h1D = Field(Centre, hdata, grid1D)
    u1D = Field(Face, udata, grid1D)
    
    @test typeof(h1D) <: Field1D{Centre}
    @test typeof(u1D) <: Field1D{Face}

    @test h1D.grid == grid1D
    @test u1D.grid == grid1D

    @test h1D.data == hdata
    @test u1D.data == udata
    
    # 2D Fields
    hdata = [sin(2π * grid2D.xC[i]) * cos(4π * grid2D.yC[j]) for i in 1:nx, j in 1:ny]
    udata = [cos(6π * grid2D.xF[i]) * cos(2π * grid2D.yC[j]) for i in 1:nx, j in 1:ny]
    vdata = [sin(8π * grid2D.xC[i]) * sin(6π * grid2D.yF[j]) for i in 1:nx, j in 1:ny]

    h2D = Field(Centre, Center, hdata, grid2D)
    u2D = Field(Face, Center, udata, grid2D)
    v2D = Field(Centre, Face, vdata, grid2D)
    
    @test typeof(h2D) <: Field2D{Centre, Centre}
    @test typeof(u2D) <: Field2D{Face, Centre}
    @test typeof(v2D) <: Field2D{Centre, Face}

    @test h2D.grid == grid2D
    @test u2D.grid == grid2D
    @test v2D.grid == grid2D

    @test h2D.data == hdata
    @test u2D.data == udata
    @test v2D.data == vdata
    
    for field in [h1D, u1D, h2D, u2D, v2D]
        @test typeof(field) <: AbstractField
    end
    
    # test for interpolate
    rtol_interpolate = 1e-2
    nx, Lx = 100, 2.0
    dx = Lx / nx
    
    func(x) = sin(4π * x / Lx)
    
    grid1D = Grid1D(nx, 0, Lx)

    @test test_interpolation(func, grid1D, rtol_interpolate)
end
