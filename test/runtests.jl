using MixedLayerThermoclineDynamics, OffsetArrays, Test

@time @testset "Grid tests" begin
    include("test_grids.jl")
    
    for Tx in (Periodic(), Bounded()), Ty in (Periodic(), Bounded())
        
        nx, ny = 10, 12
        hx, hy =  1, 2
        x0, y0 = 0.2, -0.4
        Lx, Ly = 2.0, 2.4
        
        dx, dy = Lx/nx, Ly/ny
        
        grid1D = Grid1D(Tx, nx, x0, x0 + Lx; hx=hx)        
        grid2D = Grid2D(Tx, Ty, nx, ny, x0, x0 + Lx, y0, y0 + Ly; hx=hx, hy=hy)

            
        xF = isa(Tx, Periodic) ? range(x0 - hx*dx, stop = x0 + Lx - dx + hx*dx, length = nx + 2hx) :
                                 range(x0 - hx*dx, stop = x0 + Lx + hx*dx, length = nx + 1 + 2hx)

        xC = range(x0 + dx/2 - hx*dx, stop = x0 + Lx - dx/2 + hx*dx, length = nx + 2hx)
        
        yF = isa(Ty, Periodic) ? range(y0 - hy*dy, stop = y0 + Ly - dy + hy*dy, length = ny + 2hy) :
                                 range(y0 - hy*dy, stop = y0 + Ly + hy*dy, length = ny + 1 + 2hy)
        
        yC = range(y0 + dy/2 - hy*dy, stop = y0 + Ly - dy/2 + hy*dy, length = ny + 2hy)
        
        xF = OffsetArray(xF, -hx)
        xC = OffsetArray(xC, -hx)
        yF = OffsetArray(yF, -hy)
        yC = OffsetArray(yC, -hy)
        
        @test test_dx(grid1D, dx)
        @test test_xF(grid1D, xF)
        @test test_xC(grid1D, xC)
        @test xdomain_length(grid1D, Lx)

        @test test_dx(grid2D, dx)
        @test test_dy(grid2D, dy)
        @test test_xF(grid2D, xF)
        @test test_xC(grid2D, xC)
        @test test_yF(grid2D, yF)
        @test test_yC(grid2D, yC)
        @test xdomain_length(grid2D, Lx)
        @test ydomain_length(grid2D, Ly)
    end
end

@time @testset "Field tests" begin
    nx, ny = 10, 12
    Lx, Ly = 2.0, 2.4
    
    dx, dy = Lx/nx, Ly/ny
    
    grid1D = Grid1D(Periodic(), nx, 0, Lx)
    grid2D = Grid2D(Periodic(), Periodic(), nx, ny, 0, Lx, 0, Ly)
    
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
end
