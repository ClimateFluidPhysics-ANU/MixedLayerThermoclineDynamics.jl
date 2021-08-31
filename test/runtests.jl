using MixedLayerThermoclineDynamics, OffsetArrays, Test

const rtol_interpolation = 1e-3
const rtol_derivatives = 1e-3

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
        
        @test grid1D.dx == dx
        @test test_xF(grid1D, xF)
        @test test_xC(grid1D, xC)
        @test xdomain_length(grid1D, Lx)

        @test grid2D.dx == dx
        @test grid2D.dy == dy
        @test test_xF(grid2D, xF)
        @test test_xC(grid2D, xC)
        @test test_yF(grid2D, yF)
        @test test_yC(grid2D, yC)
        @test xdomain_length(grid2D, Lx)
        @test ydomain_length(grid2D, Ly)
    end

end

@time @testset "Field tests" begin

    include("test_fields.jl")
    
    nx, ny = 100, 120
    Lx, Ly = 1.0, 1.2
    hx, hy = 2, 3

    @test_throws ErrorException("Number of halo points in x cannot be zero") Grid1D(Periodic(), nx, 0, Lx; hx = 0)
    @test_throws ErrorException("Number of halo points in x cannot be zero") Grid2D(Periodic(), Periodic(), nx, 0, Lx, ny, 0, Ly; hx = 0, hy = 1)
    @test_throws ErrorException("Number of halo points in y cannot be zero") Grid2D(Periodic(), Periodic(), nx, 0, Lx, ny, 0, Ly; hx = 1, hy = 0)

    dx, dy = Lx/nx, Ly/ny # note: this is only correct for Periodic
    
    grid1D = Grid1D(Periodic(), nx, 0, Lx; hx = hx)
    grid2D = Grid2D(Periodic(), Periodic(), nx, ny, 0, Lx, 0, Ly; hx = hx, hy = hy)
    
    # 1D Fields
    Cdata = @. sin(2π * grid1D.xC / Lx)
    Fdata = @. cos(2π * grid1D.xF / Lx)
    
    Cdata_on_F = @. sin(2π * grid1D.xF / Lx)
    Fdata_on_C = @. cos(2π * grid1D.xC / Lx)

    ∂x_Cdata = @.   2π/Lx * cos(2π * grid1D.xF / Lx)
    ∂x_Fdata = @. - 2π/Lx * sin(2π * grid1D.xC / Lx)

    ∂x_Cdata_on_C = @.   2π/Lx * cos(2π * grid1D.xC / Lx)
    ∂x_Fdata_on_F = @. - 2π/Lx * sin(2π * grid1D.xF / Lx)

    C_Field1D = Field(Centre, Cdata, grid1D)
    F_Field1D = Field(Face, Fdata, grid1D)

    Cdata_on_F_Field1D = Field(Face, Cdata_on_F, grid1D)
    Fdata_on_C_Field1D = Field(Centre, Fdata_on_C, grid1D)
    Ftest_Field1D = Field(Face, zero(Cdata), grid1D)
    Ctest_Field1D = Field(Centre, zero(Fdata), grid1D)

    ∂x_Cdata_Field1D = Field(Face, ∂x_Cdata, grid1D)
    ∂x_Fdata_Field1D = Field(Centre, ∂x_Fdata, grid1D)
    ∂x_Cdata_on_C_Field1D = Field(Centre, ∂x_Cdata_on_C, grid1D)
    ∂x_Fdata_on_F_Field1D = Field(Face, ∂x_Fdata_on_F, grid1D)
    
    C_Field1D_from_array = Field1D(Centre, Cdata[1:nx], grid1D)
    F_Field1D_from_array = Field1D(Face, Fdata[1:nx], grid1D)
    
    @test C_Field1D.grid == C_Field1D_from_array.grid
    @test F_Field1D.grid == F_Field1D_from_array.grid
    @test C_Field1D.data ≈ C_Field1D_from_array.data
    @test F_Field1D.data ≈ F_Field1D_from_array.data
    
    @test typeof(C_Field1D) <: Field1D{Centre}
    @test typeof(F_Field1D) <: Field1D{Face}
    @test typeof(C_Field1D_from_array) <: Field1D{Centre}
    @test typeof(F_Field1D_from_array) <: Field1D{Face}

    @test C_Field1D.grid == grid1D
    @test F_Field1D.grid == grid1D

    @test test_𝐼x(Fdata_on_C_Field1D, Ctest_Field1D, F_Field1D)
    @test test_𝐼x(Cdata_on_F_Field1D, Ftest_Field1D, C_Field1D)    
    
    @test test_∂x(∂x_Fdata_Field1D, Ctest_Field1D, F_Field1D)
    @test test_∂x(∂x_Cdata_Field1D, Ftest_Field1D, C_Field1D)
    @test test_∂x(∂x_Fdata_on_F_Field1D, Ftest_Field1D, F_Field1D)
    @test test_∂x(∂x_Cdata_on_C_Field1D, Ctest_Field1D, C_Field1D)

    Cdata_with_halos = OffsetArray(zeros(nx + 2hx), -hx)

    for i in 1:nx
        Cdata_with_halos[i] = Cdata[i]
    end

    for i in 1:hx
        Cdata_with_halos[nx+i] = Cdata[i]
        Cdata_with_halos[-i+1] = Cdata[nx-i+1]
    end

    Fdata_with_halos = OffsetArray(zeros(nx + 2hx), -hx)

    for i in 1:nx
        Fdata_with_halos[i] = Fdata[i]
    end

    for i in 1:hx
        Fdata_with_halos[nx+i] = Fdata[i]
        Fdata_with_halos[-i+1] = Fdata[nx-i+1]
    end
    
    @test C_Field1D.data ≈ Cdata_with_halos
    @test F_Field1D.data ≈ Fdata_with_halos
    
    #=
    # 2D Fields
    Cdata = @. [sin(2π * grid2D.xC[i]/Lx) * cos(4π * grid2D.yC[j]/Ly) for i in 1:nx, j in 1:ny]
    udata = @. [cos(6π * grid2D.xF[i]/Lx) * cos(2π * grid2D.yC[j]/Ly) for i in 1:nx, j in 1:ny]
    vdata = @. [sin(8π * grid2D.xC[i]/Lx) * sin(6π * grid2D.yF[j]/Ly) for i in 1:nx, j in 1:ny]

    𝐼hudata = @. [sin(2π * grid2D.xF[i]/Lx) * cos(4π * grid2D.yC[j]/Ly) for i in 1:nx, j in 1:ny]
    𝐼hvdata = @. [sin(2π * grid2D.xC[i]/Lx) * cos(4π * grid2D.yF[j]/Ly) for i in 1:nx, j in 1:ny]
    𝐼udata = @. [cos(6π * grid2D.xC[i]/Lx) * cos(2π * grid2D.yC[j]/Ly) for i in 1:nx, j in 1:ny]
    𝐼vdata = @. [sin(8π * grid2D.xC[i]/Lx) * sin(6π * grid2D.yC[j]/Ly) for i in 1:nx, j in 1:ny]

    ∂hudata = @. (2π/Lx) * [cos(2π * grid2D.xF[i]/Lx) * cos(4π * grid2D.yC[j]/Ly) for i in 1:nx, j in 1:ny]
    ∂hvdata = @. (-4π/Ly) * [sin(2π * grid2D.xC[i]/Lx) * sin(4π * grid2D.yF[j]/Ly) for i in 1:nx, j in 1:ny]
    ∂udata = @. (-6π/Lx) * [sin(6π * grid2D.xC[i]/Lx) * cos(2π * grid2D.yC[j]/Ly) for i in 1:nx, j in 1:ny]
    ∂vdata = @. (6π/Ly) * [sin(8π * grid2D.xC[i]/Lx) * cos(6π * grid2D.yC[j]/Ly) for i in 1:nx, j in 1:ny]

    ∂hx_on_h_data = @. (2π/Lx) * [cos(2π * grid2D.xC[i]/Lx) * cos(4π * grid2D.yC[j]/Ly) for i in 1:nx, j in 1:ny]
    ∂hy_on_h_data = @. -(4π/Ly) * [sin(2π * grid2D.xC[i]/Lx) * sin(4π * grid2D.yC[j]/Ly) for i in 1:nx, j in 1:ny]
    
    ∂u_on_u_data = @. (-6π/Lx) * [sin(6π * grid2D.xF[i]/Lx) * cos(2π * grid2D.yC[j]/Ly) for i in 1:nx, j in 1:ny]
    ∂v_on_v_data = @. (6π/Ly) * [sin(8π * grid2D.xC[i]/Lx) * cos(6π * grid2D.yF[j]/Ly) for i in 1:nx, j in 1:ny]

    h2D = Field(Centre, Center, Cdata, grid2D)
    u2D = Field(Face, Center, udata, grid2D)
    v2D = Field(Centre, Face, vdata, grid2D)

    𝐼huactual2D = Field(Face, Centre, 𝐼hudata, grid2D)
    𝐼hvactual2D = Field(Centre, Face, 𝐼hvdata, grid2D)
    𝐼uactual2D = Field(Centre, Centre, 𝐼udata, grid2D)
    𝐼vactual2D = Field(Centre, Centre, 𝐼vdata, grid2D)
    𝐼hutest2D = Field(Face, Centre, zero(𝐼hudata), grid2D)
    𝐼hvtest2D = Field(Centre, Face, zero(𝐼hvdata), grid2D)
    𝐼utest2D = Field(Centre, Centre, zero(𝐼udata), grid2D)
    𝐼vtest2D = Field(Centre, Centre, zero(𝐼vdata), grid2D)

    ∂huactual2D = Field(Face, Centre, ∂hudata, grid2D)
    ∂hvactual2D = Field(Centre, Face, ∂hvdata, grid2D)
    ∂uactual2D = Field(Centre, Centre, ∂udata, grid2D)
    ∂vactual2D = Field(Centre, Centre, ∂vdata, grid2D)
    ∂hutest2D = Field(Face, Centre, zero(∂hudata), grid2D)
    ∂hvtest2D = Field(Centre, Face, zero(∂hvdata), grid2D)
    ∂utest2D = Field(Centre, Centre, zero(∂udata), grid2D)
    ∂vtest2D = Field(Centre, Centre, zero(∂vdata), grid2D)

    ∂hx_on_h_actual2D = Field(Centre, Centre, ∂hx_on_h_data, grid2D)
    ∂hy_on_h_actual2D = Field(Centre, Centre, ∂hy_on_h_data, grid2D)
    ∂u_on_u_actual2D = Field(Face, Centre, ∂u_on_u_data, grid2D)
    ∂v_on_v_actual2D = Field(Centre, Face, ∂v_on_v_data, grid2D)

    h2D_from_outer = Field2D(Centre, Centre, Cdata[1:nx, 1:ny], grid2D)
    u2D_from_outer = Field2D(Face, Centre, udata[1:nx, 1:ny], grid2D)
    v2D_from_outer = Field2D(Centre, Face, vdata[1:nx, 1:ny], grid2D)

    @test typeof(h2D) <: Field2D{Centre, Centre}
    @test typeof(u2D) <: Field2D{Face, Centre}
    @test typeof(v2D) <: Field2D{Centre, Face}

    @test h2D.grid == grid2D
    @test u2D.grid == grid2D
    @test v2D.grid == grid2D

    @test h2D.grid == h2D_from_outer.grid
    @test u2D.grid == u2D_from_outer.grid
    @test v2D.grid == v2D_from_outer.grid
    @test h2D.data ≈ h2D_from_outer.data
    @test u2D.data ≈ u2D_from_outer.data
    @test v2D.data ≈ v2D_from_outer.data
    
    for field in [h1D, u1D, h2D, u2D, v2D]
        @test typeof(field) <: AbstractField
    end

    @test test_𝐼x(𝐼huactual2D, 𝐼hutest2D, h2D)
    @test test_𝐼y(𝐼hvactual2D, 𝐼hvtest2D, h2D)
    @test test_𝐼x(𝐼uactual2D, 𝐼utest2D, u2D)
    @test test_𝐼y(𝐼vactual2D, 𝐼vtest2D, v2D)

    @test test_∂x(∂huactual2D, ∂hutest2D, h2D)
    @test test_∂y(∂hvactual2D, ∂hvtest2D, h2D)
    @test test_∂x(∂uactual2D, ∂utest2D, u2D)
    @test test_∂y(∂vactual2D, ∂vtest2D, v2D)

    @test test_∂x(∂hx_on_h_actual2D, ∂utest2D, h2D)
    @test test_∂y(∂hy_on_h_actual2D, ∂vtest2D, h2D)
    @test test_∂x(∂u_on_u_actual2D, ∂hutest2D, u2D)
    @test test_∂y(∂v_on_v_actual2D, ∂hvtest2D, v2D)
    =#
end
