using MixedLayerThermoclineDynamics, OffsetArrays, Test

const rtol_interpolation = 8e-4
const rtol_derivatives = 8e-4

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
    
    #####
    ##### 1D Fields
    #####

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
    
    Ftest_Field1D = Field(Face, zeros(nx), grid1D)
    Ctest_Field1D = Field(Centre, zeros(nx), grid1D)

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
    
    #####
    ##### 2D Fields
    #####

    CCdata = @. [sin(2π * grid2D.xC[i] / Lx) * cos(2π * grid2D.yC[j] / Ly - π/5) for i in 1:nx, j in 1:ny]
    FCdata = @. [cos(2π * grid2D.xF[i] / Lx + π/3) * cos(2π * grid2D.yC[j] / Ly) for i in 1:nx, j in 1:ny]
    CFdata = @. [sin(2π * grid2D.xC[i] / Lx + π/4) * sin(2π * grid2D.yF[j] / Ly) for i in 1:nx, j in 1:ny]

    CCdata_on_FC = @. [sin(2π * grid2D.xF[i] / Lx) * cos(2π * grid2D.yC[j] / Ly - π/5) for i in 1:nx, j in 1:ny]
    CCdata_on_CF = @. [sin(2π * grid2D.xC[i] / Lx) * cos(2π * grid2D.yF[j] / Ly - π/5) for i in 1:nx, j in 1:ny]
    FCdata_on_CC = @. [cos(2π * grid2D.xC[i] / Lx + π/3) * cos(2π * grid2D.yC[j] / Ly) for i in 1:nx, j in 1:ny]
    CFdata_on_CC = @. [sin(2π * grid2D.xC[i] / Lx + π/4) * sin(2π * grid2D.yC[j] / Ly) for i in 1:nx, j in 1:ny]

    ∂x_CCdata_on_FC = @.   2π/Lx * [cos(2π * grid2D.xF[i] / Lx) * cos(2π * grid2D.yC[j] / Ly - π/5) for i in 1:nx, j in 1:ny]
    ∂y_CCdata_on_CF = @. - 2π/Ly * [sin(2π * grid2D.xC[i] / Lx) * sin(2π * grid2D.yF[j] / Ly - π/5) for i in 1:nx, j in 1:ny]
    ∂x_FCdata_on_CC = @. - 2π/Lx * [sin(2π * grid2D.xC[i] / Lx + π/3) * cos(2π * grid2D.yC[j] / Ly) for i in 1:nx, j in 1:ny]
    ∂y_CFdata_on_CC = @.   2π/Ly * [sin(2π * grid2D.xC[i] / Lx + π/4) * cos(2π * grid2D.yC[j] / Ly) for i in 1:nx, j in 1:ny]

    ∂x_CCdata_on_CC = @.   2π/Lx * [cos(2π * grid2D.xC[i] / Lx) * cos(2π * grid2D.yC[j] / Ly - π/5) for i in 1:nx, j in 1:ny]
    ∂y_CCdata_on_CC = @. - 2π/Ly * [sin(2π * grid2D.xC[i] / Lx) * sin(2π * grid2D.yC[j] / Ly - π/5) for i in 1:nx, j in 1:ny]
    
    ∂x_FCdata_on_FC = @. - 2π/Lx * [sin(2π * grid2D.xF[i] / Lx + π/3) * cos(2π * grid2D.yC[j] / Ly) for i in 1:nx, j in 1:ny]
    ∂y_CFdata_on_CF = @.   2π/Ly * [sin(2π * grid2D.xC[i] / Lx + π/4) * cos(2π * grid2D.yF[j] / Ly) for i in 1:nx, j in 1:ny]

    CC_Field2D = Field(Centre, Center, CCdata, grid2D)
    FC_Field2D = Field(Face, Center, FCdata, grid2D)
    CF_Field2D = Field(Centre, Face, CFdata, grid2D)

    CCdata_on_FC_Field2D = Field(Face, Centre, CCdata_on_FC, grid2D)
    CCdata_on_CF_Field2D = Field(Centre, Face, CCdata_on_CF, grid2D)
    FCdata_on_CC_Field2D = Field(Centre, Centre, FCdata_on_CC, grid2D)
    CFdata_on_CC_Field2D = Field(Centre, Centre, CFdata_on_CC, grid2D)
    
    FCtest_Field2D = Field(Face, Centre, zeros(nx, ny), grid2D)
    CFtest_Field2D = Field(Centre, Face, zeros(nx, ny), grid2D)
    CCtest_Field2D = Field(Centre, Centre, zeros(nx, ny), grid2D)

    ∂x_CCdata_on_FC_Field2D = Field(Face, Centre, ∂x_CCdata_on_FC, grid2D)
    ∂y_CCdata_on_CF_Field2D = Field(Centre, Face, ∂y_CCdata_on_CF, grid2D)
    ∂x_FCdata_on_CC_Field2D = Field(Centre, Centre, ∂x_FCdata_on_CC, grid2D)
    ∂y_CFdata_on_CC_Field2D = Field(Centre, Centre, ∂y_CFdata_on_CC, grid2D)
    ∂x_CCdata_on_CC_Field2D = Field(Centre, Centre, ∂x_CCdata_on_CC, grid2D)
    ∂y_CCdata_on_CC_Field2D = Field(Centre, Centre, ∂y_CCdata_on_CC, grid2D)
    ∂x_FCdata_on_FC_Field2D = Field(Face, Centre, ∂x_FCdata_on_FC, grid2D)
    ∂y_CFdata_on_CF_Field2D = Field(Centre, Face, ∂y_CFdata_on_CF, grid2D)

    CC_Field2D_from_array = Field2D(Centre, Centre, CCdata[1:nx, 1:ny], grid2D)
    FC_Field2D_from_array = Field2D(Face, Centre, FCdata[1:nx, 1:ny], grid2D)
    CF_Field2D_from_array = Field2D(Centre, Face, CFdata[1:nx, 1:ny], grid2D)

    @test typeof(CC_Field2D) <: Field2D{Centre, Centre}
    @test typeof(FC_Field2D) <: Field2D{Face, Centre}
    @test typeof(CF_Field2D) <: Field2D{Centre, Face}

    @test CC_Field2D.grid == grid2D
    @test FC_Field2D.grid == grid2D
    @test CF_Field2D.grid == grid2D

    @test CC_Field2D.grid == CC_Field2D_from_array.grid
    @test FC_Field2D.grid == FC_Field2D_from_array.grid
    @test CF_Field2D.grid == CF_Field2D_from_array.grid
    @test CC_Field2D.data ≈ CC_Field2D_from_array.data
    @test FC_Field2D.data ≈ FC_Field2D_from_array.data
    @test CF_Field2D.data ≈ CF_Field2D_from_array.data
    
    for field in [C_Field1D, F_Field1D, CC_Field2D, FC_Field2D, CF_Field2D]
        @test typeof(field) <: AbstractField
    end

    @test test_𝐼x(CCdata_on_FC_Field2D, FCtest_Field2D, CC_Field2D)
    @test test_𝐼y(CCdata_on_CF_Field2D, CFtest_Field2D, CC_Field2D)
    @test test_𝐼x(FCdata_on_CC_Field2D, CCtest_Field2D, FC_Field2D)
    @test test_𝐼y(CFdata_on_CC_Field2D, CCtest_Field2D, CF_Field2D)

    @test test_∂x(∂x_CCdata_on_FC_Field2D, FCtest_Field2D, CC_Field2D)
    @test test_∂y(∂y_CCdata_on_CF_Field2D, CFtest_Field2D, CC_Field2D)
    @test test_∂x(∂x_FCdata_on_CC_Field2D, CCtest_Field2D, FC_Field2D)
    @test test_∂y(∂y_CFdata_on_CC_Field2D, CCtest_Field2D, CF_Field2D)
    @test test_∂x(∂x_CCdata_on_CC_Field2D, CCtest_Field2D, CC_Field2D)
    @test test_∂y(∂y_CCdata_on_CC_Field2D, CCtest_Field2D, CC_Field2D)
    @test test_∂x(∂x_FCdata_on_FC_Field2D, FCtest_Field2D, FC_Field2D)
    @test test_∂y(∂y_CFdata_on_CF_Field2D, CFtest_Field2D, CF_Field2D)
end
