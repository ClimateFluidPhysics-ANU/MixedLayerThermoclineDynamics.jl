""" Type for location at the cell centres. """
struct Centre <: AbstractLocation end 

# US spelling alias
Center = Centre

""" Type for location at the cell faces. """
struct Face <: AbstractLocation end 

"""
    struct Field1D{LX<:AbstractLocation}

A field datatype for 1D objects.

$(TYPEDFIELDS)
"""
struct Field1D{LX<:AbstractLocation} <: AbstractField
    "Array with the values of the field."
    data
    "The grid on which the field lives."
    grid :: Grid1D   
end

function Field1D(LX, data, grid)

    return Field1D{LX}(fill_halos(data, grid, LX), grid)
end

"""
    struct Field2D{LX<:AbstractLocation, LY<:AbstractLocation}

A field datatype for 2D objects.

$(TYPEDFIELDS)
"""
struct Field2D{LX<:AbstractLocation, LY<:AbstractLocation} <: AbstractField
    "Array with the values of the field."
    data :: Array
    "The grid on which the field lives."
    grid :: Grid2D
    
    Field2D(LX, LY, data, grid) = new{LX, LY}(data, grid)    
end

#function Field2D(LX, LY, data, grid)
#    C = fill_halos(data, grid)
#    return Field2D{LX, LY}(C, grid)
#end

"""
    Field(LX, data, grid::Grid1D)

Constructs a 1D field of `data` at location `LX` on `grid`.
"""
Field(LX, data, grid::Grid1D) = Field1D(LX, data, grid)

"""
    Field(LX, LY, data, grid::Grid2D)

Constructs a 2D field of `data` at location `(LX, LY)` on `grid`.
"""
Field(LX, LY, data, grid::Grid2D) = Field2D(LX, LY, data, grid)

function fill_halos(data, grid::Grid1D{Periodic}, LX::Type{Centre})
    nx, hx, dx = grid.nx, grid.hx, grid.dx
    offsetData = zero(grid.xC)
    
    for j in 1:hx
        offsetData[nx+j] = data[j]
        offsetData[-j+1] = data[nx-j+1]
    end

    for i in 1:nx
        offsetData[i] = data[i]
    end
    return offsetData
end

function fill_halos(data, grid::Grid1D{Periodic}, LX::Type{Face})
    nx, hx, dx = grid.nx, grid.hx, grid.dx
    offsetData = zero(grid.xF)

    for j in 1:hx
        offsetData[nx+j] = data[j]
        offsetData[-j+1] = data[nx-j+1]
    end

    for i in 1:nx
        offsetData[i] = data[i]
    end
    return offsetData
end

############################################################################
#-------------------------------Interpolation------------------------------#

𝐼xᶠ(i, f::Field1D{Face})   = (f.data[i] + f.data[i+1]) / 2
𝐼xᶜ(i, f::Field1D{Centre}) = (f.data[i] + f.data[i-1]) / 2

𝐼xᶠ(i, j, f::Field2D{Face, Centre})   = (f.data[i, j] + f.data[i+1, j]) / 2
𝐼xᶜ(i, j, f::Field2D{Centre, Centre}) = (f.data[i, j] + f.data[i-1, j]) / 2
𝐼yᶠ(i, j, f::Field2D{Centre, Face})   = (f.data[i, j+1] + f.data[i, j]) / 2
𝐼yᶜ(i, j, f::Field2D{Centre, Centre}) = (f.data[i, j] + f.data[i, j-1]) / 2

"""
    𝐼x!(output::Field1D{Centre}, input::Field1D{Face})

Interpolates a 1D field of `data` from Face to Centre.
"""
function 𝐼x!(output::Field1D{Centre}, input::Field1D{Face})
    nx = input.grid.nx
        
    for i in 1:nx
        output.data[i] = 𝐼xᶠ(i, input.data)
    end
end

"""
    𝐼x!(output::Field1D{Face}, input::Field1D{Centre})

Interpolates a 1D field of `data` from Centre to Face.
"""
function 𝐼x!(output::Field1D{Face}, input::Field1D{Centre})
    nx = input.grid.nx
        
    for i in 1:nx+1
        output.data[i] = 𝐼xᶜ(i, input.data)
    end
end

"""
    𝐼x!(output::Field2D{Centre}, input::Field2D{Face})

Interpolates a 2D field of `data` from T to U grid.
"""
function 𝐼x!(output::Field2D{Face, Centre}, input::Field2D{Centre, Centre})
    nx, ny = input.grid.nx, input.grid.ny
    
    for i in 1:nx+1
        for j in 1:ny
            output.data[i, j] = 𝐼xᶠ(i, j, input.data)
        end
    end
end

"""
    𝐼x!(output::Field2D{Face}, input::Field2D{Centre})

Interpolates a 2D field of `data` from U to T grid.
"""
function 𝐼x!(output::Field2D{Centre, Centre}, input::Field2D{Face, Centre})
    nx, ny = input.grid.nx, input.grid.ny
        
    for i in 1:nx
        for j in 1:ny
            output.data[i, j] = 𝐼xᶜ(i, j, input.data)
        end
    end
end

"""
    𝐼y!(output::Field1D{Face}, input::Field2D{Centre})

Interpolates a 2D field of `data` from T to V grid.
"""
function 𝐼y!(output::Field2D{Centre, Face}, input::Field2D{Centre, Centre})
    nx, ny = input.grid.nx, input.grid.ny
        
    for i in 1:nx
        for j in 1:ny+1
            output.data[i, j] = 𝐼yᶠ(i, j, input.data)
        end
    end
end

"""
    𝐼y!(output::Field2D{Face}, input::Field2D{Centre})

Interpolates a 2D field of `data` from V to T grid.
"""
function 𝐼y!(output::Field2D{Centre, Centre}, input::Field2D{Centre, Face})
    nx, ny = input.grid.nx, input.grid.ny
        
    for i in 1:nx
        for j in 1:ny
            output.data[i, j] = 𝐼yᶜ(i, j, input.data)
        end
    end
end

############################################################################
#--------------------------------Derivatives-------------------------------#

δxᶠ(i, f::Field1D{Face})   = f.data[i+1] - f.data[i]
δxᶜ(i, f::Field1D{Centre}) = f.data[i] - f.data[i-1]

δxᶠ(i, j, f::Field2D{Face, Centre})   = f.data[i+1, j] - f.data[i, j]
δxᶜ(i, j, f::Field2D{Centre, Centre}) = f.data[i, j] - f.data[i-1, j]
δyᶠ(i, j, f::Field2D{Centre, Face})   = f.data[i, j+1] - f.data[i, j]
δyᶜ(i, j, f::Field2D{Centre, Centre}) = f.data[i, j] - f.data[i, j-1]

"""
    ∂x!(output::Field1D{Centre}, input::Field1D{Face})

Interpolates a 1D field of `data` from Face to Centre grid.
"""
function ∂x!(output::Field1D{Centre}, input::Field1D{Face})
    nx = input.grid.nx
    dx = input.grid.dx
        
    for i in 1:nx
        output.data[i] = δxᶠ(i, input.data)/dx
    end
end

"""
    ∂x!(output::Field1D{Face}, input::Field1D{Center})

Interpolates a 1D field of `data` from Centre to Face grid.
"""
function ∂x!(output::Field1D{Face}, input::Field1D{Centre})
    nx = input.grid.nx
    dx = input.grid.dx

    for i in 1:nx+1
        output.data[i] = δxᶜ(i, input.data)/dx
    end
end

"""
    ∂x!(output::Field2D{Centre, Centre}, input::Field2D{Face, Centre})

Interpolates a 1D field of `data` from U to T grid.
"""
function ∂x!(output::Field2D{Centre, Centre}, input::Field2D{Face, Centre})
    nx, ny = input.grid.nx, input.grid.ny
    dx = input.grid.dx
        
    for i in 1:nx
        for j in 1:ny
            output.data[i, j] = δxᶠ(i, j, input.data)/dx
        end
    end
end

"""
    ∂x!(output::Field2D{Face, Centre}, input::Field2D{Centre, Centre})

Interpolates a 2D field of `data` from T to U grid.
"""
function ∂x!(output::Field2D{Face, Centre}, input::Field2D{Centre, Centre})
    nx, ny = input.grid.nx, input.grid.ny
    dx = input.grid.dx

    for i in 1:nx+1
        for j in 1:ny
            output.data[i, j] = δxᶜ(i, j, input.data)/dx
        end
    end
end

"""
    ∂y!(output::Field2D{Centre, Centre}, input::Field2D{Centre, Face})

Interpolates a 2D field of `data` from V to T grid.
"""
function ∂y!(output::Field2D{Centre, Centre}, input::Field2D{Centre, Face})
    nx, ny = input.grid.nx, input.grid.ny
    dy = input.grid.dy
        
    for i in 1:nx
        for j in 1:ny
            output.data[i, j] = δyᶠ(i, j, input.data)/dy
        end
    end
end

"""
    ∂y!(output::Field2D{Centre, Face}, input::Field2D{Centre, Centre})

Interpolates a 2D field of `data` from T to V grid.
"""
function ∂y!(output::Field2D{Centre, Face}, input::Field2D{Centre, Centre})
    nx, ny = input.grid.nx, input.grid.ny
    dy = input.grid.dy

    for i in 1:nx
        for j in 1:ny+1
            output.data[i, j] = δyᶜ(i, j, input.data)/dy
        end
    end
end