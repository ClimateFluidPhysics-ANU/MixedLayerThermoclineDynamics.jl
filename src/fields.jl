""" Type for location at the cell centres. """
struct Centre <: AbstractLocation end 

# US spelling alias
Center = Centre

""" Type for location at the cell faces. """
struct Face <: AbstractLocation end 

"""
    struct Field1D{LX<:AbstractLocation, G} <: AbstractField

A field datatype for 1D objects.

$(TYPEDFIELDS)
"""
struct Field1D{LX<:AbstractLocation, G} <: AbstractField
    "Array with the values of the field."
    data :: OffsetArray
    "The grid on which the field lives."
    grid :: G
    
    Field1D(LX, data, grid::G) where G = new{LX, G}(data, grid)
end

function Field1D(LX, data::Array, grid::Grid1D)
    nx, hx = grid.nx, grid.hx
    
    data_with_halos = OffsetArray(zeros(grid.nx + 2*grid.hx), -grid.hx)
    
    @. data_with_halos[1:nx] = data
    
    field = Field1D(LX, data_with_halos, grid)

    fill_halos!(field)

    return field
end

"""
    struct Field2D{LX<:AbstractLocation, LY<:AbstractLocation, G} <: AbstractField

A field datatype for 2D objects.

$(TYPEDFIELDS)
"""
struct Field2D{LX<:AbstractLocation, LY<:AbstractLocation, G} <: AbstractField
    "Array with the values of the field."
    data :: OffsetArray
    "The grid on which the field lives."
    grid :: G
    
    Field2D(LX, LY, data, grid::G) where G = new{LX, LY, G}(data, grid)    
end

function Field2D(LX, LY, data::Array, grid::Grid2D)
    nx, hx = grid.nx, grid.hx
    ny, hy = grid.ny, grid.hy
    
    data_with_halos = OffsetArray(zeros(nx + 2hx, ny + 2hy), -hx, -hy)
    
    @. data_with_halos[1:nx, 1:ny] = data
    
    field = Field2D(LX, LY, data_with_halos, grid)

    fill_halos!(field)

    return field
end

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

#####
##### Intepolations
#####

𝐼xᶜ(i, f::Field1D{Face})   = (f.data[i] + f.data[i+1]) / 2
𝐼xᶠ(i, f::Field1D{Centre}) = (f.data[i-1] + f.data[i]) / 2

𝐼xᶜᶜ(i, j, f::Field2D{Face, Centre})   = (f.data[i, j] + f.data[i+1, j]) / 2
𝐼xᶠᶜ(i, j, f::Field2D{Centre, Centre}) = (f.data[i-1, j] + f.data[i, j]) / 2
𝐼yᶜᶜ(i, j, f::Field2D{Centre, Face})   = (f.data[i, j] + f.data[i, j+1]) / 2
𝐼yᶜᶠ(i, j, f::Field2D{Centre, Centre}) = (f.data[i, j-1] + f.data[i, j]) / 2

"""
    𝐼x!(output::Field1D{<:Any}, input::Field1D{<:Any, Grid1D{Periodic}})

Interpolates a 1D `input` field to the location where the `output` field lives for 1D grids
with periodic boundary conditions.
"""
function 𝐼x!(output::Field1D{Centre}, input::Field1D{Face, Grid1D{Periodic}})
    nx = input.grid.nx
    
    for i in 1:nx
        output.data[i] = 𝐼xᶜ(i, input)
    end
    
    fill_halos!(output)
    
    return nothing
end

function 𝐼x!(output::Field1D{Face}, input::Field1D{Centre, Grid1D{Periodic}})
    nx = input.grid.nx
    
    for i in 1:nx
        output.data[i] = 𝐼xᶠ(i, input)
    end
    
    fill_halos!(output)
    
    return nothing
end

function 𝐼x!(output::Field1D{LX}, input::Field1D{LX, Grid1D{Periodic}}) where LX
    @. output.data = input.data
    
    fill_halos!(output)

    return nothing
end

"""
    𝐼x!(output::Field2D{<:Any, <:Any}, input::Field2D{<:Any, <:Any, Grid2D{Periodic, Periodic}})

Interpolates a 2D `input` field to the location where the `output` field lives for 2D grids
with periodic boundary conditions.
"""
function 𝐼x!(output::Field2D{Face, Centre}, input::Field2D{Centre, Centre, Grid2D{Periodic, Periodic}})
    nx, ny = input.grid.nx, input.grid.ny
    
    for j in 1:ny
        output.data[1, j] = (input.data[1, j] + input.data[nx, j]) / 2
    end

    for j in 1:ny, i in 1:nx
        output.data[i, j] = 𝐼xᶠᶜ(i, j, input)
    end

    fill_halos!(output)
    
    return nothing
end

function 𝐼x!(output::Field2D{Centre, Centre}, input::Field2D{Face, Centre, Grid2D{Periodic, Periodic}})
    nx, ny = input.grid.nx, input.grid.ny
    
    for j in 1:ny
        output.data[nx, j] = (input.data[1, j] + input.data[nx, j]) / 2
    end

    for j in 1:ny, i in 1:nx
        output.data[i, j] = 𝐼xᶜᶜ(i, j, input)
    end

    fill_halos!(output)
    
    return nothing
end

function 𝐼x!(output::Field2D{LX, LY}, input::Field2D{LX, LY, Grid2D{Periodic, Periodic}}) where {LX<:AbstractLocation, LY<:AbstractLocation}
    @. output.data = input.data
    
    fill_halos!(output)

    return nothing
end

"""
    𝐼y!(output::Field2D{<:Any, <:Any}, input::Field2D{<:Any, <:Any, Grid1D{Periodic, Periodic}})

Interpolates a 2D `input` field to the location where the `output` field lives for 2D grids
with periodic boundary conditions.
"""
function 𝐼y!(output::Field2D{Centre, Face}, input::Field2D{Centre, Centre, Grid2D{Periodic, Periodic}})
    nx, ny = input.grid.nx, input.grid.ny
        
    for i in 1:nx
        output.data[i, 1] = (input.data[i, 1] + input.data[i, ny]) / 2
    end
    
    for j in 1:ny, i in 1:nx
        output.data[i, j] = 𝐼yᶜᶠ(i, j, input)
    end

    fill_halos!(output)
    
    return nothing
end

function 𝐼y!(output::Field2D{Centre, Centre}, input::Field2D{Centre, Face, Grid2D{Periodic, Periodic}})
    nx, ny = input.grid.nx, input.grid.ny
        
    for i in 1:nx
        output.data[i, ny] = (input.data[i, 1] + input.data[i, ny]) / 2
    end
    
    for j in 1:ny, i in 1:nx
        output.data[i, j] = 𝐼yᶜᶜ(i, j, input)
    end

    fill_halos!(output)
    
    return nothing
end

function 𝐼y!(output::Field2D{LX, LY}, input::Field2D{LX, LY, Grid2D{Periodic, Periodic}}) where {LX<:AbstractLocation, LY<:AbstractLocation}
    @. output.data = input.data
    
    fill_halos!(output)

    return nothing
end

#####
##### Derivatives
#####

δxᶜ(i, f::Field1D{Face})   = f.data[i+1] - f.data[i]
δxᶠ(i, f::Field1D{Centre}) = f.data[i] - f.data[i-1]

δxᶜᶜ(i, j, f::Field2D{Face, Centre})   = f.data[i+1, j] - f.data[i, j]
δxᶠᶜ(i, j, f::Field2D{Centre, Centre}) = f.data[i, j] - f.data[i-1, j]
δyᶜᶜ(i, j, f::Field2D{Centre, Face})   = f.data[i, j+1] - f.data[i, j]
δyᶜᶠ(i, j, f::Field2D{Centre, Centre}) = f.data[i, j] - f.data[i, j-1]

"""
    ∂x!(output::Field1D{Centre}, input::Field1D{Face, Grid1D{Periodic}})

Compute the derivative of a 1D `input` field onto the location where the `output` field lives
for 1D grids with periodic boundary conditions.

Compute the derivative of a 1D field of `data` from face to centre in x-direction.
"""
function ∂x!(output::Field1D{Centre}, input::Field1D{Face, Grid1D{Periodic}})
    nx, dx = input.grid.nx, input.grid.dx
    
    for i in 1:nx
        output.data[i] = δxᶜ(i, input) / dx
    end
    
    fill_halos!(output)
    
    return nothing
end

"""
    ∂x!(output::Field1D{Face}, input::Field1D{Centre, Grid1D{Periodic}})

Compute the derivative of a 1D `input` field onto the location where the `output` field lives
for 1D grids with periodic boundary conditions.

Compute the derivative of a 1D field of `data` from centre to face in x-direction.
"""
function ∂x!(output::Field1D{Face}, input::Field1D{Centre, Grid1D{Periodic}})
    nx, dx = input.grid.nx, input.grid.dx
    
    for i in 1:nx
        output.data[i] = δxᶠ(i, input) / dx
    end
    
    fill_halos!(output)
    
    return nothing
end

"""
    ∂x!(output::Field1D{Centre}, input::Field1D{Centre, Grid1D{Periodic}})

Compute the derivative of a 1D `input` field onto the location where the `output` field lives
for 1D grids with periodic boundary conditions.

Compute the derivative of a 1D field of `data` from centre to centre in x-direction.
"""
function ∂x!(output::Field1D{Centre}, input::Field1D{Centre, Grid1D{Periodic}})
    nx, dx = input.grid.nx, input.grid.dx
    
    for i in 1:nx
        output.data[i] = δxᶠ(i, input) / dx
    end
    
    fill_halos!(output)
    
    𝐼x!(output, input)
    
    fill_halos!(output)
    
    return nothing
end

"""
    ∂x!(output::Field1D{Face}, input::Field1D{Face, Grid1D{Periodic}})

Compute the derivative of a 1D `input` field onto the location where the `output` field lives
for 1D grids with periodic boundary conditions.

Compute the derivative of a 1D field of `data` from face to face in x-direction.
"""
function ∂x!(output::Field1D{Face}, input::Field1D{Face, Grid1D{Periodic}})
    nx, dx = input.grid.nx, input.grid.dx
    
    for i in 1:nx
        output.data[i] = δxᶜ(i, input) / dx
    end
    
    fill_halos!(output)
    
    𝐼x!(output, input)
    
    fill_halos!(output)
    
    return nothing
end

"""
    ∂x!(output::Field2D{Centre, Centre}, input::Field2D{Face, Centre, Grid2D{Periodic, Periodic}})

Compute the derivative of a 2D `input` field onto the location where the `output` field lives
for 2D grids with periodic boundary conditions.

Compute the derivative of a 2D field of `data` from face to centre in x-direction.
"""
function ∂x!(output::Field2D{Centre, Centre}, input::Field2D{Face, Centre, Grid2D{Periodic, Periodic}})
    nx, ny = input.grid.nx, input.grid.ny
    dx = input.grid.dx
        
    for j in 1:ny
        output.data[nx, j] = (input.data[1, j] - input.data[nx, j]) / dx
    end

    for j in 1:ny, i = 1:nx
        output.data[i, j] = δxᶜᶜ(i, j, input) / dx
    end

    fill_halos!(output)
    
    return nothing
end

"""
    ∂x!(output::Field2D{face, Centre}, input::Field2D{Face, Centre, Grid2D{Periodic, Periodic}})

Compute the derivative of a 2D `input` field onto the location where the `output` field lives
for 2D grids with periodic boundary conditions.

Compute the derivative of a 2D field of `data` from face to face in x-direction.
"""
function ∂x!(output::Field2D{Face, Centre}, input::Field2D{Face, Centre, Grid2D{Periodic, Periodic}})
    nx, ny = input.grid.nx, input.grid.ny
    dx = input.grid.dx
        
    for j in 1:ny
        output.data[nx, j] = (input.data[1, j] - input.data[nx, j]) / dx
    end

    for j in 1:ny, i = 1:nx
        output.data[i, j] = δxᶜᶜ(i, j, input) / dx
    end

    fill_halos!(output)

    𝐼x!(output, input)

    fill_halos!(output)

    return nothing
end

"""
    ∂x!(output::Field2D{Face, Centre}, input::Field2D{Centre, Centre, Grid2D{Periodic, Periodic}})

Compute the derivative of a 2D `input` field onto the location where the `output` field lives
for 2D grids with periodic boundary conditions.

Compute the derivative of a 2D field of `data` from centre to face in x-direction.
"""
function ∂x!(output::Field2D{Face, Centre}, input::Field2D{Centre, Centre, Grid2D{Periodic, Periodic}})
    nx, ny = input.grid.nx, input.grid.ny
    dx = input.grid.dx

    for j in 1:ny
        output.data[1, j] = (input.data[1, j] - input.data[nx, j]) / dx
    end

    for j in 1:ny, i in 1:nx
        output.data[i, j] = δxᶠᶜ(i, j, input)/dx
    end
    
    fill_halos!(output)

    return nothing
end

"""
    ∂x!(output::Field2D{Centre, Centre}, input::Field2D{Face, Centre, Grid2D{Periodic, Periodic}})

Compute the derivative of a 2D `input` field onto the location where the `output` field lives
for 2D grids with periodic boundary conditions.

Compute the derivative of a 2D field of `data` from centre to centre in x-direction.
"""
function ∂x!(output::Field2D{Centre, Centre}, input::Field2D{Centre, Centre, Grid2D{Periodic, Periodic}})
    nx, ny = input.grid.nx, input.grid.ny
    dx = input.grid.dx

    for j in 1:ny
        output.data[1, j] = (input.data[1, j] - input.data[nx, j]) / dx
    end

    for j in 1:ny, i in 1:nx
        output.data[i, j] = δxᶠᶜ(i, j, input)/dx
    end
    
    fill_halos!(output)

    𝐼x!(output, input)

    fill_halos!(output)

    return nothing
end

"""
    ∂y!(output::Field2D{Centre, Centre}, input::Field2D{Centre, Face, Grid2D{Periodic, Periodic}})

Compute the derivative of a 2D `input` field onto the location where the `output` field lives
for 2D grids with periodic boundary conditions.

Compute the derivative of a 2D field of `data` from face to centre in y-direction.
"""
function ∂y!(output::Field2D{Centre, Centre}, input::Field2D{Centre, Face, Grid2D{Periodic, Periodic}})
    nx, ny = input.grid.nx, input.grid.ny
    dy = input.grid.dy
        
    for i in 1:nx
        output.data[i, ny] = (input.data[i, 1] - input.data[i, ny]) / dy
    end

    for j in 1:ny, i in 1:nx
        output.data[i, j] = δyᶜᶜ(i, j, input) / dy
    end

    fill_halos!(output)
    
    return nothing
end

"""
    ∂y!(output::Field2D{Centre, Face}, input::Field2D{Centre, Face, Grid2D{Periodic, Periodic}})

Compute the derivative of a 2D `input` field onto the location where the `output` field lives
for 2D grids with periodic boundary conditions.

Compute the derivative of a 2D field of `data` from face to face in y-direction.
"""
function ∂y!(output::Field2D{Centre, Face}, input::Field2D{Centre, Face, Grid2D{Periodic, Periodic}})
    nx, ny = input.grid.nx, input.grid.ny
    dy = input.grid.dy
        
    for i in 1:nx
        output.data[i, ny] = (input.data[i, 1] - input.data[i, ny]) / dy
    end

    for j in 1:ny, i in 1:nx
        output.data[i, j] = δyᶜᶜ(i, j, input) / dy
    end

    fill_halos!(output)

    𝐼y!(output, input)

    fill_halos!(output)
    
    return nothing
end

"""
    ∂y!(output::Field2D{Centre, Face}, input::Field2D{Centre, Centre, Grid2D{Periodic, Periodic}})

Compute the derivative of a 2D `input` field onto the location where the `output` field lives
for 2D grids with periodic boundary conditions.

Compute the derivative of a 2D field of `data` from centre to face in y-direction.
"""
function ∂y!(output::Field2D{Centre, Face}, input::Field2D{Centre, Centre, Grid2D{Periodic, Periodic}})
    nx, ny = input.grid.nx, input.grid.ny
    dy = input.grid.dy

    for i in 1:nx
        output.data[i, 1] = (input.data[i, 1] - input.data[i, ny]) / dy
    end

    for j in 1:ny, i in 1:nx
        output.data[i, j] = δyᶜᶠ(i, j, input)/dy
    end

    fill_halos!(output)
    
    return nothing
end

"""
    ∂y!(output::Field2D{Centre, Centre}, input::Field2D{Centre, Centre, Grid2D{Periodic, Periodic}})

Compute the derivative of a 2D `input` field onto the location where the `output` field lives
for 2D grids with periodic boundary conditions.

Compute the derivative of a 2D field of `data` from centre to centre in y-direction.
"""
function ∂y!(output::Field2D{Centre, Centre}, input::Field2D{Centre, Centre, Grid2D{Periodic, Periodic}})
    nx, ny = input.grid.nx, input.grid.ny
    dy = input.grid.dy

    for i in 1:nx
        output.data[i, 1] = (input.data[i, 1] - input.data[i, ny]) / dy
    end

    for j in 1:ny, i in 1:nx
        output.data[i, j] = δyᶜᶠ(i, j, input)/dy
    end

    fill_halos!(output)

    𝐼y!(output, input)

    fill_halos!(output)
    
    return nothing
end

"""
    fill_halos!(field::Field1D{<:Any, Grid1D{Periodic}})

Fill halos for a 1D `field` that lives on a grid with periodic boundary conditions.
"""
function fill_halos!(field::Field1D{<:Any, Grid1D{Periodic}})
    nx, hx = field.grid.nx, field.grid.hx

    for i in 1:hx
        field.data[nx+i] = field.data[i]
        field.data[-i+1] = field.data[nx-i+1]
    end

    return nothing
 end

"""
    fill_halos!(field::Field2D{<:Any, <:Any, Grid2D{Periodic, Periodic}})

Fill halos for a 2D `field` that lives on a grid with periodic boundary conditions in both directions.
"""
function fill_halos!(field::Field2D{<:Any, <:Any, Grid2D{Periodic, Periodic}})
    nx, hx = field.grid.nx, field.grid.hx
    ny, hy = field.grid.ny, field.grid.hy

    for j in 1:ny, i in 1:hx
        field.data[nx+i, j] = field.data[i, j]
        field.data[-i+1, j] = field.data[nx-i+1, j]
    end

    for j in 1:hy, i in 1:nx
        field.data[i, ny+j] = field.data[i, j]
        field.data[i, -j+1] = field.data[i, ny-j+1]
    end
    
    return nothing
 end
