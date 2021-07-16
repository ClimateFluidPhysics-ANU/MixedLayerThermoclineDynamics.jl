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
    data :: Array
    "The grid on which the field lives."
    grid :: Grid1D
    
    Field1D(LX, data, grid) = new{LX}(data, grid)    
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

"""
"""
function interpolate!(field_output::Field1D{Face}, field_input::Field1D{Centre})
    nx = field_input.grid.nx
    
    for i in 2:nx
        field_output.data[i] = (field_input.data[i] + field_input.data[i-1]) / 2
    end
    field_output.data[1] = (field_input.data[1] + field_input.data[nx]) / 2
end

"""
"""
function interpolate!(::Field1D{Centre}, ::Field1D{Face})
    nx = field_input.grid.nx
    
    for i in 1:nx-1
        field_output.data[i] = (field_input.data[i+1] + field_input.data[i]) / 2
    end
    
    field_output.data[nx] = (field_input.data[1] + field_input.data[nx]) / 2
end

"""
"""
function ∂x!(field_output::Field1D{Face}, field_input::Field1D{Centre})
    nx, dx = field_input.grid.nx, field_input.grid.dx
    
    for i in 2:nx
        field_output.data[i] = (field_input.data[i] - field_input.data[i-1]) / dx
    end
    
    field_output.data[1] = (field_input.data[1] - field_input.data[nx]) / dx
end

"""
"""
function ∂x!(::Field1D{Centre}, ::Field1D{Face})
    nx, dx = field_input.grid.nx, field_input.grid.dx
    
    for i in 1:nx-1
        field_output.data[i] = (field_input.data[i+1] - field_input.data[i]) / dx
    end
    
    field_output.data[nx] = (field_input.data[1] - field_input.data[nx]) / dx
end

