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
struct Field1D{LX<:AbstractLocation, G} <: AbstractField
    "Array with the values of the field."
    data::AbstractArray
    "The grid on which the field lives."
    grid :: G
    
    #Field1D(LX, data, grid) = new{LX, G}(data, grid)    
end

function Field1D(LX, data::AbstractArray, grid::Grid1D)
    
    new_data = OffsetArray(zeros(grid.nx + 2*grid.hx), -grid.hx)
    field = Field1D{LX, typeof(grid)}(new_data, grid)
    fill_halos!(field, data)
    return field
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

function fill_halos!(input::Field1D{Centre, Grid1D{Periodic}}, data::AbstractArray)
    nx, hx, dx = input.grid.nx, input.grid.hx, input.grid.dx

    for i in 1:nx
        input.data[i] = data[i]
    end

    for j in 1:hx
        input.data[nx+j] = data[j]
        input.data[-j+1] = data[nx-j+1]
    end
end

function fill_halos!(input::Field1D{Face, Grid1D{Periodic}}, data::AbstractArray)
    nx, hx, dx = input.grid.nx, input.grid.hx, input.grid.dx

    for i in 1:nx
        input.data[i] = data[i]
    end

    for j in 1:hx
        input.data[nx+j] = data[j]
        input.data[-j+1] = data[nx-j+1]
    end
end
