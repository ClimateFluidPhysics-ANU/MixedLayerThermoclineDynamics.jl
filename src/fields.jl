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
    data
    "The grid on which the field lives."
    grid :: G
    
    #Field1D(LX, data, grid) = new{LX, G}(data, grid)    
end

function Field1D(LX, data, grid::Grid1D)
    
    construct_halos!(data, grid)
    return Field1D{LX, typeof(grid)}(data, grid)
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

function construct_halos!(array, grid)
    hx = grid.hx

    for i in 1:hx
        pushfirst!(array, 0.0)
        push!(array, 0.0)
    end
    OffsetArray(array, -hx)
end

function fill_halos(input::Field1D{Centre, Grid1D{Periodic}})
    nx, hx, dx = grid.nx, grid.hx, grid.dx

    for j in 1:hx
        input.data[nx+j] = input.data[j]
        input.data[-j+1] = input.data[nx-j+1]
    end
end
