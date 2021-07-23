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
    data :: OffsetArray
    "The grid on which the field lives."
    grid :: G
    
    Field1D(LX, data, grid::G) where G = new{LX, G}(data, grid)
end

function Field1D(LX, data::Array, grid::Grid1D)
    
    data_with_halos = OffsetArray(zeros(grid.nx + 2*grid.hx), -grid.hx)
    
    @. data_with_halos[1:nx] = data
    
    field = Field1D(LX, data_with_halos, grid)

    fill_halos!(field)

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
