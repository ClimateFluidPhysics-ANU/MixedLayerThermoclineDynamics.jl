""" Type for location at the cell centres. """
struct Centre <: AbstractLocation end 

# US spelling alias
Center = Centre

""" Type for location at the cell faces. """
struct Face <: AbstractLocation end 

"""
    struct Field{LX<:AbstractLocation, LY<:AbstractLocation}

A field datatype.

$(TYPEDFIELDS)
"""
struct Field{LX<:AbstractLocation, LY<:Union{AbstractLocation, Nothing}, G}
    "Array with the values of the field."
    data :: Array
    "The grid on which the field lives."
    grid :: G
    
    Field(LX, LY, data, grid::G) where G = new{LX, LY, G}(data, grid)    
end

"""
    Field1D(LX, data, grid::Grid1D)

Constructs a 1D field of `data` at location `LX` on `grid`.
"""
Field1D(LX, data, grid::Grid1D) = Field(LX, Nothing, data, grid)

"""
    Field2D(LX, LY, data, grid::Grid2D)

Constructs a 2D field of `data` at location `(LX, LY)` on `grid`.
"""
Field2D(LX, LY, data, grid::Grid2D) = Field(LX, LY, data, grid)
