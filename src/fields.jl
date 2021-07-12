include("grids.jl")

export Grid1D, Grid2D


"""
    abstract type AbstractLocation

Datatype to store location of variables. Subtypes:
1. Centre
2. Face
"""
abstract type AbstractLocation end



"""
    struct Centre <: AbstractLocation

Datatype to represent a variable on cell centre.
"""
struct Centre <: AbstractLocation end 



"""
    struct Face <: AbstractLocation

Datatype to represent a variable on cell face.
"""
struct Face <: AbstractLocation end



"""
    struct Field1D{Locx<:AbstractLocation}

Constructs a field datatype for a 1D variable.

$(TYPEDFIELDS)
"""
struct Field1D{Locx<:AbstractLocation}
	"Array storing the values of a 1D variable"
    data::Array
    "Grid properties for 1D variable"
    grid::Grid1D
end



"""
    struct Field2D{Locx<:AbstractLocation, Locy<:AbstractLocation}

Constructs a field datatype for a 2D variable.

$(TYPEDFIELDS)
"""
struct Field2D{Locx<:AbstractLocation, Locy<:AbstractLocation}
    
	"Array storing the values of a 2D variable"
    data::Array
    "Grid properties for 2D variable"
    grid::Grid2D
end