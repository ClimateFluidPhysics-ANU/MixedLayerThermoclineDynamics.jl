module MixedLayerThermoclineDynamics

using DocStringExtensions

export Grid1D, Grid2D
export Field1D, Field2D

""" Abstract supertype for grids. """
abstract type AbstractGrid end

""" Abstract supertype for location of variables. """
abstract type AbstractLocation end

include("grids.jl")
include("fields.jl")

end # module
