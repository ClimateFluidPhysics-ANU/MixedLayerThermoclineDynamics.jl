module MixedLayerThermoclineDynamics

using DocStringExtensions

export Grid1D, Grid2D
export Field1D, Field2D

"Abstract supertype for grids."
abstract type AbstractGrid end

include("grids.jl")
include("fields.jl")

end # module
