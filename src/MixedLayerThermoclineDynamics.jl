module MixedLayerThermoclineDynamics

using DocStringExtensions

export Grid1D, Grid2D

"Abstract supertype for grids."
abstract type AbstractGrid end

include("grids.jl")

end # module
