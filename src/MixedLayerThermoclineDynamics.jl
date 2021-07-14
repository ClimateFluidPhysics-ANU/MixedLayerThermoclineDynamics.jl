module MixedLayerThermoclineDynamics

using DocStringExtensions

export
  AbstractGrid,
  Grid1D,
  Grid2D,
  
  AbstractLocation,
  Face,
  Center,
  Centre,
  
  Field,
  Field1D,
  Field2D

""" Abstract supertype for grids. """
abstract type AbstractGrid end

""" Abstract supertype for location of variables. """
abstract type AbstractLocation end

include("grids.jl")
include("fields.jl")

end # module
