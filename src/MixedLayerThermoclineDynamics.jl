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
  
  AbstractField,
  Field1D,
  Field2D,
  Field,
  interpolate!

""" Abstract supertype for grids. """
abstract type AbstractGrid end

""" Abstract supertype for location of variables. """
abstract type AbstractLocation end

""" Abstract supertype for fields. """
abstract type AbstractField end

include("grids.jl")
include("fields.jl")

end # module
