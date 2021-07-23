module MixedLayerThermoclineDynamics

using
  DocStringExtensions,
  OffsetArrays

export
  AbstractGrid,
  AbstractTopology,
  Periodic,
  Bounded,
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
  fill_halos!

""" Abstract supertype for topology of grids. """
abstract type AbstractTopology end

""" Type for location at the cell centres. """
struct Periodic <: AbstractTopology end 

""" Type for location at the cell centres. """
struct Bounded <: AbstractTopology end 

""" Abstract supertype for grids. """
abstract type AbstractGrid end

""" Abstract supertype for location of variables. """
abstract type AbstractLocation end

""" Abstract supertype for fields. """
abstract type AbstractField end

include("grids.jl")
include("fields.jl")

end # module
