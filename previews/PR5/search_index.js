var documenterSearchIndex = {"docs":
[{"location":"#MixedLayerThermoclineDynamics","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics","text":"","category":"section"},{"location":"","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics","text":"Modules = [MixedLayerThermoclineDynamics]","category":"page"},{"location":"#MixedLayerThermoclineDynamics.AbstractGrid","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.AbstractGrid","text":"Abstract supertype for grids. \n\n\n\n\n\n","category":"type"},{"location":"#MixedLayerThermoclineDynamics.AbstractLocation","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.AbstractLocation","text":"Abstract supertype for location of variables. \n\n\n\n\n\n","category":"type"},{"location":"#MixedLayerThermoclineDynamics.Centre","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.Centre","text":"Type for location at the cell centres. \n\n\n\n\n\n","category":"type"},{"location":"#MixedLayerThermoclineDynamics.Face","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.Face","text":"Type for location at the cell faces. \n\n\n\n\n\n","category":"type"},{"location":"#MixedLayerThermoclineDynamics.Field","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.Field","text":"struct Field{LX<:AbstractLocation, LY<:AbstractLocation}\n\nA field datatype.\n\ndata::Array\nArray with the values of the field.\ngrid::Any\nThe grid on which the field lives.\n\n\n\n\n\n","category":"type"},{"location":"#MixedLayerThermoclineDynamics.Grid1D","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.Grid1D","text":"struct Grid1D <: AbstractGrid\n\nReturns a one-dimensional staggered grid.\n\nnx::Int64\nNumber of points in x-direction\ndx::Float64\nGrid spacing in x-direction\nLx::Float64\nDomain extent in x-direction\nxF::Vector{T} where T\nCell faces in x-direction\nxC::Vector{T} where T\nCell centres in x-direction\n\n\n\n\n\n","category":"type"},{"location":"#MixedLayerThermoclineDynamics.Grid2D","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.Grid2D","text":"struct Grid2D <: AbstractGrid\n\nReturns a two-dimensional staggered grid.\n\nnx::Int64\nNumber of points in x-direction\nny::Int64\nNumber of points in y-direction\ndx::Float64\nGrid spacing in x-direction\ndy::Float64\nGrid spacing in y-direction\nLx::Float64\nDomain extent in x-direction\nLy::Float64\nDomain extent in y-direction\nxF::Vector{T} where T\nCell faces in x-direction\nxC::Vector{T} where T\nCell centres in x-direction\nyF::Vector{T} where T\nCell faces in y-direction\nyC::Vector{T} where T\nCell centres in y-direction\n\n\n\n\n\n","category":"type"},{"location":"#MixedLayerThermoclineDynamics.Field1D-Tuple{Any, Any, Grid1D}","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.Field1D","text":"Field1D(LX, data, grid::Grid1D)\n\nConstructs a 1D field of data at location LX on grid.\n\n\n\n\n\n","category":"method"},{"location":"#MixedLayerThermoclineDynamics.Field2D-Tuple{Any, Any, Any, Grid2D}","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.Field2D","text":"Field2D(LX, LY, data, grid::Grid2D)\n\nConstructs a 2D field of data at location (LX, LY) on grid.\n\n\n\n\n\n","category":"method"}]
}
