var documenterSearchIndex = {"docs":
[{"location":"#MixedLayerThermoclineDynamics","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics","text":"","category":"section"},{"location":"","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics","text":"Modules = [MixedLayerThermoclineDynamics]","category":"page"},{"location":"#MixedLayerThermoclineDynamics.AbstractField","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.AbstractField","text":"Abstract supertype for fields. \n\n\n\n\n\n","category":"type"},{"location":"#MixedLayerThermoclineDynamics.AbstractGrid","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.AbstractGrid","text":"Abstract supertype for grids. \n\n\n\n\n\n","category":"type"},{"location":"#MixedLayerThermoclineDynamics.AbstractLocation","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.AbstractLocation","text":"Abstract supertype for location of variables. \n\n\n\n\n\n","category":"type"},{"location":"#MixedLayerThermoclineDynamics.AbstractTopology","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.AbstractTopology","text":"Abstract supertype for topology of grids. \n\n\n\n\n\n","category":"type"},{"location":"#MixedLayerThermoclineDynamics.Bounded","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.Bounded","text":"Type for location at the cell centres. \n\n\n\n\n\n","category":"type"},{"location":"#MixedLayerThermoclineDynamics.Centre","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.Centre","text":"Type for location at the cell centres. \n\n\n\n\n\n","category":"type"},{"location":"#MixedLayerThermoclineDynamics.Face","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.Face","text":"Type for location at the cell faces. \n\n\n\n\n\n","category":"type"},{"location":"#MixedLayerThermoclineDynamics.Field1D","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.Field1D","text":"struct Field1D{LX<:AbstractLocation}\n\nA field datatype for 1D objects.\n\ndata::AbstractArray\nArray with the values of the field.\ngrid::Any\nThe grid on which the field lives.\n\n\n\n\n\n","category":"type"},{"location":"#MixedLayerThermoclineDynamics.Field2D","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.Field2D","text":"struct Field2D{LX<:AbstractLocation, LY<:AbstractLocation}\n\nA field datatype for 2D objects.\n\ndata::Array\nArray with the values of the field.\ngrid::Grid2D\nThe grid on which the field lives.\n\n\n\n\n\n","category":"type"},{"location":"#MixedLayerThermoclineDynamics.Grid1D","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.Grid1D","text":"struct Grid1D{Tx<:AbstractTopology} <: AbstractGrid\n\nReturns a one-dimensional staggered grid with topology Tx.\n\nnx::Int64\nNumber of points in x-direction\nhx::Int64\nNumber of halo points in x-direction\ndx::Float64\nGrid spacing in x-direction\nLx::Float64\nDomain extent in x-direction\nxF::AbstractArray\nCell faces in x-direction\nxC::AbstractArray\nCell centres in x-direction\n\n\n\n\n\n","category":"type"},{"location":"#MixedLayerThermoclineDynamics.Grid1D-NTuple{4, Any}","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.Grid1D","text":"Grid1D(Tx, nx, x_start, x_end; hx=0)\n\nConstruct a one-dimensional staggered grid on domain x ∈ [x_start, x_end] with topology Tx, with nx interior grid points, and hx halo points.\n\n\n\n\n\n","category":"method"},{"location":"#MixedLayerThermoclineDynamics.Grid2D","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.Grid2D","text":"struct Grid2D{Tx<:AbstractTopology, Ty<:AbstractTopology} <: AbstractGrid\n\nReturns a two-dimensional staggered grid with topologies {Tx, Ty}.\n\nnx::Int64\nNumber of points in x-direction\nny::Int64\nNumber of points in y-direction\nhx::Int64\nNumber of halo points in x-direction\nhy::Int64\nNumber of halo points in y-direction\ndx::Float64\nGrid spacing in x-direction\ndy::Float64\nGrid spacing in y-direction\nLx::Float64\nDomain extent in x-direction\nLy::Float64\nDomain extent in y-direction\nxF::AbstractArray\nCell faces in x-direction\nxC::AbstractArray\nCell centres in x-direction\nyF::AbstractArray\nCell faces in y-direction\nyC::AbstractArray\nCell centres in y-direction\n\n\n\n\n\n","category":"type"},{"location":"#MixedLayerThermoclineDynamics.Grid2D-NTuple{8, Any}","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.Grid2D","text":"Grid2D(Tx, Ty, nx, ny, x_start, x_end, y_start, y_end; hx=0, hy=0)\n\nConstruct a two-dimensional staggered grid on domain (x, y) ∈ [x_start, x_end] x [y_start, y_end] with topologies {Tx, Ty}, with {nx, ny} interior grid points, and {hx, hy} halo points.\n\n\n\n\n\n","category":"method"},{"location":"#MixedLayerThermoclineDynamics.Periodic","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.Periodic","text":"Type for location at the cell centres. \n\n\n\n\n\n","category":"type"},{"location":"#MixedLayerThermoclineDynamics.Field-Tuple{Any, Any, Any, Grid2D}","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.Field","text":"Field(LX, LY, data, grid::Grid2D)\n\nConstructs a 2D field of data at location (LX, LY) on grid.\n\n\n\n\n\n","category":"method"},{"location":"#MixedLayerThermoclineDynamics.Field-Tuple{Any, Any, Grid1D}","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.Field","text":"Field(LX, data, grid::Grid1D)\n\nConstructs a 1D field of data at location LX on grid.\n\n\n\n\n\n","category":"method"},{"location":"#MixedLayerThermoclineDynamics.construct_centres-Tuple{AbstractTopology, Any, Any, Any, Any, Any}","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.construct_centres","text":"construct_centres(T::AbstractTopology, n, h, d, L, start_center)\n\nReturns the cell centers locations for a dimension with topology T, number of grid points n and h halo points, grid spacing d, and extent L.\n\n\n\n\n\n","category":"method"},{"location":"#MixedLayerThermoclineDynamics.construct_faces-Tuple{AbstractTopology, Any, Any, Any, Any, Any}","page":"MixedLayerThermoclineDynamics","title":"MixedLayerThermoclineDynamics.construct_faces","text":"construct_faces(T::AbstractTopology, n, h, d, L)\n\nReturns the face locations for a dimension with topology T, number of grid points n and h halo points, grid spacing d, and extent L.\n\n\n\n\n\n","category":"method"}]
}
