"""
    struct Grid1D{Tx<:AbstractTopology} <: AbstractGrid

Returns a one-dimensional staggered `grid` with topology `Tx`.

$(TYPEDFIELDS)
"""
struct Grid1D{Tx<:AbstractTopology} <: AbstractGrid
    "Number of points in x-direction"
    nx::Int
    "Number of halo points in x-direction"
    hx::Int
    "Grid spacing in x-direction"
    dx::Float64
    "Domain extent in x-direction"
    Lx::Float64
    "Cell faces in x-direction"
    xF::AbstractArray
    "Cell centres in x-direction"
    xC::AbstractArray
end

"""
    Grid1D(Tx, nx, x_start, x_end; hx=1)

Construct a one-dimensional staggered `grid` on domain `x ∈ [x_start, x_end]` with topology
`Tx`, with `nx` interior grid points, and `hx` halo points.

Example
=======

```jldoctest
julia> using MixedLayerThermoclineDynamics

julia> grid = Grid1D(Periodic(), 10, 0, 2.0)
1-Dimensional Grid
  ├───────── topology: Periodic
  ├─ domain extent Lx: 2.0
  ├──── resolution nx: 10
  ├── grid spacing dx: 0.2
  └─── halo points nx: 1
```
"""
function Grid1D(Tx, nx, x_start, x_end; hx=1)

    if hx == 0; throw(error("Number of halo points in x cannot be zero")); end

    Lx = x_end - x_start
    
    dx = Lx/nx
    
    xF = construct_faces(Tx, nx, hx, dx, Lx, x_start)
    xC = construct_centres(Tx, nx, hx, dx, Lx, x_start)

    return Grid1D{typeof(Tx)}(nx, hx, dx, Lx, xF, xC)
end

"""
    struct Grid2D{Tx<:AbstractTopology, Ty<:AbstractTopology} <: AbstractGrid

Returns a two-dimensional staggered `grid` with topologies `{Tx, Ty}`.

$(TYPEDFIELDS)
"""
struct Grid2D{Tx<:AbstractTopology, Ty<:AbstractTopology} <: AbstractGrid
    "Number of points in x-direction"
    nx::Int
    "Number of points in y-direction"
    ny::Int
    "Number of halo points in x-direction"
    hx::Int
    "Number of halo points in y-direction"
    hy::Int
    "Grid spacing in x-direction"
    dx::Float64
    "Grid spacing in y-direction"
    dy::Float64
    "Domain extent in x-direction"
    Lx::Float64
    "Domain extent in y-direction"
    Ly::Float64
    "Cell faces in x-direction"
    xF::AbstractArray
    "Cell centres in x-direction"
    xC::AbstractArray
    "Cell faces in y-direction"
    yF::AbstractArray
    "Cell centres in y-direction"
    yC::AbstractArray
end

"""
    Grid2D(Tx, Ty, nx, ny, x_start, x_end, y_start, y_end; hx=1, hy=1)

Construct a two-dimensional staggered `grid` on domain `(x, y) ∈ [x_start, x_end] x [y_start, y_end]`
with topologies `{Tx, Ty}`, with `{nx, ny}` interior grid points, and `{hx, hy}` halo points.

Example
=======

```jldoctest
julia> using MixedLayerThermoclineDynamics

julia> grid = Grid2D(Periodic(), Periodic(), 10, 15, 0, 2.0, 0, 3.0)
2-Dimensional Grid
  ├──── topology in x: Periodic
  ├─ domain extent Lx: 2.0
  ├──── resolution nx: 10
  ├── grid spacing dx: 0.2
  ├─── halo points nx: 1
  ├──── topology in y: Periodic
  ├─ domain extent Ly: 3.0
  ├──── resolution ny: 15
  ├── grid spacing dy: 0.2
  └─── halo points ny: 1
```
"""
function Grid2D(Tx, Ty, nx, ny, x_start, x_end, y_start, y_end; hx=1, hy=1)

    if hx == 0; throw(error("Number of halo points in x cannot be zero")); end
    if hy == 0; throw(error("Number of halo points in y cannot be zero")); end

    Lx = x_end - x_start
    Ly = y_end - y_start
    
    dx, dy = Lx/nx, Ly/ny
    
    xF = construct_faces(Tx, nx, hx, dx, Lx, x_start)
    xC = construct_centres(Tx, nx, hx, dx, Lx, x_start)
    yF = construct_faces(Ty, ny, hy, dy, Ly, y_start)
    yC = construct_centres(Ty, ny, hy, dy, Ly, y_start)
    
    return Grid2D{typeof(Tx), typeof(Ty)}(nx, ny, hx, hy, dx, dy, Lx, Ly, xF, xC, yF, yC)
end

##
# helper functions

"""
    construct_faces(T::AbstractTopology, n, h, d, L)

Returns the face locations for a dimension with topology `T`, number of grid points `n` and `h`
halo points, grid spacing `d`, and extent `L`.
"""
function construct_faces(T::AbstractTopology, n, h, d, L, start)
    start_face  = start
      end_face  = isa(T, Periodic) ? start + L - d : start + L
    total_faces = isa(T, Periodic) ? n             : n + 1
    
    F = range(start_face - h*d, stop = end_face + h*d, length = total_faces + 2h)
    
    return OffsetArray(F, -h)
end

"""
    construct_centres(T::AbstractTopology, n, h, d, L, start_center)

Returns the cell centers locations for a dimension with topology `T`, number of grid points
`n` and `h` halo points, grid spacing `d`, and extent `L`.
"""
function construct_centres(T::AbstractTopology, n, h, d, L, start)
    start_center = start + d/2
      end_center = start + L - d/2
    
    C = range(start_center - h*d, stop = end_center + h*d, length = n + 2h)
    
    return OffsetArray(C, -h)
end

show(io::IO, grid::Grid1D{Tx}) where Tx =
     print(io, "1-Dimensional Grid\n",
               "  ├───────── topology: ", Tx, '\n',
               "  ├─ domain extent Lx: ", grid.Lx, '\n',
               "  ├──── resolution nx: ", grid.nx, '\n',
               "  ├── grid spacing dx: ", grid.dx, '\n',
               "  └─── halo points nx: ", grid.hx)

show(io::IO, grid::Grid2D{Tx, Ty}) where {Tx, Ty} =
     print(io, "2-Dimensional Grid\n",
               "  ├──── topology in x: ", Tx, '\n',
               "  ├─ domain extent Lx: ", grid.Lx, '\n',
               "  ├──── resolution nx: ", grid.nx, '\n',
               "  ├── grid spacing dx: ", grid.dx, '\n',
               "  ├─── halo points nx: ", grid.hx, '\n',
               "  ├──── topology in y: ", Ty, '\n',
               "  ├─ domain extent Ly: ", grid.Ly, '\n',
               "  ├──── resolution ny: ", grid.ny, '\n',
               "  ├── grid spacing dy: ", grid.dy, '\n',
               "  └─── halo points ny: ", grid.hy)
